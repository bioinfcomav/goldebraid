'''
Created on 2013 ots 7

@author: peio
'''
from __future__ import  division
import re
from itertools import izip_longest
from Bio.Alphabet import generic_dna
from math import log10

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from goldenbraid.views.feature_views import (parse_rebase_file,
                                             get_prefix_and_suffix_index)
from goldenbraid.settings import (REBASE_FILE,
                                  DOMESTICATION_DEFAULT_MELTING_TEMP)
from goldenbraid.settings import DB
from goldenbraid.models import Feature


OLIGO_UNIVERSAL = 'GCGCCGTCTCG'
PUPD_PREFIX = 'CTCG'
MINIMUN_PCR_LENGTH = 50

CATEGORIES = OrderedDict()
CATEGORIES['01-02-03-11-12 (PROM+UTR+ATG)'] = ('PROM+UTR+ATG', 'GGAG', 'AATG')
CATEGORIES['01-02-03-11 (PROM+UTR)'] = ('PROM+UTR', 'GGAG', 'CCAT')
CATEGORIES['01-02 (OP)'] = ('OP', 'GGAG', 'TCCC')
CATEGORIES['03-11-12 (MinPROM)'] = ('MinPROM', 'TCCC', 'AATG')
CATEGORIES['01 (PROM)'] = ('PROM', 'GGAG', 'TGAC')
CATEGORIES['02 (OP)'] = ('OP', 'TGAC', 'TCCC')
CATEGORIES['01-02-03-11-12B (INTERACTION ADAPTOR'] = \
                                        ('InteractionADAPTOR', 'GGAG', 'AATG')
CATEGORIES['01-02-03-11-C (PROM+UTR+mir173)'] = \
                                        ('PROM+UTR+mir173', 'GGAG', 'CCAT')
CATEGORIES['12 (NT)'] = ('NT', 'CCAT', 'AATG')
CATEGORIES['13-14-15-16 (CDS)'] = ('CDS', 'AATG', 'GCTT')
CATEGORIES['13 (SP)'] = ('SP', 'AATG', 'AGCC')
CATEGORIES['14-15-16 (CDS)'] = ('CDS', 'AGCC', 'GCTT')
CATEGORIES['13-14-15 (CDS)'] = ('CDS', 'AATG', 'GCAG')
CATEGORIES['16 (CT)'] = ('CT', 'GCAG', 'GCTT')
CATEGORIES["12-13B (5'FS)"] = ('5FS', 'CCAT', 'GTGA')
CATEGORIES['14B-15B (Target)'] = ('Target', 'GTGA', 'TCTC')
CATEGORIES["16B (3'FS)"] = ('3FS', 'TCTC', 'GCTT')
CATEGORIES['12-13 (GOI)'] = ('goi', 'CCAT', 'GTGA')
CATEGORIES['14-15(INT)'] = ('int', 'GTGA', 'TCTC')
CATEGORIES['16 (IOG)'] = ('iog', 'TCTC', 'GCTT')
CATEGORIES['12-13-14-15-16 (GOI)'] = ('goi', 'CCAT', 'GCTT')
CATEGORIES['17-21 (TER)'] = ('TER', 'GCTT', 'CGCT')

ENZYMES_USED_IN_GOLDENBRAID = ('BsmBI', 'BsaI', 'BtgZI')


def get_ret_sites(enzymes):
    'It returns the restriction site of the given enzymes'
    sites = []
    existing_enzymes = parse_rebase_file(REBASE_FILE)
    for enzyme in enzymes:
        site = existing_enzymes[enzyme]
        if '^' in site:
            raise ValueError('We only use type IIS enzymes')
        site = site.split('(')[0]
        rev_site = str(Seq(site).reverse_complement())
        sites.extend([site, rev_site])
    return sites


def get_codontable():
    'get codontable'
    bases = ['T', 'C', 'A', 'G']
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDE'
    amino_acids += 'EGGGG'
    codon_table = dict(zip(codons, amino_acids))
    return codon_table


def domesticate(seqrec, category, prefix, suffix):
    kind = category
    seq = seqrec.seq
    min_melting_temp = DOMESTICATION_DEFAULT_MELTING_TEMP
    new_seq, rec_site_pairs, fragments = _remove_rec_sites(seq)
    segments = _get_pcr_segments(new_seq, rec_site_pairs, fragments)
    pcr_products = [str(new_seq[s['start']:s['end'] + 1]) for s in segments]
    oligos = _get_oligos(new_seq, segments, min_melting_temp)
    oligos = _add_tags_to_oligos(oligos, prefix, suffix, kind)
    pcr_products = _add_tags_to_pcrproducts(pcr_products, prefix, suffix, kind)

    vector_seq = _get_stripped_vector_seq()
    prepared_new_seq = prefix + new_seq + suffix + vector_seq

    oligo_pcrs = []
    for pcr, oligo in zip(pcr_products, oligos):
        oligo_pcrs.append({'pcr_product': pcr, 'oligo_forward': oligo[0],
                          'oligo_reverse': oligo[1]})
    return oligo_pcrs, SeqRecord(prepared_new_seq, name='domesticated_seq',
                                 id='domesticated_seq')


def _get_oligos(seq, segments, min_melting_temp):
    oligos = []
    for segment in segments:
        forw_oligo = _get_oligo(seq[segment['start']:], min_melting_temp,
                                segment.get('forward_min', None))

        reverse_min = segment.get('reverse_min', None)
        if reverse_min:
            reverse_min = segment['end'] - reverse_min
        rev_oligo = _get_oligo(seq[:segment['end'] + 1].reverse_complement(),
                               min_melting_temp, reverse_min)
        oligos.append((forw_oligo, rev_oligo))
    return oligos


def _get_pcr_segments(seq, rec_sites, fragments):
    segments = {'starts': [], 'ends': []}
    segments['starts'].append(0)

    acumulated_seq_len = 0
    for frag_5, rec_site in zip(fragments, rec_sites):

        start, end = _get_segments_from_rec_site(frag_5, rec_site,
                                                 acumulated_seq_len)
        segments['starts'].append(start)
        segments['ends'].append(end)
        acumulated_seq_len += len(frag_5) + len(rec_site['modified'])
    segments['ends'].append(len(seq))
    segments = zip(segments['starts'], segments['ends'])
    return _join_short_segments(segments)


def _join_short_segments(segments, min_length=MINIMUN_PCR_LENGTH):
    # join short segments
    len_segments = len(segments)
    joined_segments = []
    skip_segment = False
    for index, segment in enumerate(segments):

        if segment[1] - segment[0] < min_length:
            # if not last segment
            if index + 1 != len_segments:
                if skip_segment:
                    joined_segments[-1]['end'] = segments[index + 1][1]
                    joined_segments[-1]['forward_min'] = segment[1]
                    new_segment = None
                else:
                    new_segment = {'start': segment[0],
                                   'end': segments[index + 1][1],
                                   'forward_min': segment[1]}
            else:
                joined_segments[-1]['end'] = segment[1]
                joined_segments[-1]['reverse_min'] = segment[0]
                new_segment = None
            skip_segment = True

        else:
            if skip_segment:
                skip_segment = False
                continue
            new_segment = {'start': segment[0], 'end': segment[1]}
        if new_segment:
            joined_segments.append(new_segment)

    return joined_segments


def  _get_segments_from_rec_site(frag_5, rec_site, prev_seq_len):
    change_pos = 0
    for letter1, letter2 in zip(rec_site['original'], rec_site['modified']):
        if letter1 != letter2:
            break
        change_pos += 1
    change_index = prev_seq_len + len(frag_5) + change_pos
    fow_end = change_index + 1
    rev_start = fow_end - 3

    return rev_start, fow_end


def _get_stripped_vector_seq():
    pupd = Feature.objects.using(DB).get(uniquename='pUPD')
    vec_seq = pupd.residues
    prefix_start, suffix_end = get_prefix_and_suffix_index(vec_seq,
                                                        pupd.enzyme_in[0])[:2]
    if prefix_start > suffix_end:
        stripped_seq = vec_seq[prefix_start:]
        stripped_seq += vec_seq[:suffix_end + 1]
    else:
        stripped_seq = vec_seq[prefix_start:suffix_end + 1]
    return stripped_seq


def _add_tags_to_pcrproducts(pcr_products, prefix, suffix, kind):
    pcr_products_with_tags = []
    if kind in ('13-14-15-16 (CDS)', '13 (SP)', '13-14-15 (CDS)'):
        prefix = 'A'
    elif kind == '12 (NT)':
        prefix = 'CC'

    len_pcr = len(pcr_products)
    for index, pcr_product in enumerate(pcr_products):
        pcr_tag = OLIGO_UNIVERSAL
        if index == 0:
            pcr_tag += PUPD_PREFIX + prefix
        pcr_tag += pcr_product
        if index + 1 == len_pcr:
            pcr_tag += suffix + str(Seq(PUPD_PREFIX).reverse_complement())

        pcr_tag += str(Seq(OLIGO_UNIVERSAL).reverse_complement())
        pcr_products_with_tags.append(pcr_tag.upper())
    return pcr_products_with_tags


def _add_tags_to_oligos(oligos, prefix, suffix, kind):
    oligos_with_tags = []
    if kind in ('13-14-15-16 (CDS)', '13 (SP)', '13-14-15 (CDS)'):
        prefix = 'A'
    elif kind == '12 (NT)':
        prefix = 'CC'
    suffix = str(Seq(suffix).reverse_complement())
    len_oligos = len(oligos)
    for index, oligo_pair in enumerate(oligos):
        oligo_tag5 = OLIGO_UNIVERSAL
        if index == 0:
            oligo_tag5 += PUPD_PREFIX + prefix
        oligo_tag5 += oligo_pair[0]
        oligo_tag3 = OLIGO_UNIVERSAL
        if index + 1 == len_oligos:
            oligo_tag3 += PUPD_PREFIX + suffix
        oligo_tag3 += oligo_pair[1]
        oligos_with_tags.append((oligo_tag5.upper(), oligo_tag3.upper()))
    return oligos_with_tags


def _get_oligo(seq, min_melting_temp, min_length=None):
    'Giving a seq and a melting temperature it return the longest oligo'
    if not min_length:
        min_length = 20
    for index in range(min_length, len(seq)):
        oligo = seq[:index]
        if _calculate_annealing_temp(oligo) >= min_melting_temp:
            break
    return str(oligo)


def _calculate_annealing_temp(seq):
    # from  http://www.basic.northwestern.edu/biotools/oligocalc.html
    # Tm (C)= 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
    seq = seq.upper()
    len_seq = len(seq)
    return 64.9 + 41 * (seq.count('G') + seq.count('C')-16.4) / len_seq


def _remove_rec_sites(seq):
    '''It modifies all rec sites in the sequence to be able to use with
    goldenbraid pipeline'''
    rec_sites = get_ret_sites(ENZYMES_USED_IN_GOLDENBRAID)
    # regex with the sites to domesticate
    rec_sites_regex = '(' + '|'.join(rec_sites) + ')'
    rec_sites_regex = re.compile(rec_sites_regex)
    rec_sites_in_seq = []
    fragments = []
    for splitted_part in rec_sites_regex.split(str(seq)):
        if rec_sites_regex.match(splitted_part):
            rec_sites_in_seq.append(splitted_part)
        else:
            fragments.append(splitted_part)

    new_seq = Seq('', alphabet=generic_dna)
    # we can not convert a rec site in another rec site
    unusable_rec_sites = rec_sites
    _cumulative_patch = ''  # it is only used to know the frame
    rec_site_pairs = []
    for fragment, rec_site_in_seq  in izip_longest(fragments,
                                                   rec_sites_in_seq):
        new_seq += fragment
        if rec_site_in_seq is not None:
            _cumulative_patch += fragment + rec_site_in_seq
            new_rec_site = _domesticate_rec_site(rec_site_in_seq,
                                                 _cumulative_patch,
                                                 unusable_rec_sites)
            unusable_rec_sites.append(new_rec_site)
            rec_site_pairs.append({'original': rec_site_in_seq,
                                   'modified': new_rec_site})

            new_seq += new_rec_site
    if str(seq.translate()) != str(new_seq.translate()):
        msg = 'The generated sequence does not produce the same peptide'
        raise ValueError(msg)
    if  rec_sites_regex.search(str(new_seq)):
        msg = 'Not all rec_sites modified'
        raise ValueError(msg)

    return new_seq, rec_site_pairs, fragments


def _domesticate_rec_site(rec_site, patch, unusable_rec_sites):
    '''it converts a rec site in a disabled rec_site. It changes one nucleotide
    but tries not to change aa.
    It can not convert in an already unusable rec_site'''
    # get a dictionary for codon_table
    codon_table = get_codontable()

    # get the last complete codon
    lastcodon = ''
    baseindex_to_change = ''
    frame = divmod(len(patch), 3)[1] + 1
    if frame == 1:
        baseindex_to_change = -1
        lastcodon = patch[-3:]
    elif frame == 2:
        baseindex_to_change = -2
        lastcodon = patch[-4:-1]
    elif frame == 3:
        baseindex_to_change = -3
        lastcodon = patch[-5:-2]
    else:
        raise ValueError()

    # if lastcodon is Metionine, change the previous codon,
    # since Met does not have alternative codon
    if lastcodon == 'ATG':
        if rec_site == 'GCGATG':
            baseindex_to_change = -4
            lastcodon = patch[-6:-3]

    # get alternative codons (same aminoacid) for the lastcodon
    alt_codons = []
    for codon in codon_table.keys():
        if codon == lastcodon:
            continue
        if codon_table.get(codon) == codon_table.get(lastcodon):
            if codon[0] == lastcodon[0]:
                if codon[1] == lastcodon[1]:
                    alt_codons.append(codon)
    # select one alternative codon
    newsite = ''
    for alt_codon in alt_codons:
        newbase = alt_codon[2]
        newsite += rec_site[0:baseindex_to_change]
        newsite += newbase
        if baseindex_to_change < -1:
            newsite += rec_site[(baseindex_to_change + 1):]
        # check that new site is not one of the already domesticated sites
        if newsite in unusable_rec_sites:
            newsite = ''
        else:
            return newsite

    # if we reach this is because no allowed domesticated site has been found
    raise ValueError('No domestication possible for ORF site ' + rec_site)
