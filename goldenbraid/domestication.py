'''
Created on 2013 ots 7

@author: peio
'''
from __future__ import  division
import re
from itertools import izip_longest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from goldenbraid.views.feature_views import (parse_rebase_file,
                                             get_prefix_and_suffix_index)
from goldenbraid.settings import (REBASE_FILE,
                                  DOMESTICATION_DEFAULT_MELTING_TEMP,
                                  DOMESTICATION_MIN_OLIGO_LENGTH,
                                  ENZYMES_USED_IN_GOLDENBRAID, PUPD_PREFIX,
                                  OLIGO_UNIVERSAL, DOMESTICATED_SEQ,
                                  MINIMUN_PCR_LENGTH)
from goldenbraid.models import Feature, Count


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
    # coprobar que los overhangs son distintos posiciones 12-15
    forw_bin_sites = []
    for olig_for, olig_rev in oligos:
        for_bin = olig_for[11:15]
        if for_bin in forw_bin_sites:
            raise RuntimeError('Repeated overhang')
        forw_bin_sites.append(for_bin)

    # print oligos
    pcr_products = _add_tags_to_pcrproducts(pcr_products, prefix, suffix, kind)

    vector_seq = _get_stripped_vector_seq()
    prepared_new_seq = prefix + new_seq + suffix + vector_seq

    oligo_pcrs = []
    for pcr, oligo in zip(pcr_products, oligos):
        oligo_pcrs.append({'pcr_product': pcr, 'oligo_forward': oligo[0],
                          'oligo_reverse': oligo[1]})

    try:
        count = Count.objects.get(name=DOMESTICATED_SEQ)
    except Count.DoesNotExist:
        count = Count.objects.create(name=DOMESTICATED_SEQ, value=1)
    next_value = count.next
    seq_name = DOMESTICATED_SEQ + '_' + next_value
    return oligo_pcrs, SeqRecord(prepared_new_seq, name=seq_name, id=seq_name)


def _get_oligos(seq, segments, min_melting_temp):
    oligos = []
    for segment in segments:
        forward_min = segment.get('forward_min', None)
        if forward_min:
            forward_min = forward_min - segment['start'] + 1
        forw_oligo = _get_oligo(seq[segment['start']:], min_melting_temp,
                                forward_min)

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
    overhangs = []
    frag_rec_sites = zip(fragments, rec_sites)
    for index, frag_5_rec_site in enumerate(frag_rec_sites):
        frag_5 = frag_5_rec_site[0]
        rec_site = frag_5_rec_site[1]
        try:
            frag_3 = fragments[index + 1]
        except IndexError:
            frag_3 = None
        start, end, overhangs = _get_segments_from_rec_site(frag_5, frag_3,
                                                           rec_site,
                                                           acumulated_seq_len,
                                                           overhangs)
        segments['starts'].append(start)
        segments['ends'].append(end)
        acumulated_seq_len += len(frag_5) + len(rec_site['modified'])
    segments['ends'].append(len(seq) - 1)
    segments = zip(segments['starts'], segments['ends'])
    return _join_segments(segments)


def _join_segments(segments, min_length=MINIMUN_PCR_LENGTH):
    # join short segments
    segments = [{'start':s[0], 'end':s[1]} for s in segments]

    while not _all_segments_ok(segments, min_length):
        segments = _join_short_segments(segments, min_length)
    return segments


def _all_segments_ok(segments, min_length):
    for segment in segments:
        if segment['end'] - segment['start'] < min_length:
            return False
    return True


def _join_short_segments(segments, min_length):
    len_segments = len(segments)
    joined_segments = []
    skip_segment = False
    for index, segment in enumerate(segments):
        start = segment['start']
        end = segment['end']
        if end - start < min_length:
            # if not last segment
            if index + 1 != len_segments:
                if skip_segment:
                    joined_segments[-1]['end'] = segments[index + 1]['end']
                    joined_segments[-1]['forward_min'] = end + 8
                    new_segment = None
                else:
                    new_segment = {'start': start,
                                   'end': segments[index + 1]['end'],
                                   'forward_min': end + 8}
            else:
                joined_segments[-1]['end'] = end
                joined_segments[-1]['reverse_min'] = start - 8
                new_segment = None
            skip_segment = True

        else:
            if skip_segment:
                skip_segment = False
                continue
            new_segment = {'start': start, 'end': end}
        if new_segment:
            joined_segments.append(new_segment)

    return joined_segments


def  _get_segments_from_rec_site(frag_5, frag_3, rec_site, prev_seq_len,
                                 overhangs):
    change_pos = 0
    for letter1, letter2 in zip(rec_site['original'], rec_site['modified']):
        if letter1 != letter2:
            break
        change_pos += 1
    change_index = prev_seq_len + len(frag_5) + change_pos
    fow_end = change_index + 1
    rev_start = fow_end - 3
    overhang = get_overhang(rev_start, fow_end, prev_seq_len, frag_5, frag_3, rec_site)
    count = 0
    while overhang in overhangs:
        rev_start += 1
        fow_end += 1
        overhang = get_overhang(rev_start, fow_end, prev_seq_len, frag_5, frag_3, rec_site)
        if count > 10:
            msg = 'Impossible to domesticate this  sequence\n:'
            msg += 'Domesticated rec site nucleotide is too far from oligo'
            msg += ' start'
            raise RuntimeError(msg)
        count += 1

    overhangs.append(overhang)
    return rev_start, fow_end, overhangs


def get_overhang(rev_start, fow_end, prev_seq_len, frag_5, frag_3, rec_site):
    overhang_start = rev_start - prev_seq_len - len(frag_5)
    overhang_end = fow_end - prev_seq_len - len(frag_5)
    overhang = rec_site['modified'][overhang_start:overhang_end + 1]
    # si es el ultimo no pasa por aqui
    if frag_3 is not  None:
        index = 0
        while len(overhang) < 4:
            overhang += frag_3[index].upper()
    return overhang


def _get_stripped_vector_seq():
    pupd = Feature.objects.get(uniquename='pUPD')
    vec_seq = pupd.residues
    prefix_index, suffix_index, prefix_size = get_prefix_and_suffix_index(vec_seq,
                                                        pupd.enzyme_in[0])
    prefix_start = prefix_index
    suffix_end = suffix_index + prefix_size
    if prefix_start > suffix_end:
        stripped_seq = vec_seq[prefix_start:]
        stripped_seq += vec_seq[:suffix_end]
    else:
        stripped_seq = vec_seq[prefix_start:suffix_end]
    return stripped_seq


def _add_tags_to_pcrproducts(pcr_products, prefix, suffix, kind):
    pcr_products_with_tags = []
    if kind is None:
        pass
    elif kind in ('13-14-15-16 (CDS)', '13-14-15 (CDS)'):
        prefix = 'A'
    elif kind in ('13 (SP)'):
        prefix = 'A'
        suffix = 'GCAGCC'
    elif kind == '12 (NT)':
        prefix = 'CC'
    elif kind == '16 (CT)':
        prefix = 'GCAGGG'
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
    elif kind == '16 (CT)':
        prefix = 'GCAGGG'
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
    if not min_length or min_length < DOMESTICATION_MIN_OLIGO_LENGTH:
        min_length = DOMESTICATION_MIN_OLIGO_LENGTH
    oligo = []
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
    return 64.9 + 41 * (seq.count('G') + seq.count('C') - 16.4) / len_seq


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
    _cumulative_patch = ''  # it is only used to know the frame
    rec_site_pairs = []
    for fragment, rec_site_in_seq  in izip_longest(fragments,
                                                   rec_sites_in_seq):
        new_seq += fragment
        if rec_site_in_seq is not None:
            _cumulative_patch += fragment + rec_site_in_seq
            new_rec_site = _domesticate_rec_site(rec_site_in_seq,
                                                 _cumulative_patch)
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


def _domesticate_rec_site(rec_site, patch):
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
        return newsite

    # if we reach this is because no allowed domesticated site has been found
    raise ValueError('No domestication possible for ORF site ' + rec_site)
