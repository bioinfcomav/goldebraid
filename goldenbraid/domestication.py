# Copyright 2013 Diego Orzaez, Univ.Politecnica Valencia, Consejo Superior de
# Investigaciones Cientificas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import division
import re
from itertools import izip_longest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from goldenbraid.views.feature_views import get_prefix_and_suffix_index
from goldenbraid.settings import (DOMESTICATION_DEFAULT_MELTING_TEMP,
                                  DOMESTICATION_MIN_OLIGO_LENGTH,
                                  ENZYMES_USED_IN_GOLDENBRAID, PUPD_PREFIX,
                                  OLIGO_UNIVERSAL, DOMESTICATED_SEQ,
                                  MINIMUN_PCR_LENGTH, CRYSPER_SEQ)
from goldenbraid.models import Feature, Count
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from goldenbraid.tags import TARGET_MONOCOT, TARGET_DICOT
from goldenbraid.utils import get_ret_sites, has_rec_sites


def get_codontable():
    'get codontable'
    bases = ['T', 'C', 'A', 'G']
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDE'
    amino_acids += 'EGGGG'
    codon_table = dict(zip(codons, amino_acids))
    return codon_table


def domesticate_for_synthesis(seqrec, category, prefix, suffix,
                              with_intron=False):
    kind = category
    seq = seqrec.seq
    if not with_intron:
        seq = seq.upper()
    new_seq = _remove_rec_sites(seq)[0]
    seqs_for_sintesis, prefix, suffix = _add_tags_to_pcrproducts([new_seq],
                                                                 prefix,
                                                                 suffix,
                                                                 kind)

    try:
        count = Count.objects.get(name=DOMESTICATED_SEQ)
    except Count.DoesNotExist:
        count = Count.objects.create(name=DOMESTICATED_SEQ, value=1)
    next_value = count.next
    seq_name = DOMESTICATED_SEQ + '_' + next_value
    vector_seq = _get_stripped_vector_seq()
    prepared_new_seq = prefix + new_seq + suffix + vector_seq
    seq_for_synthesis = str(seqs_for_sintesis[0])
    prepared_seq = SeqRecord(prepared_new_seq, name=seq_name, id=seq_name)

    return seq_for_synthesis, prepared_seq


def domesticate(seqrec, category, prefix, suffix, with_intron=False):
    kind = category
    seq = seqrec.seq
    if not with_intron:
        seq = seq.upper()
    min_melting_temp = DOMESTICATION_DEFAULT_MELTING_TEMP
    new_seq, rec_site_pairs, fragments = _remove_rec_sites(seq)
    segments = _get_pcr_segments(new_seq, rec_site_pairs, fragments)

    pcr_products = [str(new_seq[s['start']:s['end'] + 1]) for s in segments]
    oligos = _get_oligos(new_seq, segments, min_melting_temp)
    oligos = _add_tags_to_oligos(oligos, prefix, suffix, kind)
    # coprobar que los overhangs son distintos posiciones 12-15
    forw_bin_sites = []
    for oligo in oligos:
        olig_for = oligo[0]
        for_bin = olig_for[11:15]
        if for_bin in forw_bin_sites:
            raise RuntimeError('Repeated overhang')
        forw_bin_sites.append(for_bin)

    # print oligos
    pcr_products, prefix, suffix = _add_tags_to_pcrproducts(pcr_products,
                                                            prefix, suffix,
                                                            kind)

    oligo_pcrs = []
    for pcr, oligo in zip(pcr_products, oligos):
        oligo_pcrs.append({'pcr_product': pcr, 'oligo_forward': oligo[0],
                          'oligo_reverse': oligo[1]})

    try:
        count = Count.objects.get(name=DOMESTICATED_SEQ)
    except Count.DoesNotExist:
        count = Count.objects.create(name=DOMESTICATED_SEQ, value=1)
    next_value = count.next

    vector_seq = _get_stripped_vector_seq()
    prepared_new_seq = prefix + new_seq + suffix + vector_seq
    seq_name = DOMESTICATED_SEQ + '_' + next_value
    new_seq_record = SeqRecord(prepared_new_seq, name=seq_name, id=seq_name)

    if with_intron:
        cds = _get_cds_from_seq(seq, prefix)
        new_seq_record.features.append(cds)

    return oligo_pcrs, new_seq_record


def _get_cds_from_seq(seq, prefix):

    cds_locs = []
    for match in re.finditer('[A-Z]+', str(seq)):
        cds_locs.append(FeatureLocation(match.start() + len(prefix),
                                        match.end() + len(prefix), strand=1))

    qualifiers = {'translation': Seq(_get_upper_nucls(seq)).translate()}
    cds = SeqFeature(CompoundLocation(cds_locs), type='CDS', strand=1,
                     qualifiers=qualifiers)
    return cds


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
    segments = [{'start': s[0], 'end': s[1]} for s in segments]

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


def is_dna_palindrome(seq):
    nucl_palindromes = [('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')]
    seq_pals = []
    if divmod(len(seq), 2)[1] != 0:
        return False
    for num in range(0, len(seq)):
        rev_num = len(seq) - 1 - num
        nucl_a = seq[num]
        nucl_b = seq[rev_num]
        if (nucl_a, nucl_b) in nucl_palindromes:
            seq_pals.append(True)
        else:
            seq_pals.append(False)
    return all(seq_pals)


def _get_segments_from_rec_site(frag_5, frag_3, rec_site, prev_seq_len,
                                overhangs):
    change_pos = 0
    for letter1, letter2 in zip(rec_site['original'], rec_site['modified']):
        if letter1 != letter2:
            break
        change_pos += 1
    change_index = prev_seq_len + len(frag_5) + change_pos
    fow_end = change_index + 1
    rev_start = fow_end - 3

    overhang = get_overhang(rev_start, fow_end, prev_seq_len, frag_5, frag_3,
                            rec_site)
    if is_dna_palindrome(overhang):
        fow_end = change_index + 2
        rev_start = fow_end - 3
        overhang = get_overhang(rev_start, fow_end, prev_seq_len, frag_5,
                                frag_3, rec_site)
    count = 0
    while overhang in overhangs:
        rev_start += 1
        fow_end += 1
        overhang = get_overhang(rev_start, fow_end, prev_seq_len, frag_5,
                                frag_3, rec_site)
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
    if frag_3 is not None:
        index = 0
        while len(overhang) < 4:
            overhang += frag_3[index].upper()
    return overhang


def _get_stripped_vector_seq():
    pupd = Feature.objects.get(uniquename='pUPD')
    vec_seq = pupd.residues
    pre_suf_size = get_prefix_and_suffix_index(vec_seq, pupd.enzyme_in[0])
    prefix_index, suffix_index, prefix_size = pre_suf_size
    prefix_start = prefix_index
    suffix_end = suffix_index + prefix_size
    if prefix_start > suffix_end:
        stripped_seq = vec_seq[prefix_start:]
        stripped_seq += vec_seq[:suffix_end]
    else:
        stripped_seq = vec_seq[prefix_start:suffix_end]
    return stripped_seq


def _guess_prefix_suffix_tag(kind, prefix, suffix):
    '''It select the needed prefix and suffix to add  to oligos and
    pcr_products for especial category cases'''

    if kind in ('CDS (B3-B4-B5)', 'CDS (B3-B4)'):
        prefix = 'A'
    elif kind == 'SP (B3)':
        prefix = 'A'
        suffix = 'GCAGCC'
    elif kind == 'NTAG (B2)':
        prefix = 'CC'
        suffix = 'TCAATG'
    elif kind == 'CTAG (B5)':
        prefix = 'GCAGGG'
    return prefix, suffix


def _add_tags_to_pcrproducts(pcr_products, prefix, suffix, kind):
    pcr_products_with_tags = []
    prefix, suffix = _guess_prefix_suffix_tag(kind, prefix, suffix)
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
    return pcr_products_with_tags, prefix, suffix


def _add_tags_to_oligos(oligos, prefix, suffix, kind):
    oligos_with_tags = []
    prefix, suffix = _guess_prefix_suffix_tag(kind, prefix, suffix)

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
    rec_sites_regex = re.compile(rec_sites_regex, flags=re.IGNORECASE)
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
    for fragment, rec_site_in_seq in izip_longest(fragments, rec_sites_in_seq):
        new_seq += fragment
        if rec_site_in_seq is not None:
            _cumulative_patch += fragment + rec_site_in_seq
            new_rec_site = _domesticate_rec_site(rec_site_in_seq,
                                                 _cumulative_patch,
                                                 rec_sites_regex)
            rec_site_pairs.append({'original': rec_site_in_seq,
                                   'modified': new_rec_site})

            new_seq += new_rec_site
    coding_seq = Seq(_get_upper_nucls(seq))
    new_coding_seq = Seq(_get_upper_nucls(new_seq))
    if str(coding_seq.translate()) != str(new_coding_seq.translate()):
        msg = 'The generated sequence does not produce the same peptide'
        raise ValueError(msg)
    if rec_sites_regex.search(str(new_seq)):
        msg = 'Not all rec_sites modified'
        raise ValueError(msg)
    return new_seq, rec_site_pairs, fragments


def _get_upper_nucls(seq):
    return ''.join([nucl for nucl in seq if nucl.isupper()])


def change_nucl_in_intron_rec_site(rec_site, rec_sites_regex):
    for index, nucl in enumerate(rec_site):
        if nucl.islower():
            for new_nucl in ('a', 't', 'c', 'g'):
                new_rec_site = rec_site[:index] + new_nucl
                if index < (len(rec_site) - 1):
                    new_rec_site += rec_site[index + 1:]

                if not rec_sites_regex.match(new_rec_site):
                    return new_rec_site


def _domesticate_rec_site(rec_site, patch, rec_sites_regex):
    '''it converts a rec site in a disabled rec_site. It changes one nucleotide
    but tries not to change aa.
    It can not convert in an already unusable rec_site'''
    with_intron = False
    for letter in rec_site:
        if letter.islower():
            with_intron = True
    if with_intron:
        return change_nucl_in_intron_rec_site(rec_site, rec_sites_regex)
    # get a dictionary for codon_table
    codon_table = get_codontable()
    # get the last complete codon
    lastcodon = ''
    baseindex_to_change = ''
    coding_patch = _get_upper_nucls(patch)
    frame = divmod(len(coding_patch), 3)[1] + 1
    if frame == 1:
        baseindex_to_change = -1
        lastcodon = coding_patch[-3:]
    elif frame == 2:
        baseindex_to_change = -2
        lastcodon = coding_patch[-4:-1]
    elif frame == 3:
        baseindex_to_change = -3
        lastcodon = coding_patch[-5:-2]
    else:
        raise ValueError()

    # if lastcodon is Metionine, change the previous codon,
    # since Met does not have alternative codon
    if lastcodon == 'ATG':
        if rec_site == 'GCGATG':
            baseindex_to_change = -4
            lastcodon = coding_patch[-6:-3]

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


# cryspers
def domestication_crysper(seq, category=None, prefix=None, suffix=None):
    if len(seq) != 20:
        raise ValueError('Seq length different 20')
    if category == TARGET_DICOT and str(seq[0]) != 'G':
        raise ValueError('First nucleotide must be G for target dicot category')
    if category == TARGET_MONOCOT and str(seq[0]) != 'A':
        raise ValueError('First nucleotide must be G for target monocot category')

    if has_rec_sites(seq):
        msg = 'This secuence can not be domesticated. It has internal restriction sites'
        raise ValueError(msg)

    if category:
        prefix = prefix[:3]

    try:
        count = Count.objects.get(name=CRYSPER_SEQ)
    except Count.DoesNotExist:
        count = Count.objects.create(name=CRYSPER_SEQ, value=1)
    next_value = count.next

    prepared_seq = Seq(prefix + seq + suffix)
    seq_name = CRYSPER_SEQ + '_' + next_value
    new_seq_record = SeqRecord(prepared_seq, name=seq_name, id=seq_name)
    return new_seq_record
