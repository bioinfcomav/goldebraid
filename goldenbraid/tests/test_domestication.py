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

from django.test import TestCase
import os.path
import re

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import goldenbraid
from goldenbraid.domestication import (domesticate, _join_segments,
                                       _get_pcr_segments,
                                       _get_segments_from_rec_site,
                                       domesticate_for_synthesis,
                                       get_ret_sites,
                                       change_nucl_in_intron_rec_site,
                                       is_dna_palindrome,
                                       domestication_crispr)
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD
from goldenbraid.settings import MANDATORY_DOMEST_ENZYMES

TEST_DATA = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                         'goldenbraid', 'tests', 'data')


class DomesticationTest(TestCase):
    fixtures = ['auth.json'] + FIXTURES_TO_LOAD
    multi_db = True

    def test_domestication(self):
        seq = 'aggctgactatgtcagctagctgacgatcgatgctagctagctgactagctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq))
        category = 'CDS (B3-B4-B5)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[0] == {'pcr_product': 'GCGCCGTCTCGCTCGAAGGCTGACTATGTCAGCTAGCTGACGATCGATGCTAGCTAGCTGACTAGCTAGCAGGTGCTAGGAGACTGCGAGACGGCGC',
                                 'oligo_reverse': 'GCGCCGTCTCGCAGTCTCCTAGCACCTGCTA',
                                 'oligo_forward': 'GCGCCGTCTCGCTCGAAGGCTGACTATGTCAGCTAG'}

        assert oligo_pcrs[1] == {'pcr_product': 'GCGCCGTCTCGACTGGGTCATGCTAGCTTCAGCTAGCTGATCGATCGACTAGCTGATCGATCTGATCGATGCTAGCTAGCTGTACGGAGACTGGCGAGACGGCGC',
                                 'oligo_reverse': 'GCGCCGTCTCGCCAGTCTCCGTACAGCTAGC',
                                 'oligo_forward': 'GCGCCGTCTCGACTGGGTCATGCTAGCTTCA'}

        domesticated_seq = domesticate(seq, category, 'CCAT', 'AATG')[1]
        assert ('TCGCATGCTCCCGGCCGCCATGGCGGCCGCGGGAATTCGATGGGCGATGAGTGGTCTCG') in domesticated_seq.seq
        assert ('TCGCATGCTCCCGGCCGCCATGGCGGCCGCGGGAATTCGATGGGCGATGAGTGGTCTCGT') not in domesticated_seq.seq

        # prefix-suffix mod in tags
        seq = SeqIO.read(os.path.join(TEST_DATA, 'gfp.txt'), 'fasta')

        new_seq = domesticate(seq, 'CTAG (B5)', 'GCAG', 'GCTT')[1]
        result_seq = SeqIO.read(os.path.join(TEST_DATA, 'CTresult.gb'), 'gb')
        assert str(new_seq.seq).upper() == str(result_seq.seq).upper()

    def test_domestication_short_segments(self):
        seq = 'aggctgactatgtcagctaGAGACCgctgacgatcgatgctagctagctgactagctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq))
        category = 'CDS (B3-B4-B5)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG', with_intron=False)[0]

        assert oligo_pcrs[0]['pcr_product'] == 'GCGCCGTCTCGCTCGAAGGCTGACTATGTCAGCTAGAGATCGCTGACGATCGATGCTAGCTAGCTGACTAGCTAGCAGGTGCTAGGAGACTGCGAGACGGCGC'
        assert oligo_pcrs[0]['oligo_reverse'] == 'GCGCCGTCTCGCAGTCTCCTAGCACCTGCTA'
        assert oligo_pcrs[0]['oligo_forward'] == 'GCGCCGTCTCGCTCGAAGGCTGACTATGTCAGCTAGAGATCGCTGACGAT'

        seq = 'aggctgactatgtcagctagctgacgatcgatgctagctagctgactatagaggaaacccgtaacgctacgtacgCGTCTCgctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq.upper()))
        category = 'CDS (B3-B4-B5)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[1] == {'pcr_product': 'GCGCCGTCTCGCTTGCTAGCAGGTGCTAGGAGACTGGGTCATGCTAGCTTCAGCTAGCTGATCGATCGACTAGCTGATCGATCTGATCGATGCTAGCTAGCTGTACGGAGACTGGCGAGACGGCGC',
                                 'oligo_reverse': 'GCGCCGTCTCGCCAGTCTCCGTACAGCTAGC',
                                 'oligo_forward': 'GCGCCGTCTCGCTTGCTAGCAGGTGCTAGGAGACTGGGTCATGC'}

        seq = 'aggctgactatgtcagctagctgacgatcgatgctagctagctgactatagaggaaacccgtaacgctacgtacggctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgctgcgtcgggcgtctgatgctagctagctCGTCTCgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq.upper()))
        category = 'CDS (B3-B4-B5)'
        
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[2] == {'pcr_product': 'GCGCCGTCTCGCTGGGTCATGCTAGCTGATCTGATCGCTGCGTCGGGCGTCTGATGCTAGCTAGCTCGTTTCGTACGTCATCTTTTCAGTCGATCTAAATGTGAGCGAGACGGCGC',
                                 'oligo_reverse': 'GCGCCGTCTCGCTCACATTTAGATCGACTGAAAAGATGACGTACGAAACGAGCTAG',
                                 'oligo_forward': 'GCGCCGTCTCGCTGGGTCATGCTAGCTGATC'}

        gb_path = os.path.join(TEST_DATA, 'GB_DOMEST_15.gb')
        seq = SeqIO.read(gb_path, 'gb')
        category = 'OP (A2)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[1]['oligo_reverse'] == 'GCGCCGTCTCGCTCACATTCGGGACCACTCATCGCCCAT'
        assert oligo_pcrs[1]['oligo_forward'] == 'GCGCCGTCTCGACTACTCATCGCCCATCACT'

        gb_path = os.path.join(TEST_DATA, 'CHS.gb')
        seq = SeqIO.read(gb_path, 'gb')
        category = 'CDS (B3-B4-B5)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[0]['oligo_reverse'] == 'GCGCCGTCTCGCTCACATTTCAAGCACCCACACTGTGAAGCACAACAGTCTCAACA'
        assert oligo_pcrs[0]['oligo_forward'] == 'GCGCCGTCTCGCTCGAATGGTGACCGTCGAAGAAGT'

    def test_domestication_with_intron(self):
        seq = SeqIO.read(os.path.join(TEST_DATA, 'seq_with_intron.fasta'),
                         'fasta')
        oligo_pcrs = domesticate(seq, None, 'CCAT', 'AATG', with_intron=True)
        assert len(oligo_pcrs[0]) == 4

    def test_domestication_multiple_sites(self):
        seq = SeqIO.read(os.path.join(TEST_DATA, 'seq_with_rep_rec_sites.fasta'), 'fasta')
        category = 'CDS (B3-B4-B5)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')

    def test_get_segments_from_rec_site(self):
        frag5 = 'TATCGATCGATCGATGCTAGCTGATCGATCGAATCTACTACTACTACTAC'
        frag3 = 'TTATATATAT'
        rec_site = {'original': 'CGTCTC', 'modified': 'CGTATC'}
        prev_seq_len = 0
        end, start, overhangs = _get_segments_from_rec_site(frag5, frag3,
                                                            rec_site,
                                                            prev_seq_len,
                                                            overhangs=[])
        assert (51, 54) == (end, start)

    def test_domestication_palindrome(self):
        frag5 = 'AAGCTTGATACTGCAGTTCTTCAGTTCCAAAATCTTAAGCTATCACAAAAGCTAGAGGCTCAGCAGGTTGAGTGTTCTATTCTTGAGGATAAACTCTCTCAGATCAAGGAA'
        frag3 = 'AAACAATTACCATACAACTCCAGTTTGAAGACTGTCCATAAGTCTTGGGAAAAGCTTACAGCTTCAGTGGAATCATGCTCTGTTCGTGTGAGTGATTCAAGCAGCGGAGCTCATAGGTTTGTAAACAAGGAGGATGGGTCTTCTCCAGCCGTGAAAAACGATTTCATCAACCGGCTACTTGAAACTGGTGCTACTGAGAGCTCCTCATCCAATATCTGCTCGAATCAGATGGAAGAAAATGGAGTGAATACGTCAAGCCAGATGACGCAAACCTTGTATAATCTAGTAGCCGCGACA'
        rec_site = {'original': 'GGTCTC', 'modified': 'GGTCTT'}
        prev_seq_len = 120
        overhangs = ['CTCG', 'CTCA', 'CTTA']
        end, start, overhangs = _get_segments_from_rec_site(frag5, frag3,
                                                            rec_site,
                                                            prev_seq_len,
                                                           overhangs=overhangs)
        for overhang in overhangs:
            assert not is_dna_palindrome(overhang)


    def test_get_pcr_segments(self):

        fragments = ['TATCGATCGATCGATGCTAGCTGATCGATCGAATCTACTACTACTACTAC',
                     'TTGCATGCTAGCTTTTATTTCGGGGTACTGGGATCTACTACTACTACTACTACTA']
        rec_sites = [{'original':'CGTCTC', 'modified':'CGTATC'}]
        seq = fragments[0] + rec_sites[0]['modified'] + fragments[1]
        segments = _get_pcr_segments(seq, rec_sites, fragments)
        assert segments == [{'start': 0, 'end': 54}, {'start': 51, 'end': 110}]

        seq = 'aggctgactatCGTCTCgtcagctagctgacgatcgatgctagctagctgactatagaggaaacccgtaacgctacgtacggctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq))

        rec_sites = [{'original':'CGTCTC', 'modified':'CGTATC'},
                     {'original':'GAGACC', 'modified':'GAGACA'},
                     {'original':'GAGACC', 'modified':'GAGACT'}]
        fragments = ['aggctgactat',
                     'gtcagctagctgacgatcgatgctagctagctgactatagaggaaacccgtaacgctacgtacggctagcaggtgctag',
                     'gggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg',
                     'gggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta']
        segments = _get_pcr_segments(seq, rec_sites, fragments)
        assert segments == [{'start': 0, 'end': 102, 'forward_min': 23},
                            {'start': 99, 'end': 180},
                            {'start': 177, 'end': 241}]

        # test with repeated overhangs
        seq = 'ATGGCGGCTGCTGCCTCACCATCTCCATGTTTCTCCAAAACCCTACCTCCATCTTCCTCCAAATCTTCCACCATTCTACCTAGATCTACCTTCTCTTTCCACAATCACCCACAAAAAGCCTCACCCCTTCATCTCATCCACGCTCAACATAATCGTCGTGGTTTTGCCGTTGCCAATGTCGTCATATCCACTACCACCCATAACGACGTTTCTGAACCTGAAACATTCGTTTCCCGTTTCGCCCCTGACGAACCCAGAAAGGGTTGTGATGTTCTTGTGGAGGCACTTGAAAGGGAAGGTGTTACGGATGTATTTGCATACCCAGGAGGTGCTTCTATGGAGATTCATCAAGCTTTGACACGTTCGAATATTATTCGTAATGTGCTACCACGTCATGAGCAAGGTGGTGTGTTTGCTGCAGAGGGTTACGCACGGGCTACTGGGTTCCCTGGTGTTTGCATTGCTACCTCTGGTCCCGGAGCTACAAATCTTGTTAGTGGTCTTGCGGATGCTTTGTTAGATAGTATTCCGATTGTTGCTATTACAGGTCAAGTGCCAAGGAGGATGATTGGTACTGATGCGTTCCAGGAAACGCCTATTGTTGAGGTAACGAGATCTATTACGAAGCATAATTATCTTGTTATGGATGTAGAAGATATTCCTAGGGTTGTTCGTGAAGCATTTTTTCTTGCGAAATCGGGACGGCCTGGCCCAGTTTTGATTGATGTACCTAAGGATATTCAGCAACAATTGGTGATACCTAATTGGGATCAGCCAATGAGGTTGCCTGGTTACATGTCTAGGTTACCTAAATTGCCTAATGAAATGCTTTTGGAACAAATTGTTAGGCTGATTTCCGAGTCGAAGAAGCCTGTTTTGTATGTGGGTGGTGGGTGTTCGCAATCAAGTGAGGAGCTGAGGCGATTTGTGGAGCTTACAGGTATTCCTGTAGCGAGTACTTTGATGGGTCTTGGAGCTTTTCCAACTGGGGATGAGCTTTCACTTCAAATGTTGGGTATGCATGGAACTGTGTATGCTAATTATGCTGTGGATAGTAGTGATTTGTTGCTTGCATTTGGGGTGAGGTTTGATGATCGAGTTACTGGTAAATTGGAAGCTTTTGCTAGTCGAGCGAAAATTGTCCACATTGATATTGATTCGGCAGAGATTGGAAAAAACAAGCAACCTCATGTTTCCATTTGTGCAGATATCAAGTTGGCATTACAGGGTTTGAATTCCATATTGGAGGGTAAAGAAGGTAAGATGAAGTTAGATTTTTCTGCCTGGAGGCAGGAGTTAACGGAGCAGAAGATGAAGTACCCACTGAATTTTAAGACTTTTGGTGATGCCATCCCTCCACAATATGCTATTCAGGTTCTTGATGAGTTAACTAACGGAAATGCCATTATTAGTACTGGTGTGGGGCAACACCAGATGTGGGCTGCCCAATACTATAAGTACAAAAAGCCACGCCAATGGTTGACATCTGGTGGATTAGGAGCAATGGGATTTGGTTTGCCTGCTGCTATAGGTGCGGCTGTTGGGAGGCCGGGTGAGATTGTGGTTGACATTGACGGTGATGGGAGTTTTATCATGAATGTGCAAGAGTTAGCAACAATTAAGGTGGAGAATCTCCCAGTTAAGATTATGTTGCTGAATAATCAACACTTGGGAATGGTGGTTCAATGGGAGGATCGATTCTATAAAGCTAACAGAGCACACACTTACTTGGGTGACCCTTCTAACGAGGAAGAGATCTTCCCTAATATGTTGAAATTTGCAGAGGCTTGTGGCGTACCTGCTGCAAGAGTGTCACACAGGGATGATCTTAGAGCTGCCATTCAAAAGATGTTAGACACTCCTGGGCCATACTTGTTGGATGTGATTGTACCTCATCAGGAGCACGTTCTACCTATGATTCCCAGCGGTGGTGCTTTCAAAGATGTGATCACGGAGGGCGACGGGAGATGTTCCTATTGA'
        rec_sites = [{'original': 'GAGACG', 'modified': 'GAGGCG'}, {'original': 'GAGACC', 'modified': 'GAGGCC'}, {'original': 'GCGATG', 'modified': 'GCGACG'}]
        fragments = ['ATGGCGGCTGCTGCCTCACCATCTCCATGTTTCTCCAAAACCCTACCTCCATCTTCCTCCAAATCTTCCACCATTCTACCTAGATCTACCTTCTCTTTCCACAATCACCCACAAAAAGCCTCACCCCTTCATCTCATCCACGCTCAACATAATCGTCGTGGTTTTGCCGTTGCCAATGTCGTCATATCCACTACCACCCATAACGACGTTTCTGAACCTGAAACATTCGTTTCCCGTTTCGCCCCTGACGAACCCAGAAAGGGTTGTGATGTTCTTGTGGAGGCACTTGAAAGGGAAGGTGTTACGGATGTATTTGCATACCCAGGAGGTGCTTCTATGGAGATTCATCAAGCTTTGACACGTTCGAATATTATTCGTAATGTGCTACCACGTCATGAGCAAGGTGGTGTGTTTGCTGCAGAGGGTTACGCACGGGCTACTGGGTTCCCTGGTGTTTGCATTGCTACCTCTGGTCCCGGAGCTACAAATCTTGTTAGTGGTCTTGCGGATGCTTTGTTAGATAGTATTCCGATTGTTGCTATTACAGGTCAAGTGCCAAGGAGGATGATTGGTACTGATGCGTTCCAGGAAACGCCTATTGTTGAGGTAACGAGATCTATTACGAAGCATAATTATCTTGTTATGGATGTAGAAGATATTCCTAGGGTTGTTCGTGAAGCATTTTTTCTTGCGAAATCGGGACGGCCTGGCCCAGTTTTGATTGATGTACCTAAGGATATTCAGCAACAATTGGTGATACCTAATTGGGATCAGCCAATGAGGTTGCCTGGTTACATGTCTAGGTTACCTAAATTGCCTAATGAAATGCTTTTGGAACAAATTGTTAGGCTGATTTCCGAGTCGAAGAAGCCTGTTTTGTATGTGGGTGGTGGGTGTTCGCAATCAAGTGAGGAGCT', 'ATTTGTGGAGCTTACAGGTATTCCTGTAGCGAGTACTTTGATGGGTCTTGGAGCTTTTCCAACTGGGGATGAGCTTTCACTTCAAATGTTGGGTATGCATGGAACTGTGTATGCTAATTATGCTGTGGATAGTAGTGATTTGTTGCTTGCATTTGGGGTGAGGTTTGATGATCGAGTTACTGGTAAATTGGAAGCTTTTGCTAGTCGAGCGAAAATTGTCCACATTGATATTGATTCGGCAGAGATTGGAAAAAACAAGCAACCTCATGTTTCCATTTGTGCAGATATCAAGTTGGCATTACAGGGTTTGAATTCCATATTGGAGGGTAAAGAAGGTAAGATGAAGTTAGATTTTTCTGCCTGGAGGCAGGAGTTAACGGAGCAGAAGATGAAGTACCCACTGAATTTTAAGACTTTTGGTGATGCCATCCCTCCACAATATGCTATTCAGGTTCTTGATGAGTTAACTAACGGAAATGCCATTATTAGTACTGGTGTGGGGCAACACCAGATGTGGGCTGCCCAATACTATAAGTACAAAAAGCCACGCCAATGGTTGACATCTGGTGGATTAGGAGCAATGGGATTTGGTTTGCCTGCTGCTATAGGTGCGGCTGTTGG', 'GGGTGAGATTGTGGTTGACATTGACGGTGATGGGAGTTTTATCATGAATGTGCAAGAGTTAGCAACAATTAAGGTGGAGAATCTCCCAGTTAAGATTATGTTGCTGAATAATCAACACTTGGGAATGGTGGTTCAATGGGAGGATCGATTCTATAAAGCTAACAGAGCACACACTTACTTGGGTGACCCTTCTAACGAGGAAGAGATCTTCCCTAATATGTTGAAATTTGCAGAGGCTTGTGGCGTACCTGCTGCAAGAGTGTCACACAGGGATGATCTTAGAGCTGCCATTCAAAAGATGTTAGACACTCCTGGGCCATACTTGTTGGATGTGATTGTACCTCATCAGGAGCACGTTCTACCTATGATTCCCAGCGGTGGTGCTTTCAAAGATGTGATCACGGAGG', 'GGAGATGTTCCTATTGA']
        segments = _get_pcr_segments(seq, rec_sites, fragments)
        assert segments == [{'start': 0, 'end': 921},
                            {'start': 918, 'end': 1550},
                            {'start': 1547, 'reverse_min': 1951, 'end': 1979}]

    def test_join_short_segments(self):
        segments = [(0, 10), (8, 33), (31, 60)]
        min_length = 11
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 33, 'forward_min': 18},
                                {'start': 31, 'end': 60}]

        segments = [(0, 34), (31, 40), (38, 80)]
        min_length = 11
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 34},
                                {'start': 31, 'end': 80, 'forward_min': 48}]
        segments = [(0, 30), (31, 60), (61, 70)]
        min_length = 11
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 30},
                                {'start': 31, 'reverse_min': 53, 'end': 70}]

        segments = [(0, 10), (11, 20), (21, 70)]
        min_length = 11
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 70, 'forward_min': 28}]

        segments = [(0, 50), (51, 60), (61, 70)]
        min_length = 11
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 50},
                                {'start': 51, 'reverse_min': 53, 'end': 70,
                                 'forward_min': 68}]

        segments = [(0, 4083), (4080, 4092), (4089, 7116), (7113, 7125), (7122, 7126)]
        min_length = 50
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 4083},
                                {'start': 4080, 'reverse_min': 7105,
                                 'end': 7126}]

    def test_rec_site_domestication(self):
        seq = SeqIO.read(os.path.join(TEST_DATA, 'CH3H.gb'), 'gb')
        category = '13-14-15-16 (CDS)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')

    def test_domestication_sp_category(self):
        seq = 'ATGATGATGGGCACTTCCTCTGTTTTTCTACTATTCCTTCTTTCTTTTCTTCTCCTTCTCCCG'
        seq += 'TCCCTCCTA'

        seqrec = SeqRecord(Seq(seq))
        category = 'CDS1 (B3)'
        oligo_pcrs = domesticate(seqrec, category, 'AATG', 'AGCC')[0]
        assert oligo_pcrs[0]['oligo_reverse'] == 'GCGCCGTCTCGCTCAGGCTGCTAGGAGGGACGGGAGAAGGA'

    def test_domestication_sintesis(self):
        seq = 'ATGATGATGGGCACTTCCTCTGTTTTTCTACTATTCCTTCTTTCTTTTCTTCTCCTTCTCCCG'
        seq += 'TCCCTCCTA'

        seqrec = SeqRecord(Seq(seq))
        category = 'SP (B3)'
        prefix = 'AATG'
        suffix = 'AGCC'
        enzymes = ()
        seq_sin, prep_seq = domesticate_for_synthesis(seqrec, category, prefix,
                                                      suffix, enzymes)
        result = 'GCGCCGTCTCGCTCGAATGATGATGATGGGCACTTCCTCTGTTTTTCTACTATTCCTTCT'
        result += 'TTCTTTTCTTCTCCTTCTCCCGTCCCTCCTAAGCCTGAGCGAGACGGCGC'
        assert seq_sin == result
        assert 'AATGATGATGATGGGCACTTCCTCTGTTTTTCTACTATTCCTTC' in str(prep_seq.seq)

        # with rec site
        seq = 'ATGATGATGGGCACTTCCTCTGTTTGGTCTCTATTCCTTCTTTCTTTTCTTCTCCTTCTCCCG'
        seq += 'TCCCTCCTA'

        seqrec = SeqRecord(Seq(seq))
        category = 'SP (B3)'
        seq_sin = domesticate_for_synthesis(seqrec, category, 'AATG', 'AGCC', enzymes)[0]
        result = 'GCGCCGTCTCGCTCGAATGATGATGATGGGCACTTCCTCTGTTTGGTCCCTATTCC'
        result += 'TTCTTTCTTTTCTTCTCCTTCTCCCGTCCCTCCTAAGCCTGAGCGAGACGGCGC'
        assert seq_sin == result

    def test_change_nucl_in_intron_rec_site(self):
        rec_sites = get_ret_sites(MANDATORY_DOMEST_ENZYMES)
        # regex with the sites to domesticate
        rec_sites_regex = '(' + '|'.join(rec_sites) + ')'
        rec_sites_regex = re.compile(rec_sites_regex, flags=re.IGNORECASE)
        rec_site = 'cGTCTC'
        res = change_nucl_in_intron_rec_site(rec_site, rec_sites_regex)
        assert res == 'aGTCTC'
        rec_site = 'CGTCTc'
        res = change_nucl_in_intron_rec_site(rec_site, rec_sites_regex)
        assert res == 'CGTCTa'

    def test_check_crispr_domestication_error(self):
        seq = "aatataatatgggacgagac".upper()
        category = 'D Target (B3c-B4-B5c)'
        prefix = "ATTG"
        suffix = "GTTT"
        self.assertRaises(ValueError, domestication_crispr, seq, category, prefix, suffix)

    def test_crispr_multiplexing_domestication(self):
        seq = "CGCTAGGATCGGACTGACGA"
        category = "CRISPR Multiplexing"
        prefix = ""
        suffix = ""
        seq = domestication_crispr(seq, category, prefix, suffix)
