'''
Created on 2013 ots 7

@author: peio
'''
from django.test import TestCase
import os.path

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import goldenbraid
from goldenbraid.domestication import (domesticate, _join_segments,
                                       _get_pcr_segments,
                                       _get_segments_from_rec_site)
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD
TEST_DATA = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                         'goldenbraid', 'tests', 'data')


class DomesticationTest(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_domestication(self):
        seq = 'aggctgactatgtcagctagctgacgatcgatgctagctagctgactagctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq))
        category = '13-14-15-16 (CDS)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[0] == {'pcr_product': 'GCGCCGTCTCGCTCGAAGGCTGACTATGTCAGCTAGCTGACGATCGATGCTAGCTAGCTGACTAGCTAGCAGGTGCTAGGAGACAGCGAGACGGCGC',
                                'oligo_reverse': 'GCGCCGTCTCGCTGTCTCCTAGCACCTGCTA',
                                'oligo_forward': 'GCGCCGTCTCGCTCGAAGGCTGACTATGTCAGCTAG'}
        assert oligo_pcrs[1] == {'pcr_product': 'GCGCCGTCTCGACAGGGTCATGCTAGCTTCAGCTAGCTGATCGATCGACTAGCTGATCGATCTGATCGATGCTAGCTAGCTGTACGGAGACTGCGAGACGGCGC',
                                 'oligo_reverse': 'GCGCCGTCTCGCAGTCTCCGTACAGCTAGCT',
                                 'oligo_forward': 'GCGCCGTCTCGACAGGGTCATGCTAGCTTCA'}
        domesticated_seq = domesticate(seq, category, 'CCAT', 'AATG')[1]
        assert ('TCGCATGCTCCCGGCCGCCATGGCGGCCGCGGGAATTCGATGGGCGATGAGTGGTCTCG') in domesticated_seq.seq
        assert ('TCGCATGCTCCCGGCCGCCATGGCGGCCGCGGGAATTCGATGGGCGATGAGTGGTCTCGT') not in domesticated_seq.seq

    def test_domestication_short_segments(self):
        seq = 'aggctgactatgtcagctaGAGACCgctgacgatcgatgctagctagctgactagctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq))
        category = '13-14-15-16 (CDS)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[0] == {'pcr_product': 'GCGCCGTCTCGCTCGAAGGCTGACTATGTCAGCTAGAGATCGCTGACGATCGATGCTAGCTAGCTGACTAGCTAGCAGGTGCTAGGAGACAGCGAGACGGCGC',
                                 'oligo_reverse': 'GCGCCGTCTCGCTGTCTCCTAGCACCTGCTA',
                                 'oligo_forward': 'GCGCCGTCTCGCTCGAAGGCTGACTATGTCAGCTAGAGATC'}

        seq = 'aggctgactatgtcagctagctgacgatcgatgctagctagctgactatagaggaaacccgtaacgctacgtacgCGTCTCgctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq))
        category = '13-14-15-16 (CDS)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[1] == {'pcr_product': 'GCGCCGTCTCGCTTGCTAGCAGGTGCTAGGAGACAGGGTCATGCTAGCTTCAGCTAGCTGATCGATCGACTAGCTGATCGATCTGATCGATGCTAGCTAGCTGTACGGAGACTGCGAGACGGCGC',
                                 'oligo_reverse': 'GCGCCGTCTCGCAGTCTCCGTACAGCTAGCT',
                                 'oligo_forward': 'GCGCCGTCTCGCTTGCTAGCAGGTGCTAGGAGACAG'}

        seq = 'aggctgactatgtcagctagctgacgatcgatgctagctagctgactatagaggaaacccgtaacgctacgtacggctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgctgcgtcgggcgtctgatgctagctagctCGTCTCgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq))
        category = '13-14-15-16 (CDS)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[2] == {'pcr_product': 'GCGCCGTCTCGACTGGGTCATGCTAGCTGATCTGATCGCTGCGTCGGGCGTCTGATGCTAGCTAGCTCGTATCGTACGTCATCTTTTCAGTCGATCTAAATGCGAGCGAGACGGCGC',
                                 'oligo_reverse': 'GCGCCGTCTCGCTCGCATTTAGATCGACTGAAAAGATGACGTACGATA',
                                 'oligo_forward': 'GCGCCGTCTCGACTGGGTCATGCTAGCTGAT'}

        gb_path = os.path.join(TEST_DATA, 'GB_DOMEST_15.gb')
        seq = SeqIO.read(gb_path, 'gb')
        category = '02 (OP)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[1]['oligo_reverse'] == 'GCGCCGTCTCGCTCGCATTCGCGACCACTCGTCGCCCAT'
        assert oligo_pcrs[1]['oligo_forward'] == 'GCGCCGTCTCGACAACTCATCGACCATCACTA'

    def test_get_segments_from_rec_site(self):
        frag5 = 'TATCGATCGATCGATGCTAGCTGATCGATCGAATCTACTACTACTACTAC'
        rec_site = {'original': 'CGTCTC', 'modified': 'CGTATC'}
        prev_seq_len = 0
        end, start = _get_segments_from_rec_site(frag5, rec_site, prev_seq_len)
        assert (51, 54) == (end, start)

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
        assert segments == [{'start': 0, 'end': 102, 'forward_min': 15},
                            {'start': 99, 'end': 180},
                            {'start': 177, 'end': 241}]

    def test_join_short_segments(self):
        segments = [(0, 10), (8, 33), (31, 60)]
        min_length = 11
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 33, 'forward_min': 10},
                                {'start': 31, 'end': 60}]

        segments = [(0, 34), (31, 40), (38, 80)]
        min_length = 11
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 34},
                                {'start': 31, 'end': 80, 'forward_min': 40}]
        segments = [(0, 30), (31, 60), (61, 70)]
        min_length = 11
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 30},
                                {'start': 31, 'reverse_min': 61, 'end': 70}]

        segments = [(0, 10), (11, 20), (21, 70)]
        min_length = 11
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 70, 'forward_min': 20}]

        segments = [(0, 50), (51, 60), (61, 70)]
        min_length = 11
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 50},
                                {'start': 51, 'reverse_min': 61, 'end': 70,
                                 'forward_min': 60}]

        segments = [(0, 4083), (4080, 4092), (4089, 7116), (7113, 7125), (7122, 7126)]
        min_length = 50
        new_segments = _join_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 4083},
                                {'start': 4080, 'reverse_min': 7113,
                                 'end': 7126}]
