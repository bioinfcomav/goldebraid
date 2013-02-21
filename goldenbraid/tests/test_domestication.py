'''
Created on 2013 ots 7

@author: peio
'''
from django.test import TestCase
from Bio.Seq import Seq

from goldenbraid.domestication import domesticate, _join_short_segments, \
    _get_pcr_segments, _get_segments_from_rec_site
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD
from Bio.SeqRecord import SeqRecord


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

    def test_domestication_short_segments(self):
        seq = 'aggctgactatgtcagctaGAGACCgctgacgatcgatgctagctagctgactagctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq))
        category = '13-14-15-16 (CDS)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        assert oligo_pcrs[0] == {'pcr_product': 'GCGCCGTCTCGCTCGAAGGCTGACTATGTCAGCTAGAGATCGCTGACGATCGATGCTAGCTAGCTGACTAGCTAGCAGGTGCTAGGAGACAGCGAGACGGCGC',
                                 'oligo_reverse': 'GCGCCGTCTCGCTGTCTCCTAGCACCTGCTA',
                                 'oligo_forward': 'GCGCCGTCTCGCTCGAAGGCTGACTATGTCAGCTAGAGAT'}

        seq = 'aggctgactatCGTCTCgtcagctagctgacgatcgatgctagctagctgactatagaggaaacccgtaacgctacgtacggctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq))
        category = '13-14-15-16 (CDS)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        # print oligo_pcrs
        seq = 'aggctgactatgtcagctagctgacgatcgatgctagctagctgactatagaggaaacccgtaacgctacgtacgCGTCTCgctagcaggtgctag'
        seq += 'GAGACCgggtcatgctagcttcagctagctgatcgatcgactagctgatcgatctgatcgatgctagctagctgtacg'
        seq += 'GAGACCgggtcatgctagctgatctgatcgatgctagctagctgtacgtcatcttttcagtcgatcta'
        seq = SeqRecord(Seq(seq))
        category = '13-14-15-16 (CDS)'
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG')[0]
        print oligo_pcrs






    def test_get_segments_from_rec_site(self):
        frag5 = 'TATCGATCGATCGATGCTAGCTGATCGATCGAATCTACTACTACTACTAC'
        rec_site = {'original':'CGTCTC', 'modified':'CGTATC'}
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
        new_segments = _join_short_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 33, 'forward_min': 10},
                                {'start': 31, 'end': 60}]

        segments = [(0, 34), (31, 40), (38, 80)]
        min_length = 11
        new_segments = _join_short_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 34},
                                {'start': 31, 'end': 80, 'forward_min': 40}]
        segments = [(0, 30), (31, 60), (61, 70)]
        min_length = 11
        new_segments = _join_short_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 30},
                                {'start': 31, 'reverse_min': 61, 'end': 70}]

        segments = [(0, 10), (11, 20), (21, 70)]
        min_length = 11
        new_segments = _join_short_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 70, 'forward_min': 20}]

        segments = [(0, 50), (51, 60), (61, 70)]
        min_length = 11
        new_segments = _join_short_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 50},
                                {'start': 51, 'reverse_min': 61, 'end': 70,
                                 'forward_min': 60}]
