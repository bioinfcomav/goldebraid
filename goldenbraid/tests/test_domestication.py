'''
Created on 2013 ots 7

@author: peio
'''
from django.test import TestCase
from Bio.Seq import Seq

from goldenbraid.domestication import domesticate, _join_short_segments
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
        oligo_pcrs = domesticate(seq, category, 'CCAT', 'AATG', 50)[0]
        assert oligo_pcrs[0] == {'pcr_product': 'GCGCCGTCTCTCTCGAAGGCTGACTATGTCAGCTAGCTGACGATCGATGCTAGCTAGCTGACTAGCTAGCAGGTGCTAGGAGACAGAGAGACGGCGC',
                             'oligo_reverse': 'GCGCCGTCTCTCTGTCTCCTAGCACCTGCT',
                        'oligo_forward': 'GCGCCGTCTCTCTCGAAGGCTGACTATGTCAGCTA'}
        assert oligo_pcrs[1] == {'pcr_product': 'GCGCCGTCTCTACAGGGTCATGCTAGCTTCAGCTAGCTGATCGATCGACTAGCTGATCGATCTGATCGATGCTAGCTAGCTGTACGGAGACTGAGAGACGGCGC',
                            'oligo_reverse': 'GCGCCGTCTCTCAGTCTCCGTACAGCTAGC',
                            'oligo_forward': 'GCGCCGTCTCTACAGGGTCATGCTAGCTTC'}

    def test_join_short_segments(self):
        segments = [(0, 10), (11, 30), (31, 60)]
        min_length = 11
        new_segments = _join_short_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 30, 'forward_min': 10},
                                {'start': 31, 'end': 60}]

        segments = [(0, 30), (31, 40), (41, 80)]
        min_length = 11
        new_segments = _join_short_segments(segments, min_length)
        assert new_segments == [{'start': 0, 'end': 30},
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
