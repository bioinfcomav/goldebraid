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

import os.path

from django.test import TestCase
from django.conf import settings as proj_settings

from Bio import SeqIO
from Bio.Seq import Seq

import goldenbraid
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD
from goldenbraid.views.feature import get_prefix_and_suffix, add_feature, \
    _get_pref_suff_from_index
from goldenbraid.utils import (_choose_rec_sites,
                               _pref_suf_index_from_rec_sites,
                               get_prefix_and_suffix_index)
from goldenbraid.tags import (VECTOR_TYPE_NAME, ENZYME_IN_TYPE_NAME,
                              ENZYME_OUT_TYPE_NAME, TU_TYPE_NAME, DERIVES_FROM)
from goldenbraid.models import Feature, Cv, Cvterm, FeatureRelationship
from goldenbraid.domestication import domesticate
from Bio.SeqRecord import SeqRecord
from django.contrib.auth.models import User

TEST_DATA = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                         'goldenbraid', 'tests', 'data')


class FeatureTest(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_get_prefix_and_suffix(self):
        'it tests get a suffix and prefix test'
        gb_path = os.path.join(TEST_DATA, 'pAn11_uniq.gb')
        seq = SeqIO.read(gb_path, 'gb')
        seq = seq.seq
        assert ('AATG', 'GCTT') == get_prefix_and_suffix(seq, 'BsaI')

        fasta_path = os.path.join(TEST_DATA, 'seq.fasta')
        seq = SeqIO.read(fasta_path, 'fasta')
        seq = seq.seq
        assert ('TGGA', 'AATG') == get_prefix_and_suffix(seq, 'BsaI')

        # no rec_sites
        fasta_path = os.path.join(TEST_DATA, 'seq2.fasta')
        seq = SeqIO.read(fasta_path, 'fasta')
        seq = seq.seq
        assert (None, None) == get_prefix_and_suffix(seq, 'BsaI')

        # no rec_sites
        fasta_path = os.path.join(TEST_DATA, 'seq_with_intron.fasta')
        seqrec = SeqIO.read(fasta_path, 'fasta')
        dom_seq = domesticate(seqrec, None, 'CCAT', 'AATG')[1]
        assert ('CCAT', 'AATG') == get_prefix_and_suffix(dom_seq.seq, 'BsaI')

    def xtest_daniels_promoter(self):
        fasta_path = os.path.join(TEST_DATA, 'prVvbHLH35v2.fa')
        seqrec = SeqIO.read(fasta_path, 'fasta')

        dom_seq = domesticate(seqrec, '01-02-03-11-12 (PROM+UTR+ATG)',
                              'GGAG', 'AATG')[1]
        print len(dom_seq), seqrec.seq.upper()
        print len(dom_seq), dom_seq.seq
        assert ('GGAG', 'AATG') == get_prefix_and_suffix(dom_seq.seq, 'BsaI')

    def test_choose_rec_sites(self):
        'it tests choose rec_sites func'
        forw_sites = [4083]
        rev_sites = [1039, 2472]
        result = (4083, 1039)

        assert _choose_rec_sites(forw_sites, rev_sites) == result

        forw_sites = [1039, 2472]
        rev_sites = [4083]
        result = (2472, 4083)

        assert _choose_rec_sites(forw_sites, rev_sites) == result

        forw_sites = [1039, 2472]
        rev_sites = [1500]
        result = (1039, 1500)

        assert _choose_rec_sites(forw_sites, rev_sites) == result

        forw_sites = [1500]
        rev_sites = [1039, 2472]
        result = (1500, 2472)

        assert _choose_rec_sites(forw_sites, rev_sites) == result

    def test_pref_suf_from_rec_sites(self):
        'it tests _pref_suf_from_rec_sites'
        #                    1         2         3
        #          0123456789012345678901234567890123456789
        seq = Seq('atctgcatcgactgactgactgatcgactgatcgatcgat')

        forw_site = 33
        rev_site = 20
        rec_site = 'aaaaaa'
        forw_cut_delta = 1
        rev_cut_delta = 5
        pref_size = rev_cut_delta - forw_cut_delta
        result = ('atct', 'ctga')
        p_idx, s_idx = _pref_suf_index_from_rec_sites(seq, forw_site, rev_site,
                                                      rec_site,
                                                      forw_cut_delta,
                                                      rev_cut_delta)

        assert _get_pref_suff_from_index(seq, p_idx, s_idx, pref_size) == \
                                                                        result

        forw_site = 3
        rev_site = 30
        rec_site = 'aaaaaa'
        forw_cut_delta = 1
        rev_cut_delta = 5
        pref_size = rev_cut_delta - forw_cut_delta
        result = ('actg', 'gact')
        p_idx, s_idx = _pref_suf_index_from_rec_sites(seq, forw_site, rev_site,
                                                      rec_site,
                                                      forw_cut_delta,
                                                      rev_cut_delta)
        assert _get_pref_suff_from_index(seq, p_idx, s_idx, pref_size) == \
                                                                        result

        forw_site = 25
        rev_site = 3
        rec_site = 'aaaaaa'
        forw_cut_delta = 1
        rev_cut_delta = 5
        pref_size = rev_cut_delta - forw_cut_delta
        result = ('cgat', 'atat')
        p_idx, s_idx = _pref_suf_index_from_rec_sites(seq, forw_site, rev_site,
                                                      rec_site,
                                                      forw_cut_delta,
                                                      rev_cut_delta)
        assert _get_pref_suff_from_index(seq, p_idx, s_idx, pref_size) == \
                                                                        result

        forw_site = 4082
        rev_site = 1039
        seq = 'CTCGAATGGAGAATTCAAGTCAAGAATCGCATCTCCGATCTGAAAATTCCGTTACATATGACTCCTCTTACCCGATCTACGCTATGGCTTTTTCATCCTTCACTTCTTCCCTCACAAACCGCCGCCGTCGACTTGCCGTCGGAAGCTTTATCGAAGAGTTCAACAATCGGGTTGATATTCTCTCTTTCGACGAAGATACCCTAACCCTTAAGCCCGTTCCAAATCTCTCTTTCGAACACCCTTATCCACCAACAAAGCTCATGTTTCATCCTAATCCTTCTGCTTCTCTCAAGACTAATGATATTCTTGCCTCTTCCGGCGACTACCTCCGGCTCTGGGATGTTACTGATACTTCCATTGAACCACTTTTCACTCTCAGTAACAATAAAACCAGTGAATACTGTGCTCCTTTGACGTCTTTTGATTGGAATGAAGTGGAGCCGAGAAGAATTGGTACTTCTAGTATAGACACTACTTGTACCATCTGGGATGTTGAAAAAGGAGTTGTGGAAACTCAATTGATAGCACATGACAAAGAGGTTTACGATATAGCTTGGGGTGAAGCTGGGGTTTTTGCGTCTGTTTCTGCTGATGGATCCGTTAGGATTTTTGATTTGAGAGATAAGGAACACTCGACGATTATTTATGAGAGCCCGAAACCGGATACGCCATTGTTGAGGTTGGCTTGGAACAAACAGGATTTGAGATACATGGCTACCATATTGATGGATAGCAACAAGATTGTGATCTTAGATATTAGATCTCCAGCAATGCCGGTGGCTGAACTGGAAAGGCATCAGGCGAGTGTGAATGCTATTGCTTGGGCTCCGCAGAGCTGTAGACATATTTGTTCTGGTGGGGATGACGGACAGGCGCTCATTTGGGAGTTGCCAACTGTTGCAGGGCCTAATGGGATTGATCCCATGTCAGTGTACACCGCCGGAGCTGAGATTAATCAAATTCAGTGGTCTGCTGCACAGCGTGATTGGATTGCCATTACGTTTTCTAACAAGTTGCAGCTGCTTAAAGTATGAGCTTCGAGACCACTCATCGCCCATCACTAGTGAATTCGCGGCCGCCTGCAGGTCGACCATATGGGAGAGCTCCCAACGCGTTGGATGCATAGCTTGAGTATTCTATAGTGTCACCTAAATAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGAAATTGTAAGCGTTAATATTTTGTTAAAATTCGCGTTAAATTTTTGTTAAATCAGCTCATTTTTTAACCAATAGGCCGAAATCGGCAAAATCCCTTATAAATCAAAAGAATAGACCGAGATAGGGTTGAGTGTTGTTCCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGGCCCACTACGTGAACCATCACCCTAATCAAGTTTTTTGGGGTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAGAAAGCGAAAGGAGCGGGCGCTAGGGCGCTGGCAAGTGTAGCGGTCACGCTGCGCGTAACCACCACACCCGCCGCGCTTAATGCGCCGCTACAGGGCGCGTCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTGTAATACGACTCACTATAGGGCGAATTGGGCCCGACGTCGCATGCTCCCGGCCGCCATGGCGGCCGCGGGAATTCGATGGGCGATGAGTGGT'
        rec_site = 'GGTCTC'
        forw_cut_delta = 1
        rev_cut_delta = 5
        pref_size = rev_cut_delta - forw_cut_delta
        result = ('AATG', 'GCTT')
        p_idx, s_idx = _pref_suf_index_from_rec_sites(seq, forw_site, rev_site,
                                                      rec_site,
                                                      forw_cut_delta,
                                                      rev_cut_delta)
        assert _get_pref_suff_from_index(seq, p_idx, s_idx, pref_size) == \
                                                                        result

    def test_add_feature(self):
        'It tests the add_feature function'
        name = 'test_name'
        type_name = VECTOR_TYPE_NAME
        vector = None
        genbank = os.path.join(TEST_DATA, 'pANT1_uniq.gb')
        props = {ENZYME_IN_TYPE_NAME: ['BsaI'],
                 ENZYME_OUT_TYPE_NAME: ['BsaI']}
        try:
            add_feature(name, type_name, vector, open(genbank), props,
                        'test_user')
            self.fail()
        except RuntimeError as error:
            if 'the given user does not exist' in str(error):
                pass
            else:
                raise
        User.objects.create_user(username='admin', email='admin@upv.es',
                                 password='top_secret')
        add_feature(name, type_name, vector, open(genbank), props, 'admin')

        feat = Feature.objects.get(uniquename='pANT1_uniq')
        assert FeatureRelationship.objects.filter(object=feat).count() == 0
        assert feat.uniquename == "pANT1_uniq"
        os.remove(os.path.join(proj_settings.MEDIA_ROOT,
                               feat.genbank_file.name))

        # add relationships to feature
        # first we need to add derives from to db as it is not in the fixtures
        name = 'GB_UA_17A'
        vector = 'pDGB1_alpha1'
        type_name = TU_TYPE_NAME
        genbank = os.path.join(TEST_DATA, 'GB_UA_17A.gb')
        props = {ENZYME_IN_TYPE_NAME: ['BsaI'],
                 ENZYME_OUT_TYPE_NAME: ['BsaI']}
        # add_feature(name, type_name, vector, open(genbank), props, 'admin')
        feat = Feature.objects.get(uniquename='GB_UA_17A')
        assert FeatureRelationship.objects.filter(object=feat).count() == 2
