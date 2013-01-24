import os.path

from django.test import TestCase, Client
from django.core.urlresolvers import reverse
from django.core.files.uploadedfile import SimpleUploadedFile

from Bio import SeqIO
from Bio.Seq import Seq

import goldenbraid
from goldenbraid.views.feature_views import (FeatureForm,
                                             _get_prefix_an_suffix,
                                            _choose_rec_sites,
    _pref_suf_from_rec_sites)
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD
from goldenbraid.models import Feature
from goldenbraid.settings import DB
from goldenbraid.tags import VECTOR_TYPE_NAME

TEST_DATA = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                                 'goldenbraid', 'tests', 'data')


class FeatureTestViews(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_add_feature_form(self):
        test_data = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                                 'goldenbraid', 'tests', 'data')
        # test of the form
        gb_path = os.path.join(test_data, 'pAn11.gb')
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_in': 'AagI',
                     'enzyme_out': 'AaaI',
                     'resistance': 'pepe'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        form.is_valid()
        print form.errors
        self.assertTrue(form.is_valid())

        # test of the form with blanck values
        gb_path = os.path.join(test_data, 'pAn11.gb')
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_out': '',
                     'enzyme_in': 'vector1_enz_in', 'resistance': ''}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertFalse(form.is_valid())

        # test of the form with wrong type
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': 'vecto', 'enzyme_out': 'vector1_enz_out'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertFalse(form.is_valid())
        assert form.errors.get('type')

        # vector does not exist
        # test of the form with wrong type
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_out': 'vector1_enz_out',
                     'vector': 'vector1'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertFalse(form.is_valid())
        assert form.errors.get('vector')

        # enzyme_in not added. fails because is required
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_out': 'vector1_enz_out'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertFalse(form.is_valid())
        enz_in = 'enzyme_in'
        assert 'vector must have a enzyme in' in str(form.errors.get(enz_in))

        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_out': 'vector1_enz_out',
                     'enzyme_in': 'no_exist'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertFalse(form.is_valid())
        enz_in = 'enzyme_in'
        assert 'This enzyme: no_exist is not a' in str(form.errors.get(enz_in))

        # enzyme_out with two enzymes
        gb_path = os.path.join(test_data, 'pAn11.gb')
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_in': 'AagI',
                     'enzyme_out': 'AamI,AauI', 'resistance': 'aa'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertTrue(form.is_valid())

        # resistance not added to a vector
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_out': 'vector1_enz_out',
                     'enzyme_in': 'vector1_enz_in', 'resistance': ''}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertFalse(form.is_valid())
        err_str = form.errors.get('resistance')
        assert 'A vector must have a resistance' in err_str

    def test_add_feature_view(self):
        # test of the form page
        # test of the form
        gb_path = os.path.join(TEST_DATA, 'pAn11.gb')
        client = Client()
        url = reverse('add_feature')
        response = client.post(url, {'name': 'vector1',
                                     'type': VECTOR_TYPE_NAME,
                                     'description': 'vector1 desc',
                                     'reference': 'vector1 ref',
                                     'enzyme_in': 'AagI',
                            'enzyme_out': 'AamI,AauI',
                            'resistance': 'vector1_resistance',
                            'gbfile': open(gb_path)})
        feat = Feature.objects.using(DB).get(uniquename='pAn11')
        assert feat.name == 'vector1'
        assert  feat.props == {u'Enzyme_in': [u'AagI'],
                               u'Enzyme_out': [u'AamI', u'AauI'],
                               u'Description': [u'vector1 desc'],
                               u'Reference': [u'vector1 ref'],
                               u'Resistance': [u'vector1_resistance']}

    def test_get_prefix_and_suffix(self):
        'it tests get a suffix and prefix test'
        gb_path = os.path.join(TEST_DATA, 'pAn11.gb')
        seq = SeqIO.read(gb_path, 'gb')
        seq = seq.seq
        assert ('AATG', 'GCTT') == _get_prefix_an_suffix(seq, 'BsaI')

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
        result = ('atct', 'ctga')
        assert _pref_suf_from_rec_sites(seq, forw_site, rev_site, rec_site,
                                       forw_cut_delta, rev_cut_delta) == result

        forw_site = 3
        rev_site = 30
        rec_site = 'aaaaaa'
        forw_cut_delta = 1
        rev_cut_delta = 5
        result = ('actg', 'gact')
        assert _pref_suf_from_rec_sites(seq, forw_site, rev_site, rec_site,
                                       forw_cut_delta, rev_cut_delta) == result

        forw_site = 25
        rev_site = 3
        rec_site = 'aaaaaa'
        forw_cut_delta = 1
        rev_cut_delta = 5
        result = ('cgat', 'atat')
        assert _pref_suf_from_rec_sites(seq, forw_site, rev_site, rec_site,
                                       forw_cut_delta, rev_cut_delta) == result

        forw_site = 4082
        rev_site = 1039
        seq = 'CTCGAATGGAGAATTCAAGTCAAGAATCGCATCTCCGATCTGAAAATTCCGTTACATATGACTCCTCTTACCCGATCTACGCTATGGCTTTTTCATCCTTCACTTCTTCCCTCACAAACCGCCGCCGTCGACTTGCCGTCGGAAGCTTTATCGAAGAGTTCAACAATCGGGTTGATATTCTCTCTTTCGACGAAGATACCCTAACCCTTAAGCCCGTTCCAAATCTCTCTTTCGAACACCCTTATCCACCAACAAAGCTCATGTTTCATCCTAATCCTTCTGCTTCTCTCAAGACTAATGATATTCTTGCCTCTTCCGGCGACTACCTCCGGCTCTGGGATGTTACTGATACTTCCATTGAACCACTTTTCACTCTCAGTAACAATAAAACCAGTGAATACTGTGCTCCTTTGACGTCTTTTGATTGGAATGAAGTGGAGCCGAGAAGAATTGGTACTTCTAGTATAGACACTACTTGTACCATCTGGGATGTTGAAAAAGGAGTTGTGGAAACTCAATTGATAGCACATGACAAAGAGGTTTACGATATAGCTTGGGGTGAAGCTGGGGTTTTTGCGTCTGTTTCTGCTGATGGATCCGTTAGGATTTTTGATTTGAGAGATAAGGAACACTCGACGATTATTTATGAGAGCCCGAAACCGGATACGCCATTGTTGAGGTTGGCTTGGAACAAACAGGATTTGAGATACATGGCTACCATATTGATGGATAGCAACAAGATTGTGATCTTAGATATTAGATCTCCAGCAATGCCGGTGGCTGAACTGGAAAGGCATCAGGCGAGTGTGAATGCTATTGCTTGGGCTCCGCAGAGCTGTAGACATATTTGTTCTGGTGGGGATGACGGACAGGCGCTCATTTGGGAGTTGCCAACTGTTGCAGGGCCTAATGGGATTGATCCCATGTCAGTGTACACCGCCGGAGCTGAGATTAATCAAATTCAGTGGTCTGCTGCACAGCGTGATTGGATTGCCATTACGTTTTCTAACAAGTTGCAGCTGCTTAAAGTATGAGCTTCGAGACCACTCATCGCCCATCACTAGTGAATTCGCGGCCGCCTGCAGGTCGACCATATGGGAGAGCTCCCAACGCGTTGGATGCATAGCTTGAGTATTCTATAGTGTCACCTAAATAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGAAATTGTAAGCGTTAATATTTTGTTAAAATTCGCGTTAAATTTTTGTTAAATCAGCTCATTTTTTAACCAATAGGCCGAAATCGGCAAAATCCCTTATAAATCAAAAGAATAGACCGAGATAGGGTTGAGTGTTGTTCCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGGCCCACTACGTGAACCATCACCCTAATCAAGTTTTTTGGGGTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAGAAAGCGAAAGGAGCGGGCGCTAGGGCGCTGGCAAGTGTAGCGGTCACGCTGCGCGTAACCACCACACCCGCCGCGCTTAATGCGCCGCTACAGGGCGCGTCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTGTAATACGACTCACTATAGGGCGAATTGGGCCCGACGTCGCATGCTCCCGGCCGCCATGGCGGCCGCGGGAATTCGATGGGCGATGAGTGGT'
        rec_site = 'GGTCTC'
        forw_cut_delta = 1
        rev_cut_delta = 5
        result = ('AATG', 'GCTT')
        assert ('AATG', 'GCTT') == _pref_suf_from_rec_sites(seq, forw_site,
                                                            rev_site, rec_site,
                                                            forw_cut_delta,
                                                            rev_cut_delta)


class MultipartiteTestViews(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_empty_type(self):
        client = Client()
        url = reverse('multipartite_view', kwargs={'multi_type': ''})
        response = client.get(url)
        assert "<div id='main'>" in response.content

    def test_basic_type(self):
        'It tests the basic typo of the form'
        client = Client()
        url = reverse('multipartite_view', kwargs={'multi_type': 'basic'})
        response = client.get(url)
        print response

