import os.path
from cStringIO import StringIO

from django.test import TestCase, Client
from django.core.urlresolvers import reverse
from django.core.files.uploadedfile import SimpleUploadedFile
from django.conf import settings as proj_settings

from Bio import SeqIO
from Bio.Seq import Seq

import goldenbraid
from goldenbraid.views.feature_views import (FeatureForm,
                                             get_prefix_and_suffix,
                                            _choose_rec_sites,
                                            _pref_suf_index_from_rec_sites,
                                            _get_pref_suff_from_index,
                                            add_feature)
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD
from goldenbraid.models import Feature
from goldenbraid.settings import DB
from goldenbraid.tags import (VECTOR_TYPE_NAME, ENZYME_IN_TYPE_NAME,
                              ENZYME_OUT_TYPE_NAME)

TEST_DATA = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                         'goldenbraid', 'tests', 'data')


class FeatureTestViews(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_add_feature_form(self):
        test_data = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                                 'goldenbraid', 'tests', 'data')
        # test of the form
        gb_path = os.path.join(test_data, 'pAn11_uniq.gb')
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_in': 'AagI',
                     'enzyme_out': 'AaaI',
                     'resistance': 'pepe'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        form.is_valid()
        self.assertTrue(form.is_valid())

        # test of the form with blanck values
        gb_path = os.path.join(test_data, 'pAn11_uniq.gb')
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
        gb_path = os.path.join(test_data, 'pAn11_uniq.gb')
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
        gb_path = os.path.join(TEST_DATA, 'pAn11_uniq.gb')
        client = Client()
        url = reverse('add_feature')
        response = client.post(url, {'name': 'vector1',
                                     'type': VECTOR_TYPE_NAME,
                                     'description': 'vector1 desc',
                                     'reference': 'vector1 ref',
                                     'enzyme_in': 'BsaI',
                            'enzyme_out': 'AamI,AauI',
                            'resistance': 'vector1_resistance',
                            'gbfile': open(gb_path)})
        assert response.status_code == 200
        # TODO url to genbank file
        # response = client.get('/media/genbank_files/pAn11.gb')

        feat = Feature.objects.using(DB).get(uniquename='pAn11_uniq')
        assert feat.name == 'vector1'
        assert  feat.props == {u'Enzyme_in': [u'BsaI'],
                               u'Enzyme_out': [u'AamI', u'AauI'],
                               u'Description': [u'vector1 desc'],
                               u'Reference': [u'vector1 ref'],
                               u'Resistance': [u'vector1_resistance']}

        os.remove(os.path.join(proj_settings.MEDIA_ROOT,
                               feat.genbank_file.name))


    def test_get_prefix_and_suffix(self):
        'it tests get a suffix and prefix test'
        gb_path = os.path.join(TEST_DATA, 'pAn11_uniq.gb')
        seq = SeqIO.read(gb_path, 'gb')
        seq = seq.seq
        assert ('AATG', 'GCTT') == get_prefix_and_suffix(seq, 'BsaI')

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
        add_feature(DB, name, type_name, vector, open(genbank), props)
        feat = Feature.objects.using(DB).get(uniquename='pANT1_uniq')
        assert feat.uniquename == "pANT1_uniq"
        os.remove(os.path.join(proj_settings.MEDIA_ROOT,
                               feat.genbank_file.name))


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
        response = client.post(url)
        # print response
        assert """<p><label for="id_TER">Ter:</label>""" in str(response)
        assert """<select name="TER" id="id_TER">""" in str(response)

#        # 'It tests the basic typo of the form'
#        for uniq in ('pPE8', 'pANT1', 'pTnos', 'pDGB1_alpha1', 'pDGB1_alpha1R'):
#            feat = Feature.objects.using(DB).get(uniquename=uniq)
#            feat.genbank_file = File(open(os.path.join(TEST_DATA,
#                                                       '{0}.gb'.format(uniq))))
#            feat.save()

        client = Client()
        url = reverse('multipartite_view', kwargs={'multi_type': 'basic'})
        response = client.post(url, {"PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTnos',
                                     'Vector':'pDGB1_alpha1'})

        # print response
        assert 'error' not in response
        assert response.status_code == 200

        client = Client()
        url = reverse('multipartite_view_genbank', kwargs={'multi_type': 'basic'})
        response = client.post(url, {"PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTnos',
                                     'Vector': 'pDGB1_alpha1'})
        assert "LOCUS" in  str(response)
        client = Client()
        url = reverse('multipartite_view', kwargs={'multi_type': 'basic'})
        response = client.post(url, {"PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTno'})

        err1 = """<ul class="errorlist"><li>This field is required.</li></ul"""
        assert err1 in str(response)
        err2 = """<ul class="errorlist"><li>This feature does not exist in"""
        assert err2 in str(response)

        # reverse vector
        url = reverse('multipartite_view_genbank', kwargs={'multi_type': 'basic'})
        response = client.post(url, {"PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTnos',
                                     'Vector': 'pDGB1_alpha1R'})
        seqrec = SeqIO.read(StringIO(str(response)), 'gb')
        # seqrec = seqIO.read(tu fichero, 'genbank')

    def test_protocol_view(self):
        'it test that the protocol file is generated'
        client = Client()
        url = reverse('multipartite_protocol_view')
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'assembled_seq':'aaa',
                                    'multi_type':'basic',
                                    "PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTnos',
                                     'Vector':'pDGB1_alpha1'})
        assert "75 ng of pPE8" in str(response)
