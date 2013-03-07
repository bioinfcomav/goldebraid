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
                              ENZYME_OUT_TYPE_NAME, MODULE_TYPE_NAME)

TEST_DATA = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                         'goldenbraid', 'tests', 'data')


class FeatureTestViews(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_feature_page(self):
        client = Client()
        url = reverse('feature_view', kwargs={'uniquename': 'pAn11'})
        response = client.get(url)
        assert response.status_code == 200
        assert "Feature pAn11" in str(response)

    def test_add_feature_form(self):
        test_data = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                                 'goldenbraid', 'tests', 'data')
        # test of the form
        gb_path = os.path.join(test_data, 'pAn11_uniq.gb')
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': 'CDS', 'vector':'pDGB1_alpha1'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertTrue(form.is_valid())

        # test of the form with blanck values
        gb_path = os.path.join(test_data, 'pAn11_uniq.gb')
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': 'CDS', 'vector':'pDGB1_alpha1'}
        uploaded_fhand = open(gb_path)
        file_dict = {}
        form = FeatureForm(post_dict, file_dict)
        self.assertFalse(form.is_valid())

        # test of the form with wrong type
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': 'vecto'}
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

    def test_add_feature_view(self):
        # test of the form page
        # test of the form
        gb_path = os.path.join(TEST_DATA, 'pAn11_uniq.gb')
        client = Client()
        url = reverse('add_feature')
        response = client.post(url, {'name': 'vector1',
                                     'type': MODULE_TYPE_NAME,
                                     'description': 'vector1 desc',
                                     'reference': 'vector1 ref',
                                     'vector': 'pDGB1_omega1R',
                                     'gbfile': open(gb_path)})
        assert response.status_code == 200
        # TODO url to genbank file
        # response = client.get('/media/genbank_files/pAn11.gb')

        feat = Feature.objects.using(DB).get(uniquename='pAn11_uniq')
        assert feat.name == 'vector1'
        assert  feat.props == {u'Description': [u'vector1 desc'],
                               u'Reference': [u'vector1 ref']}

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

    def test_search_feature(self):
        client = Client()
        url = reverse('search_features')
        response = client.get(url)
        assert response.status_code == 200
        assert  "<option value=" in str(response)

        response = client.post(url, {'name_or_description': 'pAn11'})
        assert response.status_code == 302

        response = client.post(url, {'kind': 'TER'})
        assert response.status_code == 200
        assert  "<td>Agrobacterium tumefaciens terminator" in str(response)


class MultipartiteFreeTestViews(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_view(self):
        client = Client()
        url = reverse('multipartite_view_free')
        response = client.get(url)
        assert "pDGB2_alpha1R" in str(response)

        url = reverse('multipartite_view_free', kwargs={'form_num': '1'})

        response = client.post(url, {'vector': 'pDGB2_alpha1R',
                                     'part_1': 'pP2A11'})
        assert "An11" in str(response)

        url = reverse('multipartite_view_free', kwargs={'form_num': '2'})
        response = client.post(url, {'vector': 'pDGB2_alpha1R',
                                     'part_1': 'pP2A11',
                                     'part_2': 'pLuciferas'})
        assert 'feature does not exist' in str(response)

        response = client.post(url, {'vector': 'pDGB2_alpha1R',
                                     'part_1': 'pP2A11',
                                     'part_2': 'pLuciferase'})
        assert "pT35S" in str(response)

        response = client.post(url, {'vector': 'pDGB2_alpha1R',
                                     'part_1': 'pP2A11',
                                     'part_2': 'pLuciferase',
                                     'part_3': 'pT35S'})

        assert  "<p>You have assembled in the GoldenBraid" in str(response)

          # reverse vector
        url = reverse('multipartite_view_free_genbank')
        response = client.post(url, {'part_1': 'pP2A11',
                                     'part_2': 'pMYB12',
                                     'part_3': 'pTerm2A11',
                                     'vector': 'pDGB1_alpha1R'})

        assert response.status_code == 200

        seqrec1 = SeqIO.read(StringIO(str(response)), 'gb')
        multipartite_free_seq1 = str(seqrec1.seq)
        gb_path = os.path.join(TEST_DATA, 'pEGBMybrev_uniq.gb')
        seqrec2 = SeqIO.read(gb_path, 'gb')
        multipartite_free_seq2 = str(seqrec2.seq)[4:]
        multipartite_free_seq2 += str(seqrec2.seq)[:4]

        assert multipartite_free_seq1 == multipartite_free_seq2



    def test_genbank_view(self):
        'it test that the genbank file is generated'
        client = Client()
        url = reverse('multipartite_view_free_genbank')
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'assembled_seq':'aaa',
                                     'vector':'pDGB1_omega1',
                                     'part_1': 'pPE8',
                                     'part_2': 'pANT1',
                                     'part_3': 'pTnos'})

        assert  'LOCUS' in str(response)


    def test_protocol_view(self):
        'it test that the protocol file is generated'
        client = Client()
        url = reverse('multipartite_view_free_protocol')
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'assembled_seq':'aaa',
                                     'vector':'pDGB1_omega1',
                                     'part_1': 'pPE8',
                                     'part_2': 'pANT1',
                                     'part_3': 'pTnos'})

        assert "75 ng of pPE8" in str(response)




class MultipartiteTestViews(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_empty_type(self):
        client = Client()
        url = reverse('multipartite_view', kwargs={'multi_type': ''})
        response = client.get(url)
        assert "/do/multipartite/basic" in response.content

    def test_basic_type(self):
        'It tests the basic typo of the form'
        client = Client()
        url = reverse('multipartite_view', kwargs={'multi_type': 'basic'})
        response = client.post(url)
        assert """<p><label for="id_TER">Ter:</label>""" in str(response)
        assert """<select name="TER" id="id_TER">""" in str(response)
        assert """<option value="pDGB1_alpha1R">pDGB1_alpha""" in str(response)
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


        # forward vector
        url = reverse('multipartite_view_genbank', kwargs={'multi_type': 'basic'})
        response = client.post(url, {"PROM+UTR+ATG": 'pP35S',
                                     "CDS": 'pMYB12',
                                     "TER": 'pTnos',
                                     'Vector': 'pDGB1_omega2'})

        seqrec1 = SeqIO.read(StringIO(str(response)), 'gb')
        multipartite_seq1 = str(seqrec1.seq)
        gb_path = os.path.join(TEST_DATA, 'pEGBMyb_uniq.gb')
        seqrec2 = SeqIO.read(gb_path, 'gb')
        multipartite_seq2 = str(seqrec2.seq)
        assert multipartite_seq1 == multipartite_seq2


        # reverse vector
        url = reverse('multipartite_view_genbank', kwargs={'multi_type': 'basic'})
        response = client.post(url, {"PROM+UTR+ATG": 'pP2A11',
                                     "CDS": 'pMYB12',
                                     "TER": 'pTerm2A11',
                                     'Vector': 'pDGB1_alpha1R'})

        assert response.status_code == 200

        seqrec1 = SeqIO.read(StringIO(str(response)), 'gb')
        multipartite_seq1 = str(seqrec1.seq)
        gb_path = os.path.join(TEST_DATA, 'pEGBMybrev_uniq.gb')
        seqrec2 = SeqIO.read(gb_path, 'gb')
        multipartite_seq2 = str(seqrec2.seq)[4:]
        multipartite_seq2 += str(seqrec2.seq)[:4]

        assert multipartite_seq1 == multipartite_seq2

    def test_protocol_view(self):
        'it test that the protocol file is generated'
        client = Client()
        url = reverse('multipartite_view_protocol')
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'assembled_seq':'aaa',
                                    'multi_type':'basic',
                                    "PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTnos',
                                     'Vector':'pDGB1_alpha1'})
        assert "75 ng of pPE8" in str(response)

    def test_genbank_view(self):
        'it test that the protocol file is generated'
        client = Client()
        url = reverse('multipartite_view_genbank', kwargs={'multi_type':
                                                           'basic'})
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'assembled_seq':'aaa',
                                    'multi_type':'basic',
                                    "PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTnos',
                                     'Vector':'pDGB1_alpha1'})
        assert  'LOCUS' in str(response)


class BipartiteViewTest(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_bipartite(self):
        client = Client()
        # do initial
        url = reverse('bipartite_view')
        response = client.get(url)
        assert """<option value="GB0125">GB0125</option>""" in str(response)

        # do page 1
        url = reverse('bipartite_view', kwargs={'form_num': '1'})
        response = client.post(url, {'part_1': 'GB0125'})
        longline1 = """<input name="part_1" value="GB0125" readonly="True" maxlength="100" """
        longline1 += """type="text" id="id_part_1" />"""
        assert longline1 in str(response)
        assert """<p><label for="id_part_2">Part 2:</label>""" in str(response)

        # do page 2
        url = reverse('bipartite_view', kwargs={'form_num': '2'})
        response = client.post(url, {'part_1': 'GB0125', 'part_2': 'GB0126'})
        longline2 = """<input name="part_2" value="GB0126" readonly="True" """
        longline2 += """maxlength="100" type="text" id="id_part_2" />"""
        assert longline2 in str(response)
        assert "pDGB1_omega1" in str(response)

        # do page 3
        url = reverse('bipartite_view', kwargs={'form_num': '3'})
        response = client.post(url, {'part_1': 'GB0125', 'part_2': 'GB0126',
                                     'Vector': 'pDGB1_omega1'})
        assert """<INPUT type="hidden" name="Vector" value="pDGB1_omega1">""" in str(response)
        assert """ <p>The resulted sequence of the assembly is""" in str(response)

        # forward vector
        url = reverse('bipartite_view_genbank')
        response = client.post(url, {'part_1': 'GB0129',
                                     'part_2': 'GB0131',
                                     'Vector':'pDGB1_alpha1'})

        assert response.status_code == 200

        seqrec1 = SeqIO.read(StringIO(str(response)), 'gb')
        bipartite_seq1 = str(seqrec1.seq)
        gb_path = os.path.join(TEST_DATA, 'pEGBRosDelMyb.gb')
        seqrec2 = SeqIO.read(gb_path, 'gb')
        bipartite_seq2 = str(seqrec2.seq)
        assert bipartite_seq1 == bipartite_seq2

        # check bipartite_view_genbank
    def test_genbank_view(self):
        'it test that the genbank file is generated'
        client = Client()
        url = reverse('bipartite_view_genbank')
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'assembled_seq':'aaa',
                                     'part_1': 'GB0125',
                                     'part_2': 'GB0126',
                                     'Vector':'pDGB1_omega1'})
        assert  'LOCUS' in str(response)


    # check bipartite_view_protocol
    def test_protocol_view(self):
        'it test that the protocol file is generated'
        client = Client()
        url = reverse('bipartite_view_protocol')
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'assembled_seq':'aaa',
                                     'part_1': 'GB0125',
                                     'part_2': 'GB0126',
                                     'Vector':'pDGB1_omega1'})
        assert "75 ng of GB0125" in str(response)


class DomesticationViewTest(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_domestication(self):
        client = Client()
        # do initial
        url = reverse('domestication_view')
        response = client.get(url)
        assert ("""<option value="12 (NT)">12 (NT)</option>""") in str(response)

        # send data to formulary to test validations
        gb_path = os.path.join(TEST_DATA, 'domseq.gb')

        # add seq and category
        response = client.post(url, {'seq': open(gb_path),
                                     'category': '12 (NT)'})
        assert """<ul class="errorlist"><li>This field is required.</li></ul>""" not in str(response)

        # not add a sequence
        response = client.post(url, {'seq': '',
                                     'category': '12 (NT)'})
        assert """<ul class="errorlist"><li>This field is required.</li></ul>""" in str(response)

        # add category, prefix and suffix
        response = client.post(url, {'seq': open(gb_path),
                                     'prefix': 'ggac', 'suffix': 'cgtc', 'category': '12 (NT)'})
        assert """<ul class="errorlist"><li>Can not use category and prefix/suffix simoultaneously</li></ul>"""in str(response)

        # add category and suffix
        response = client.post(url, {'seq': open(gb_path),
                                     'prefix': '', 'suffix': 'cgtc', 'category': '12 (NT)'})
        assert """<ul class="errorlist"><li>Can not use category and prefix/suffix simoultaneously</li></ul>"""in str(response)

        # add suffix
        response = client.post(url, {'seq': open(gb_path),
                                     'prefix': '', 'suffix': 'cgtc', 'category': ''})
        assert """<ul class="errorlist"><li>You must provide prefix and suffix together</li></ul>""" in str(response)

        # not add category nor prefix and suffix
        response = client.post(url, {'seq': open(gb_path),
                                     'prefix': '', 'suffix': '', 'category': ''})
        assert """<ul class="errorlist"><li>At least we need category or prefix/suffix pair</li></ul>""" in str(response)

        # check that uses validators
        response = client.post(url, {'seq': open(gb_path),
                                     'category': '13-14-15-16 (CDS)'})
        assert 'The provided seq must start with start' in str(response)

        response = client.post(url, {'seq': open(gb_path),
                                     'category': '12-13 (GOI)'})
        assert 'The provided seq must have less' in str(response)

        # sequence start with atg
        fasta_path = os.path.join(TEST_DATA, 'domseqatg.fasta')
        response = client.post(url, {'seq': open(fasta_path),
                                     'category': '13 (SP)'})
        assert 'The provided seq must start with start' not in str(response)


    def test_genbank_view(self):
        'it test that the genbank file is generated'
        client = Client()
        url = reverse('domestication_view_genbank')
        response = client.get(url)
        assert response.status_code == 400
        response = client.post(url, {'seq': 'gagaggggggggagagagattcccctctccccccccccccccccccccccccccccccccccccctttgacctcgaaacgccccc',
                                     'prefix': 'ggag',
                                     'suffix': 'aatg',
                                     'category': '01-02-03-11-12 (PROM+UTR+ATG)'})
        assert  'LOCUS' in str(response)

    # check bipartite_view_protocol
    def test_protocol_view(self):
        'it test that the protocol file is generated'
        client = Client()
        url = reverse('domestication_view_protocol')
        response = client.get(url)
        assert response.status_code == 400
        response = client.post(url, {'seq': 'gagaggggggggagagagattcccctctccccccccccccccccctccccccccccccccccccccccccccctttgacctcgaaacgccccc',
                                     'prefix': 'ggag',
                                     'suffix': 'aatg',
                                     'category': '01-02-03-11-12 (PROM+UTR+ATG)'})
        assert "Oligo forward: GCGCCGTCTCGCTCGGGAGGAGAGGGGGGGGAGAGAGAT" in str(response)
