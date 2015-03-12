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
from cStringIO import StringIO

from django.test import TestCase, Client
from django.core.urlresolvers import reverse
from django.core.files.uploadedfile import SimpleUploadedFile
from django.conf import settings as proj_settings

from Bio import SeqIO

import goldenbraid
from goldenbraid.views.feature_views import FeatureForm
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD
from goldenbraid.models import Feature
from goldenbraid.tags import VECTOR_TYPE_NAME, MODULE_TYPE_NAME

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
                     'type': 'CDS', 'vector': 'pDGB1_alpha1'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertTrue(form.is_valid())

        # test of the form with blanck values
        gb_path = os.path.join(test_data, 'pAn11_uniq.gb')
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': 'CDS', 'vector': 'pDGB1_alpha1'}
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

        # no login, no access
        response = client.post(url, {'name': 'vector1',
                                     'type': MODULE_TYPE_NAME,
                                     'description': 'vector1 desc',
                                     'reference': 'vector1 ref',
                                     'vector': 'pDGB1_omega1R',
                                     'gbfile': open(gb_path)})
        assert response.status_code == 302

        client.login(username='admin', password='password')
        # show form
        response = client.get(url)
        assert "pDGB1_alpha1" in str(response)

        # add a feature
        url = reverse('add_feature')
        response = client.post(url, {'name': 'vector1',
                                     'type': MODULE_TYPE_NAME,
                                     'description': 'vector1 desc',
                                     'reference': 'vector1 ref',
                                     'vector': 'pDGB1_omega1R',
                                     'gbfile': open(gb_path)})
        assert response.status_code == 302
        # TODO url to genbank file
        # response = client.get('/media/genbank_files/pAn11.gb')

        feat = Feature.objects.get(uniquename='pAn11_uniq')
        assert feat.name == 'vector1'
        assert feat.props == {u'Description': [u'vector1 desc'],
                              u'Reference': [u'vector1 ref']}

        # add a feature
        url = reverse('add_feature')
        gb_path = os.path.join(TEST_DATA, 'GB_DOMEST_15.gb')
        response = client.post(url, {'name': 'vector1',
                                     'type': 'TU',
                                     'description': 'vector1 desc',
                                     'reference': 'vector1 ref',
                                     'vector': 'pDGB1_alpha2',
                                     'gbfile': open(gb_path)})

        assert response.status_code == 200

        os.remove(os.path.join(proj_settings.MEDIA_ROOT,
                               feat.genbank_file.name))

    def test_search_feature(self):
        client = Client()
        url = reverse('search_features')
        response = client.get(url)
        assert response.status_code == 200
        assert "<option value=" in str(response)

        response = client.post(url, {'name_or_description': 'pAn11'})
        assert response.status_code == 302

        response = client.post(url, {'kind': 'TER'})
        assert response.status_code == 200
        assert "<td>This is a pGreen destiny vector of the" in str(response)

        client.login(username='test', password='testpass')
        response = client.post(url, {'only_user': True})
        assert response.status_code == 200
        assert 'pEGB 2A11:Myb12:Tnos' in str(response)


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

        assert "<p>You have assembled in the GoldenBraid" in str(response)

        # reverse vector
        url = reverse('multipartite_view_free_genbank')
        response = client.post(url, {'part_1': 'pP2A11',
                                     'part_2': 'pMYB12',
                                     'part_3': 'pTerm2A11',
                                     'vector': 'pDGB1_alpha1R'})

        assert response.status_code == 200

        seqrec1 = SeqIO.read(StringIO(str(response)), 'gb')
        assert seqrec1.name == 'GB_UA_8'
        multipartite_free_seq1 = str(seqrec1.seq)
        gb_path = os.path.join(TEST_DATA, 'pEGBMybrev_uniq.gb')
        seqrec2 = SeqIO.read(gb_path, 'gb')
        multipartite_free_seq2 = str(seqrec2.seq)[4:]
        multipartite_free_seq2 += str(seqrec2.seq)[:4]

        assert multipartite_free_seq1 == multipartite_free_seq2

        # with more than one part of the same type
        url = reverse('multipartite_view_free', kwargs={'form_num': '5'})
        response = client.post(url, {'part_1': 'pP2A11',
                                     'part_2': 'GB0365',
                                     'part_3': 'GB0653',
                                     'part_4': 'GB0655',
                                     'part_5': 'pT35S',
                                     'vector': 'pDGB1_alpha1'})
        assert "<p>Other.2:<a href='/feature/GB0655'>GB0655</a></p>" in  str(response)

    def test_genbank_view(self):
        'it test that the genbank file is generated'
        client = Client()
        url = reverse('multipartite_view_free_genbank')
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'assembled_seq': 'aaa',
                                     'vector': 'pDGB1_omega1',
                                     'part_1': 'pPE8',
                                     'part_2': 'pANT1',
                                     'part_3': 'pTnos'})
        assert 'GB_UA_8' in str(response)
        assert 'LOCUS' in str(response)

        response = client.post(url, {'assembled_seq': 'aaa',
                                     'vector': 'pDGB1_omega1',
                                     'part_1': 'pPE8',
                                     'part_2': 'pANT1',
                                     'part_3': 'pTnos'})
        assert 'GB_UA_9' in str(response)
        assert 'LOCUS' in str(response)

        # with more than one part of the same type
        response = client.post(url, {'part_1': 'pP2A11',
                                     'part_2': 'GB0365',
                                     'part_3': 'GB0653',
                                     'part_4': 'GB0655',
                                     'part_5': 'pT35S',
                                     'vector': 'pDGB1_alpha1'})
        assert '(pP2A11,GB0365,GB0653,GB0655,pT35S)pDGB1_alpha1' in str(response)

    def test_protocol_view(self):
        'it test that the protocol file is generated'
        client = Client()
        url = reverse('multipartite_view_free_protocol')
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'assembled_seq': 'aaa',
                                     'vector': 'pDGB1_omega1',
                                     'part_1': 'pPE8',
                                     'part_2': 'pANT1',
                                     'part_3': 'pTnos'})

        assert "75 ng of pPE8" in str(response)
        # with more than one part of the same type
        response = client.post(url, {'part_1': 'pP2A11',
                                     'part_2': 'GB0365',
                                     'part_3': 'GB0653',
                                     'part_4': 'GB0655',
                                     'part_5': 'pT35S',
                                     'vector': 'pDGB1_alpha1'})
        assert "75 ng of GB0653" in str(response)

    def test_mantras_bug(self):
        'it test that the protocol file is generated'
        client = Client()
        client.login(username='admin', password='password')
        url = reverse('multipartite_view_add')
        response = client.get(url)

        assert response.status_code == 200

        response = client.post(url, {'Other': 'GB_UD_186',
                                     'Other.2': 'GB_UD_188',
                                     'Vector': 'pDGB1_alpha1',
                                     'category': 'free',
                                     'name': 'aa',
                                     'description': '',
                                     'reference': 'aa',
                                     'order': 'Other:Other.2'})


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
        assert """<select id="id_TER" maxlength="100" name="TER">""" in str(response)
        assert """<option value="pDGB1_alpha1R">pDGB1_alpha""" in str(response)
        client = Client()
        url = reverse('multipartite_view', kwargs={'multi_type': 'basic'})
        response = client.post(url, {"PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTnos',
                                     'Vector': 'pDGB1_alpha1'})

        # print response
        assert 'error' not in response
        assert response.status_code == 200

        client = Client()
        url = reverse('multipartite_view_genbank',
                      kwargs={'multi_type': 'basic'})
        response = client.post(url, {"PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTnos',
                                     'Vector': 'pDGB1_alpha1'})

        assert "LOCUS" in str(response)
        client = Client()
        url = reverse('multipartite_view',
                      kwargs={'multi_type': 'basic'})
        response = client.post(url, {"PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTno'})

        err1 = """<ul class="errorlist"><li>This field is required.</li></ul"""
        assert err1 in str(response)
        err2 = """<ul class="errorlist"><li>This feature does not exist in"""
        assert err2 in str(response)

        # forward vector
        url = reverse('multipartite_view_genbank',
                      kwargs={'multi_type': 'basic'})
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
        url = reverse('multipartite_view_genbank',
                      kwargs={'multi_type': 'basic'})
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

        response = client.post(url, {'assembled_seq': 'aaa',
                                     'multi_type': 'basic',
                                     "PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTnos',
                                     'Vector': 'pDGB1_alpha1'})
        assert "75 ng of pPE8" in str(response)

    def test_genbank_view(self):
        'it test that the protocol file is generated'
        client = Client()
        url = reverse('multipartite_view_genbank', kwargs={'multi_type':
                                                           'basic'})
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'assembled_seq': 'aaa',
                                     'multi_type': 'basic',
                                     "PROM+UTR+ATG": 'pPE8',
                                     "CDS": 'pANT1',
                                     "TER": 'pTnos',
                                     'Vector': 'pDGB1_alpha1'})
        assert 'LOCUS' in str(response)


class BipartiteViewTest(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_bipartite(self):
        client = Client()
        # do initial
        url = reverse('bipartite_view')
        response = client.get(url)
        assert """<option value="GB0125">GB0125 - pEGB 35S:Rosea:Tnos</option>""" in str(response)

        # do page 1
        url = reverse('bipartite_view', kwargs={'form_num': '1'})
        response = client.post(url, {'part_1': 'GB0125'})
        assert 'readonly="True"' in str(response)
        assert 'value="GB0125"' in str(response)
        assert """<p><label for="id_part_2">Part 2:</label>""" in str(response)

        # do page 2
        url = reverse('bipartite_view', kwargs={'form_num': '2'})
        response = client.post(url, {'part_1': 'GB0125', 'part_2': 'GB0126'})
        assert 'value="GB0126"' in str(response)
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
                                     'Vector': 'pDGB1_alpha1'})

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
                                     'Vector': 'pDGB1_omega1'})
        assert 'LOCUS' in str(response)

    # check bipartite_view_protocol
    def test_protocol_view(self):
        'it test that the protocol file is generated'
        client = Client()
        url = reverse('bipartite_view_protocol')
        response = client.get(url)
        assert response.status_code == 400

        response = client.post(url, {'name': 'kk',
                                     'Description': 'desc',
                                     'Reference': 'ref',
                                     'assembled_seq': 'aaa',
                                     'part_1': 'GB0125',
                                     'part_2': 'GB0126',
                                     'Vector': 'pDGB1_omega1'})

        assert 'Bipartite Assembly Protocol' in str(response)

    # check bipartite_view_add
    def test_add_view(self):
        'it test that the protocol file is generated'
        client = Client()
        client.login(username='admin', password='password')
        url = reverse('bipartite_view_add')
        response = client.get(url)

        assert response.status_code == 200

        response = client.post(url, {'assembled_seq': 'aaa',
                                     'part_1': 'GB0125',
                                     'part_2': 'GB0126',
                                     'Vector': 'pDGB1_omega1',
                                     'name': 'aa',
                                     'description': '',
                                     'reference': 'aa'})
        assert response.status_code == 302


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
        # print str(response)
        assert """<ul class="errorlist"><li>The provided s""" in str(response)

        # not add a sequence
        response = client.post(url, {'seq': '',
                                     'category': '12 (NT)'})
        assert """<ul class="errorlist"><li>Fasta or genbank File Required</li></ul>""" in str(response)

        # add category, prefix and suffix

        response = client.post(url, {'seq': open(gb_path),
                                     'prefix': 'ggac', 'suffix': 'cgtc',
                                     'category': '17-21 (TER)'})
        assert """<ul class="errorlist"><li>Can not use category and prefix/suffix simoultaneously</li></ul>"""in str(response)

        # add category and suffix
        response = client.post(url, {'seq': open(gb_path),
                                     'prefix': '', 'suffix': 'cgtc',
                                     'category': '17-21 (TER)'})
        assert """<ul class="errorlist"><li>Can not use category and prefix/suffix simoultaneously</li></ul>"""in str(response)

        # add suffix
        response = client.post(url, {'seq': open(gb_path),
                                     'prefix': '', 'suffix': 'cgtc',
                                     'category': ''})
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

        # domesticate with prefix and suffix
        response = client.post(url, {'seq': open(gb_path),
                                     'suffix': 'ACCT', 'prefix': 'TTCC'})
        assert  "<p>Prefix:TTCC</p>" in str(response)

        residues = str(SeqIO.read(open(gb_path), format='gb').seq)
        response = client.post(url, {'residues': residues,
                                     'category': '13-14-15-16 (CDS)'})
        assert 'The provided seq must start with start' in str(response)

    def test_genbank_view(self):
        'it test that the genbank file is generated'
        client = Client()
        url = reverse('domestication_view_genbank')
        response = client.get(url)
        assert response.status_code == 400
        response = client.post(url, {'seq': 'gagaggggggggagagagattcccctctccccccccccccccccccccccccccccccccccccctttgacctcgaaacgccccc',
                                     'prefix': 'ggag',
                                     'suffix': 'aatg',
                                     'category': '01-02-03-11-12 (PROM+UTR+ATG)',
                                     'seq_name': 'test',
                                     'with_intron': '0'})
        assert 'LOCUS' in str(response)

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
                                     'category': '01-02-03-11-12 (PROM+UTR+ATG)',
                                     'seq_name': 'test',
                                     'with_intron': '0'})
        assert "Oligo forward: GCGCCGTCTCGCTCGGGAGGAGAGGGGGGGGAGAGAGAT" in str(response)
