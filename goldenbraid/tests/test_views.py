import os.path

from django.test import TestCase, Client
from django.core.urlresolvers import reverse
from django.core.files.uploadedfile import SimpleUploadedFile

import goldenbraid
from goldenbraid.views.feature_views import FeatureForm
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD
from goldenbraid.models import Feature
from goldenbraid.settings import DB
from goldenbraid.tags import VECTOR_TYPE_NAME


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
        test_data = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                                 'goldenbraid', 'tests', 'data')
        # test of the form
        gb_path = os.path.join(test_data, 'pAn11.gb')
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





