import os.path

from django.test import TestCase, Client
from django.core.urlresolvers import reverse
from django.core.files.uploadedfile import SimpleUploadedFile

import goldenbraid
from goldenbraid.views import FeatureForm
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD
from goldenbraid.models import Feature
from goldenbraid.settings import DB
from goldenbraid.tags import VECTOR_TYPE_NAME, DESCRIPTION_TYPE_NAME


class FeatureTestViews(TestCase):
    fixtures = FIXTURES_TO_LOAD

    def test_add_feature_form(self):
        test_data = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                                 'goldenbraid', 'tests', 'data')
        # test of the form
        gb_path = os.path.join(test_data, 'pAn11.gb')
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_in': 'vector1_enz_in',
                     'enzyme_out': 'vector1_enz_out'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertTrue(form.is_valid())

        # test of the form with blanck values
        gb_path = os.path.join(test_data, 'pAn11.gb')
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_out': 'vector1_enz_out',
                     'enzyme_in': 'vector1_enz_in'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertTrue(form.is_valid())

        # test of the form with wrong type
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': 'vecto', 'enzyme_out': 'vector1_enz_out'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertFalse(form.is_valid())
        assert form._errors.get('type')

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
        assert form._errors.get('vector')

        # enzyme_in not added wrong
        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_out': 'vector1_enz_out'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertFalse(form.is_valid())
        enz_in = 'enzyme_in'
        assert 'A vector must have a enzyme i' in str(form._errors.get(enz_in))

        post_dict = {'uniquename': 'vector1', 'name': 'vector1',
                     'type': VECTOR_TYPE_NAME, 'enzyme_out': 'vector1_enz_out',
                     'enzyme_in': 'no_exist'}
        uploaded_fhand = open(gb_path)
        file_dict = {'gbfile': SimpleUploadedFile(uploaded_fhand.name,
                                                  uploaded_fhand.read())}
        form = FeatureForm(post_dict, file_dict)
        self.assertFalse(form.is_valid())
        enz_in = 'enzyme_in'

        assert 'iven enzyme in does not exist' in str(form._errors.get(enz_in))


    def test_add_feature_view(self):
        # test of the form page
        test_data = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                                 'goldenbraid', 'tests', 'data')
        # test of the form
        gb_path = os.path.join(test_data, 'pAn11.gb')
        client = Client()
        url = reverse('add_feature')
        response = client.post(url, {'uniquename': 'vector1',
                                     'name': 'vector1',
                                     'type': VECTOR_TYPE_NAME,
                                     'description': 'vector1 desc',
                                     'enzyme_in': 'vector1_enz_in',
                                     'enzyme_out': 'vector1_enz_out',
                                     'gbfile': open(gb_path)})
        feat = Feature.objects.using(DB).get(uniquename='vector1')
        assert feat.name == 'vector1'
        assert feat.props[DESCRIPTION_TYPE_NAME] == 'vector1 desc'
        print response




