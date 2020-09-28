'''
Created on 2015 eka. 1

@author: peio
'''
from django.test.testcases import TestCase
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD
from django.urls import reverse
from django.test.client import Client


class ApiViewTest(TestCase):
    fixtures = ['auth.json'] + FIXTURES_TO_LOAD

    def test_keyfeat(self):
        client = Client()
        url = reverse('api_feature_key_elements')

        response = client.get('api_feature_key_elements',
                              data={'features': ['GB0130']})
        print (str(response))
