from django.test import TestCase
from goldenbraid.models import Db

FIXTURES_TO_LOAD = ['initial_data.json']


class TestFixtures(TestCase):
    'It tests that we can load the fixtures'
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_initial(self):
        'It loads the fixtures.'
        assert Db.objects.all()
