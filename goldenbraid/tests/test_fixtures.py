from django.test import TestCase
from goldenbraid.models import Db
from goldenbraid.settings import DB

FIXTURES_TO_LOAD = ['db', 'cv', 'cvterm', 'dbxref', 'feature', 'featureprop']
#                    'feature_cvterm', ,
#                    'feature_dbxref',
#                    ]

FIXTURES_TO_LOAD = ['gb_%s.json' % table for table in FIXTURES_TO_LOAD]


class TestFixtures(TestCase):
    'It tests that we can load the fixtures'
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_initial(self):
        'It loads the fixtures.'
        assert Db.objects.all()
