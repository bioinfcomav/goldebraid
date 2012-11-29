from django.test import TestCase
from goldenbraid.models import Db

FIXTURES_TO_LOAD = ['db', 'cv', 'cvterm', 'dbxref', 'feature']
#                    'feature_cvterm', 'featureprop',
#                    'feature_dbxref',
#                    ]

FIXTURES_TO_LOAD = ['gb_%s.json' % table for table in FIXTURES_TO_LOAD]


class TestFixtures(TestCase):
    'It tests that we can load the fixtures'
    fixtures = FIXTURES_TO_LOAD

    def test_initial(self):
        'It loads the fixtures.'
        assert Db.objects.all()
