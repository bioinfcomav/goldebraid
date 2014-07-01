
from django.test.testcases import TestCase

from gb_genome_domestication.tests.test_fixtures import FIXTURES_TO_LOAD
from gb_genome_domestication.models import Feature
from gb_genome_domestication.settings import DB, DB_URLPREFIX


class FeatureTestModels(TestCase):
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_featureperm(self):
        feature = Feature.objects.using(DB).get(uniquename='AT1G01020.1')
        assert feature.primary_dbxref.url == DB_URLPREFIX + 'AT1G01020.1'
        assert feature.species == 'arabidopsis_thaliana'

        assert feature.num_rec_sites == 2

        dbxrefs = list(feature.secondary_dbxrefs)
        dbxref = dbxrefs[0]
        assert dbxref.db.name == 'plant_ensembl_arabidoposis'
        assert dbxref.accession == feature.uniquename



