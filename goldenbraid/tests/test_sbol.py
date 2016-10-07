import os
from django.test import TestCase
from Bio.SeqIO import read

import goldenbraid
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD, FIXTURES_TO_LOAD5
from goldenbraid.sbol import convert_to_sbol
# from goldenbraid.models import Feature

TEST_DATA = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                         'goldenbraid', 'tests', 'data')


class TestSBOL(TestCase):
    'It tests that we can load the fixtures'
    # fixtures = FIXTURES_TO_LOAD5
    multi_db = True

    def test_sbol(self):
        # feat = os.path.join(TEST_DATA, 'pEGBAn11.gb')
        # print feat.genbank_file
        seq = read(os.path.join(TEST_DATA, 'pEGBAn11.gb'), 'gb')
        fhand = open('/tmp/sbol.xml', 'w')
        fhand.write(convert_to_sbol(seq))
        fhand.close()
