'''
tests the goldenbraid models
'''

from django.test import TestCase
from goldenbraid import settings
from goldenbraid.models import (Db, Dbxref, Cv, Cvterm, Feature, Featureprop,
                                Contact, Stock, Stockcollection)

DB = settings.DB


class FeatureTest(TestCase):

    def test_create(self):
        'can we create a feature?'
        db = Db(name='goldebraid', description='testdb', urlprefix='/',
                url='localhost/')
        db.save()
        part1_dbxref = Dbxref(db=db, accession='part11')
        part1_dbxref.save()
        cv = Cv(name='goldenbraid', definition='goldenbraid control voc')
        cv.save()
        vector_cvterm = Cvterm(cv=cv, name='vector', definition='vector type')
        vector_cvterm.save()
        vector_feat = Feature(uniquename='vector1', name='vector1',
                              type=vector_cvterm, residues='ATTTAGGCTC',
                              dbxref=part1_dbxref)

        vector_feat.save()
        selected_vec_feat = Feature.objects.using(DB).get(uniquename='vector1')

        assert selected_vec_feat.name == 'vector1'

        enzyme_cvterm = Cvterm(cv=cv, name='Enzyme', definition='enzyme')
        enzyme_cvterm.save()

        feat_prop = Featureprop(feature=vector_feat, type=enzyme_cvterm,
                                value='BSA1')
        feat_prop.save()

        assert vector_feat.props['Enzyme'] == 'BSA1'
        part2_dbxref = Dbxref(db=db, accession='part2')
        part2_dbxref.save()
        feat = Feature(uniquename='part1', name='part1',
                              type=vector_cvterm, residues='ACTC',
                              dbxref=part1_dbxref, vector=vector_feat)
        feat.save()

        assert feat.vector.props['Enzyme'] == 'BSA1'

        # add a stock
        contact = Contact(name='pepito', email='pepito@casa.net')
        contact.save()
        stock = Stock(name='vector1_stock', uniquename='stock1_vector',
              description='vector1_stock description')
        stock.save()
        collection = Stockcollection(contact=contact, name='ibmcp',
                                     uniquename='ibmcp')
        collection.save()

        collection.stocks.add(stock)

        assert collection.contact.name == 'pepito'

        feat.stocks.add(stock)

        assert 'pepito' == stock.in_stockcollections.all()[0].contact.name



