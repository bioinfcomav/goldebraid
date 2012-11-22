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
        ibmcp_contact = Contact(name='pepito', email='pepito@ibmcp.org')
        ibmcp_contact.save()
        ibmcp_collection = Stockcollection(contact=ibmcp_contact, name='ibmcp',
                                           uniquename='ibmcp')
        ibmcp_collection.save()

        ibmcp_vec_stock = Stock(name='vector1_stock_ibmcp',
                                uniquename='stock_vector1_ibmcp',
                                stockcollection=ibmcp_collection,
                              description='vector1_stock in ibmcp description')
        ibmcp_vec_stock.save()
        # another collection
        comav_contact = Contact(name='pepito', email='pepito@comav.org')
        comav_contact.save()
        comav_collection = Stockcollection(contact=comav_contact, name='comav',
                                           uniquename='comav')
        comav_collection.save()
        comav_vec_stock = Stock(name='vector1_stock_comav',
                                uniquename='stock_vector1_comav',
                                stockcollection=comav_collection,
                              description='vector1_stock in comav description')
        comav_vec_stock.save()
        assert ibmcp_collection.contact.name == 'pepito'

        # There are various stocks of the vector feature
        selected_vec_feat.stocks.add(ibmcp_vec_stock)
        selected_vec_feat.stocks.add(comav_vec_stock)
        stocks = selected_vec_feat.stocks.all()
        stock_names = [stock.name for stock in stocks]
        assert 'vector1_stock_ibmcp' in stock_names
        assert 'vector1_stock_comav' in stock_names
