'''
tests the goldenbraid models
'''

from django.test import TestCase
from goldenbraid import settings
from goldenbraid.models import (Db, Dbxref, Cv, Cvterm, Feature, Featureprop,
                                Contact, Stock, Stockcollection)
from goldenbraid.tags import (ENZYME_IN_TYPE_NAME, ENZYME_OUT_TYPE_NAME,
                              VECTOR_TYPE_NAME, RESISTANCE_TYPE_NAME)
from django.db.utils import IntegrityError

DB = settings.DB


class FeatureTestModels(TestCase):

    def test_create(self):
        'can we create a feature?'

        db = Db.objects.using(DB).create(name='goldenbraid',
                                         description='testdb', urlprefix='/',
                                         url='localhost/')

        part1_dbxref = Dbxref.objects.using(DB).create(db=db, accession='part11')

        cv = Cv.objects.using(DB).create(name='goldenbraid',
                                         definition='goldenbraid control voc')

        vector_cvterm = Cvterm.objects.using(DB).create(cv=cv, name=VECTOR_TYPE_NAME,
                                                        definition='vector type')


        # create a feature that is a vector
        vector_feat = Feature.objects.using(DB).create(uniquename='vector1',
                                                       name='vector1',
                                                       type=vector_cvterm,
                                                       residues='ATTTAGGCTC',
                                                       dbxref=part1_dbxref,
                                                       prefix='ATCT',
                                                       suffix='tttt')

        selected_vec_feat = Feature.objects.using(DB).get(uniquename='vector1')
        assert selected_vec_feat.name == 'vector1'

        # try to add again the same feat
        try:
            vector_feat = Feature.objects.using(DB).create(uniquename='vector1',
                                                       name='vector1',
                                                       type=vector_cvterm,
                                                       residues='ATTTAGGCTC',
                                                       dbxref=part1_dbxref,
                                                       prefix='ATCT',
                                                       suffix='tttt')
            self.fail()
        except IntegrityError:
            pass

        # add the properties to the feature
        enzyme_out_cvterm = Cvterm.objects.using(DB).create(cv=cv,
                                                     name=ENZYME_OUT_TYPE_NAME,
                                                     definition='enzyme out')

        Featureprop.objects.using(DB).create(feature=vector_feat,
                                                         type=enzyme_out_cvterm,
                                                         value='BSA1', rank=0)
        Featureprop.objects.using(DB).create(feature=vector_feat,
                                                         type=enzyme_out_cvterm,
                                                         value='BSA2', rank=1)
        assert vector_feat.enzyme_out == ['BSA1', 'BSA2']


        resistance_cvterm = Cvterm.objects.using(DB).create(cv=cv,
                                                            name=RESISTANCE_TYPE_NAME,
                                                            definition='resistance')

        Featureprop.objects.using(DB).create(feature=vector_feat,
                                                type=resistance_cvterm,
                                                value='Ampicillin', rank=0)
        assert vector_feat.resistance == ['Ampicillin']

        # create a second feature
        part2_dbxref = Dbxref.objects.using(DB).create(db=db, accession='part2')

        promoter_cvterm = Cvterm.objects.using(DB).create(cv=cv,
                                                     name='promoter',
                                                     definition='promoter')
        feat = Feature.objects.using(DB).create(uniquename='part1',
                                                name='part1',
                                                type=promoter_cvterm,
                                                residues='ACTC',
                                                dbxref=part1_dbxref,
                                                vector=vector_feat,
                                                prefix='ATCT',
                                                suffix='tttt')

        assert feat.enzyme_in is None
        assert feat.enzyme_out == ['BSA1', 'BSA2']
        assert feat.resistance == ['Ampicillin']

        # create a second vector
        part3_dbxref = Dbxref.objects.using(DB).create(db=db, accession='part3')

        vector2_feat = Feature.objects.using(DB).create(uniquename='vector2',
                                                       name='vector2', type=vector_cvterm,
                                                       residues='ATTCATTAGGCTC',
                                                       dbxref=part3_dbxref,
                                                       prefix='ATCT',
                                                       suffix='tttt')

        enzyme_in_cvterm = Cvterm.objects.using(DB).create(cv=cv,
                                                            name=ENZYME_IN_TYPE_NAME,
                                                            definition='enzyme in')

        Featureprop.objects.using(DB).create(feature=vector2_feat,
                                                         type=enzyme_in_cvterm,
                                                         value='BSA1', rank=0)
        Featureprop.objects.using(DB).create(feature=vector2_feat,
                                                         type=enzyme_out_cvterm,
                                                         value='BSA2', rank=0)

        Featureprop.objects.using(DB).create(feature=vector2_feat,
                                                type=resistance_cvterm,
                                                value='Kan', rank=0)
        assert vector2_feat.enzyme_in == ['BSA1']

        # add a stock
        ibmcp_contact = Contact.objects.using(DB).create(name='pepito',
                                                         email='pepito@ibmcp.org')

        ibmcp_collection = Stockcollection.objects.using(DB).create(
                                                contact=ibmcp_contact, name='ibmcp',
                                                uniquename='ibmcp')

        ibmcp_vec_stock = Stock.objects.using(DB).create(name='vector1_stock_ibmcp',
                                uniquename='stock_vector1_ibmcp',
                                stockcollection=ibmcp_collection,
                                description='vector1_stock in ibmcp description')

        # a second stock in the ibmcp stockcollection
        ibmcp_part_stock = Stock.objects.using(DB).create(name='part1_stock_ibmcp',
                                uniquename='stock_part1_ibmcp',
                                stockcollection=ibmcp_collection,
                                description='part1_stock in ibmcp description')

        # another collection
        comav_contact = Contact.objects.using(DB).create(name='pepito',
                                                         email='pepito@comav.org')

        comav_collection = Stockcollection.objects.using(DB).create(
                                            contact=comav_contact, name='comav',
                                            uniquename='comav')

        comav_vec_stock = Stock.objects.using(DB).create(name='vector1_stock_comav',
                                uniquename='stock_vector1_comav',
                                stockcollection=comav_collection,
                              description='vector1_stock in comav description')

        assert ibmcp_collection.contact.name == 'pepito'
        assert comav_collection.contact.name == 'pepito'


        # There are various stocks of the vector feature
        selected_vec_feat.stocks.add(ibmcp_vec_stock)
        selected_vec_feat.stocks.add(comav_vec_stock)
        stocks = selected_vec_feat.stocks.all()
        stock_names = [stock.name for stock in stocks]
        assert 'vector1_stock_ibmcp' in stock_names
        assert 'vector1_stock_comav' in stock_names


        # There are various stocks in a stock collection
        assert ibmcp_vec_stock.stockcollection.name == 'ibmcp'
        assert ibmcp_part_stock.stockcollection.name == 'ibmcp'








