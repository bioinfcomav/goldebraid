# Copyright 2013 Diego Orzaez, Univ.Politecnica Valencia, Consejo Superior de
# Investigaciones Cientificas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
'''
tests the goldenbraid models
'''
import os
from django.test import TestCase
from django.db.utils import IntegrityError
from django.core.files import File
from django.conf import settings as proj_settings
from django.db import transaction

import goldenbraid
from goldenbraid.utils import has_rec_sites
from goldenbraid.models import (Db, Dbxref, Cv, Cvterm, Feature, Featureprop,
                                Contact, Stock, Stockcollection, Count,
                                FeaturePerm, FeatureRelationship, Experiment,
                                ExperimentPropExcel)
from goldenbraid.tags import (ENZYME_IN_TYPE_NAME, ENZYME_OUT_TYPE_NAME,
                              VECTOR_TYPE_NAME, RESISTANCE_TYPE_NAME,
                              DERIVES_FROM)
from goldenbraid.tests.test_fixtures import FIXTURES_TO_LOAD4

TEST_DATA = os.path.join(os.path.split(goldenbraid.__path__[0])[0],
                         'goldenbraid', 'tests', 'data')


class FeatureTestModels(TestCase):
    fixtures = FIXTURES_TO_LOAD4

    def test_create(self):
        'can we create a feature?'
        gb_file = File(open(os.path.join(TEST_DATA, 'pAn11_uniq.gb')))
        db = Db.objects.get(name='goldenbraid')
        part1_dbxref = Dbxref.objects.create(db=db, accession='part11')
        cv = Cv.objects.get(name='goldenbraid')

        vector_cvterm = Cvterm.objects.get(cv=cv, name=VECTOR_TYPE_NAME)

        # create a feature that is a vector
        vector_feat = Feature.objects.create(uniquename='vector1',
                                             name='vector1',
                                             type=vector_cvterm,
                                             residues='ATTTAGGCTC',
                                             dbxref=part1_dbxref,
                                             prefix='ATCT',
                                             suffix='tttt',
                                             genbank_file=gb_file)
        selected_vec_feat = Feature.objects.get(uniquename='vector1')
        assert selected_vec_feat.name == 'vector1'

        # try to add again the same feat
        try:
            with transaction.atomic():
                vector_feat = Feature.objects.create(uniquename='vector1',
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
        enzyme_out_cvterm = Cvterm.objects.get(cv=cv,
                                               name=ENZYME_OUT_TYPE_NAME)

        Featureprop.objects.create(feature=vector_feat, type=enzyme_out_cvterm,
                                   value='BSA1', rank=0)
        Featureprop.objects.create(feature=vector_feat, type=enzyme_out_cvterm,
                                   value='BSA2', rank=1)
        assert vector_feat.enzyme_out == ['BSA1', 'BSA2']

        resistance_cvterm = Cvterm.objects.get(cv=cv,
                                               name=RESISTANCE_TYPE_NAME)

        Featureprop.objects.create(feature=vector_feat,
                                   type=resistance_cvterm,
                                   value='Ampicillin', rank=0)
        assert vector_feat.resistance == ['Ampicillin']

        # create a second feature
        part2_dbxref = Dbxref.objects.create(db=db, accession='part2')

        promoter_cvterm = Cvterm.objects.create(cv=cv, name='promoter',
                                                definition='promoter')
        feat = Feature.objects.create(uniquename='part1', name='part1',
                                      type=promoter_cvterm,
                                      residues='ACTC',
                                      dbxref=part1_dbxref,
                                      vector=vector_feat,
                                      prefix='ATCT',
                                      suffix='tttt')
        os.remove(os.path.join(proj_settings.MEDIA_ROOT,
                               vector_feat.genbank_file.name))

        assert feat.enzyme_in is None
        assert feat.enzyme_out == ['BSA1', 'BSA2']
        assert feat.resistance == ['Ampicillin']

        # create a second vector
        part3_dbxref = Dbxref.objects.create(db=db, accession='part3')

        vector2_feat = Feature.objects.create(uniquename='vector2',
                                              name='vector2', type=vector_cvterm,
                                              residues='ATTCATTAGGCTC',
                                              dbxref=part3_dbxref,
                                              prefix='ATCT', suffix='tttt')

        enzyme_in_cvterm = Cvterm.objects.get(cv=cv, name=ENZYME_IN_TYPE_NAME)

        Featureprop.objects.create(feature=vector2_feat, type=enzyme_in_cvterm,
                                   value='BSA1', rank=0)
        Featureprop.objects.create(feature=vector2_feat, rank=0,
                                   type=enzyme_out_cvterm, value='BSA2')

        Featureprop.objects.create(feature=vector2_feat,
                                   type=resistance_cvterm,
                                   value='Kan', rank=0)
        assert vector2_feat.enzyme_in == ['BSA1']

        # add a stock
        ibmcp_contact = Contact.objects.create(name='pepito',
                                               email='pepito@ibmcp.org')

        ibmcp_collection = Stockcollection.objects.create(uniquename='ibmcp',
                                                          contact=ibmcp_contact, name='ibmcp')

        ibmcp_vec_stock = Stock.objects.create(name='vector1_stock_ibmcp',
                                               uniquename='stock_vector1_ibmcp',
                                               stockcollection=ibmcp_collection,
                                               description='vector1_stock in ibmcp description')

        # a second stock in the ibmcp stockcollection
        ibmcp_part_stock = Stock.objects.create(name='part1_stock_ibmcp',
                                                uniquename='stock_part1_ibmcp',
                                                stockcollection=ibmcp_collection,
                                                description='part1_stock in ibmcp description')

        # another collection
        comav_contact = Contact.objects.create(name='pepito',
                                               email='pepito@comav.org')

        comav_collection = Stockcollection.objects.create(
                                            contact=comav_contact, name='comav',
                                            uniquename='comav')

        comav_vec_stock = Stock.objects.create(name='vector1_stock_comav',
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

    def test_counter(self):
        count = Count.objects.create(name='assembled_seq', value=1)

        assert count.next == '1'
        assert count.next == '2'
        assert count.next == '3'
        assert count.value == 4

    def test_featureperm(self):

        feature = Feature.objects.get(uniquename='pAn11')
        assert feature.is_public
        featureperm = FeaturePerm.objects.get(feature=feature)
        featureperm.is_public = False
        featureperm.save()
        assert not feature.is_public

    def test_featurelevel(self):
        feature = Feature.objects.get(uniquename='pAn11')
        assert feature.level == '0'
        f1 = Feature.objects.get(uniquename='pEGBMyb_rev1')
        assert f1.level == '1-alpha'

    def test_featureversion(self):
        feature = Feature.objects.get(uniquename='pAn11')
        assert feature.gb_version == '3.0'

    def test_moclo_compatible(self):
        # alpha and omega  vectors have recsites for this enzymes
        feature = Feature.objects.get(uniquename='pDGB1_alpha1')
        assert has_rec_sites(feature.residues, enzymes=('BpiI', 'BsaI'))
        feature = Feature.objects.get(uniquename='pUPD')
        assert has_rec_sites(feature.residues, enzymes=('BpiI', 'BsaI'))

        feature = Feature.objects.get(uniquename='pAn11')
        assert feature.moclo_compatible
        feature = Feature.objects.get(uniquename='pDGB1_alpha1')
        assert feature.moclo_compatible == 'not_evaluable'

        f1 = Feature.objects.get(uniquename='pDelila')
        assert not f1.moclo_compatible

    def test_gb_category(self):
        feature = Feature.objects.get(uniquename='pAn11')
        assert feature.gb_category_sections == "(B3-B4-B5)"
        feature = Feature.objects.get(uniquename='GB0125')
        assert feature.gb_category_sections is None

        feature = Feature.objects.get(uniquename='pUPD')
        assert feature.gb_category_sections is None

    def test_feature_relationship(self):

        cvterm = Cvterm.objects.get(name=DERIVES_FROM)
        f1 = Feature.objects.get(feature_id=47)
        f2 = Feature.objects.get(uniquename="GB0365")
        fet_rel = FeatureRelationship.objects.create(type=cvterm, object=f1,
                                                     subject=f2)
        assert fet_rel.type.name == DERIVES_FROM
        assert f1.children[0].uniquename == "GB0365"

    def test_feature_images(self):
        f1 = Feature.objects.get(uniquename='GB0125')
        assert f1.experiment_images == [(u'/api/excel_graph/1', u'excel'),
                                        (u'/api/excel_graph/2', u'columns')]
        f1 = Feature.objects.get(uniquename='GB0129')
        print f1.experiment_images
        assert f1.experiment_images == [(u'/api/excel_graph/8', u'columns'),
                                        (u'/api/excel_graph/9', u'columns')]

    def test_sbol_image(self):
        f1 = Feature.objects.get(uniquename='GB_UA_15')
        assert f1.sbol_images == ['prom_5utr_ntag.png', 'cds.png',
                                  '3utr_term.png', '3utr_term.png', 'cds.png',
                                  'prom_5utr_ntag.png']


class ExperimentTests(TestCase):
    fixtures = FIXTURES_TO_LOAD4

    def test_experiment_excels(self):
        exp = Experiment.objects.get(uniquename='GB_EXP_2B')
        excel_exp = ExperimentPropExcel.objects.get(experiment=exp)
        fhand = open(os.path.join(TEST_DATA, 'scatter.xlsx'))
        excel_exp.excel = File(fhand)
        excel_exp.save()
        excel_exp.excel.close()
        fhand.close()

        assert excel_exp.drawed_image[1] == 'image/png'

        # image url from excel
        assert exp.image_urls == [(u'/api/excel_graph/1', u'excel',)]
        # image_url_from_image
        exp = Experiment.objects.get(uniquename='GB_EXP_2A')
        assert not exp.image_urls
        exp = Experiment.objects.get(uniquename='GB_EXP_29')
        assert exp.image_urls == [('/media/result_files/hayedo.jpg', 'Image1'),
                                  ('/media/result_files/escocia.jpg',
                                   'Image2')]

        # print open(excelExp.excel).read()
