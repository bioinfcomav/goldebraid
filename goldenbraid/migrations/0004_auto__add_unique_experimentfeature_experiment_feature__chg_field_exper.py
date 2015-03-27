# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding unique constraint on 'ExperimentFeature', fields ['experiment', 'feature']
        db.create_unique(u'experimentfeature', ['experiment_id', 'feature_id'])


        # Changing field 'ExperimentPerm.experiment'
        db.alter_column(u'experimentperm', 'experiment_id', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['goldenbraid.Experiment'], unique=True, primary_key=True))
        # Adding unique constraint on 'ExperimentKeySubFeature', fields ['experiment', 'feature']
        db.create_unique(u'experimentkeysubfeature', ['experiment_id', 'feature_id'])


    def backwards(self, orm):
        # Removing unique constraint on 'ExperimentKeySubFeature', fields ['experiment', 'feature']
        db.delete_unique(u'experimentkeysubfeature', ['experiment_id', 'feature_id'])

        # Removing unique constraint on 'ExperimentFeature', fields ['experiment', 'feature']
        db.delete_unique(u'experimentfeature', ['experiment_id', 'feature_id'])


        # Changing field 'ExperimentPerm.experiment'
        db.alter_column(u'experimentperm', 'experiment_id', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['goldenbraid.Feature'], unique=True, primary_key=True))

    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Group']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Permission']"}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'goldenbraid.contact': {
            'Meta': {'object_name': 'Contact', 'db_table': "u'contact'"},
            'contact_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'email': ('django.db.models.fields.EmailField', [], {'unique': 'True', 'max_length': '75'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '255'})
        },
        u'goldenbraid.count': {
            'Meta': {'object_name': 'Count', 'db_table': "u'count'"},
            'count_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'value': ('django.db.models.fields.IntegerField', [], {})
        },
        u'goldenbraid.cv': {
            'Meta': {'object_name': 'Cv', 'db_table': "u'cv'"},
            'cv_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'definition': ('django.db.models.fields.TextField', [], {}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '255'})
        },
        u'goldenbraid.cvterm': {
            'Meta': {'object_name': 'Cvterm', 'db_table': "u'cvterm'"},
            'cv': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Cv']"}),
            'cvterm_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'definition': ('django.db.models.fields.TextField', [], {}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '1024'})
        },
        u'goldenbraid.db': {
            'Meta': {'object_name': 'Db', 'db_table': "u'db'"},
            'db_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '255'}),
            'url': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'urlprefix': ('django.db.models.fields.CharField', [], {'max_length': '255'})
        },
        u'goldenbraid.dbxref': {
            'Meta': {'object_name': 'Dbxref', 'db_table': "u'dbxref'"},
            'accession': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'db': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Db']"}),
            'dbxref_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'goldenbraid.experiment': {
            'Meta': {'object_name': 'Experiment', 'db_table': "'experiment'"},
            'chasis_1': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'chasis_2': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'dbxref': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Dbxref']"}),
            'description': ('django.db.models.fields.TextField', [], {'max_length': '255'}),
            'experiment_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'timecreation': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Cvterm']"}),
            'uniquename': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '255'})
        },
        u'goldenbraid.experimentfeature': {
            'Meta': {'unique_together': "(('experiment', 'feature'),)", 'object_name': 'ExperimentFeature', 'db_table': "u'experimentfeature'"},
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Experiment']"}),
            'experiment_feature_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'feature': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Feature']"}),
            'role': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Cvterm']"})
        },
        u'goldenbraid.experimentkeysubfeature': {
            'Meta': {'unique_together': "(('experiment', 'feature'),)", 'object_name': 'ExperimentKeySubFeature', 'db_table': "u'experimentkeysubfeature'"},
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Experiment']"}),
            'experiment_key_subfeature_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'feature': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Feature']"})
        },
        u'goldenbraid.experimentperm': {
            'Meta': {'object_name': 'ExperimentPerm', 'db_table': "u'experimentperm'"},
            'experiment': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['goldenbraid.Experiment']", 'unique': 'True', 'primary_key': 'True'}),
            'is_public': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        },
        u'goldenbraid.experimentpropexcel': {
            'Meta': {'object_name': 'ExperimentPropExcel', 'db_table': "u'experimentpropexcel'"},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Experiment']"}),
            'experiment_prop_excel_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'image': ('django.db.models.fields.files.FileField', [], {'max_length': '100'})
        },
        u'goldenbraid.experimentpropimage': {
            'Meta': {'object_name': 'ExperimentPropImage', 'db_table': "u'experimentpropimage'"},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Experiment']"}),
            'experiment_prop_image_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'image': ('django.db.models.fields.files.ImageField', [], {'max_length': '100'})
        },
        u'goldenbraid.experimentpropnumeric': {
            'Meta': {'object_name': 'ExperimentPropNumeric', 'db_table': "u'experimentpropnumeric'"},
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Experiment']"}),
            'experiment_prop_numeric_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Cvterm']"}),
            'value': ('django.db.models.fields.FloatField', [], {})
        },
        u'goldenbraid.experimentproptext': {
            'Meta': {'object_name': 'ExperimentPropText', 'db_table': "u'experimentproptext'"},
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Experiment']"}),
            'experiment_prop_text_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'value': ('django.db.models.fields.TextField', [], {'max_length': '255'})
        },
        u'goldenbraid.feature': {
            'Meta': {'object_name': 'Feature', 'db_table': "u'feature'"},
            'dbxref': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Dbxref']"}),
            'feature_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'genbank_file': ('django.db.models.fields.files.FileField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'prefix': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'residues': ('django.db.models.fields.TextField', [], {}),
            'suffix': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'timecreation': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'timelastmodified': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Cvterm']"}),
            'uniquename': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '255'}),
            'vector': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Feature']", 'null': 'True'})
        },
        u'goldenbraid.featureperm': {
            'Meta': {'object_name': 'FeaturePerm', 'db_table': "u'featureperm'"},
            'feature': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['goldenbraid.Feature']", 'unique': 'True', 'primary_key': 'True'}),
            'is_public': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        },
        u'goldenbraid.featureprop': {
            'Meta': {'object_name': 'Featureprop', 'db_table': "u'featureprop'"},
            'feature': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Feature']"}),
            'featureprop_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'rank': ('django.db.models.fields.IntegerField', [], {}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Cvterm']"}),
            'value': ('django.db.models.fields.TextField', [], {})
        },
        u'goldenbraid.featurerelationship': {
            'Meta': {'object_name': 'FeatureRelationship', 'db_table': "u'feature_relationship'"},
            'featurerelationship_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'object': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'object'", 'to': u"orm['goldenbraid.Feature']"}),
            'subject': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'subject'", 'to': u"orm['goldenbraid.Feature']"}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Cvterm']"})
        },
        u'goldenbraid.stock': {
            'Meta': {'object_name': 'Stock', 'db_table': "u'stock'"},
            'description': ('django.db.models.fields.TextField', [], {}),
            'feature': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'stocks'", 'null': 'True', 'to': u"orm['goldenbraid.Feature']"}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'stock_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'stockcollection': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Stockcollection']"}),
            'uniquename': ('django.db.models.fields.TextField', [], {})
        },
        u'goldenbraid.stockcollection': {
            'Meta': {'object_name': 'Stockcollection', 'db_table': "u'stockcollection'"},
            'contact': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['goldenbraid.Contact']"}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'stockcollection_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'uniquename': ('django.db.models.fields.TextField', [], {})
        }
    }

    complete_apps = ['goldenbraid']