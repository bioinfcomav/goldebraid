# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Db'
        db.create_table(u'db', (
            ('db_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=255)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('urlprefix', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('url', self.gf('django.db.models.fields.CharField')(max_length=255)),
        ))
        db.send_create_signal(u'goldenbraid', ['Db'])

        # Adding model 'Cv'
        db.create_table(u'cv', (
            ('cv_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=255)),
            ('definition', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'goldenbraid', ['Cv'])

        # Adding model 'Dbxref'
        db.create_table(u'dbxref', (
            ('dbxref_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('db', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['goldenbraid.Db'])),
            ('accession', self.gf('django.db.models.fields.CharField')(max_length=255)),
        ))
        db.send_create_signal(u'goldenbraid', ['Dbxref'])

        # Adding model 'Cvterm'
        db.create_table(u'cvterm', (
            ('cvterm_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('cv', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['goldenbraid.Cv'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=1024)),
            ('definition', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'goldenbraid', ['Cvterm'])

        # Adding model 'Count'
        db.create_table(u'count', (
            ('count_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('value', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'goldenbraid', ['Count'])

        # Adding model 'Contact'
        db.create_table(u'contact', (
            ('contact_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('email', self.gf('django.db.models.fields.EmailField')(unique=True, max_length=75)),
        ))
        db.send_create_signal(u'goldenbraid', ['Contact'])

        # Adding model 'Stockcollection'
        db.create_table(u'stockcollection', (
            ('stockcollection_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('contact', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['goldenbraid.Contact'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('uniquename', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'goldenbraid', ['Stockcollection'])

        # Adding model 'Stock'
        db.create_table(u'stock', (
            ('stock_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('uniquename', self.gf('django.db.models.fields.TextField')()),
            ('description', self.gf('django.db.models.fields.TextField')()),
            ('stockcollection', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['goldenbraid.Stockcollection'])),
            ('feature', self.gf('django.db.models.fields.related.ForeignKey')(related_name='stocks', null=True, to=orm['goldenbraid.Feature'])),
        ))
        db.send_create_signal(u'goldenbraid', ['Stock'])

        # Adding model 'Feature'
        db.create_table(u'feature', (
            ('feature_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('uniquename', self.gf('django.db.models.fields.CharField')(unique=True, max_length=255)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['goldenbraid.Cvterm'])),
            ('residues', self.gf('django.db.models.fields.TextField')()),
            ('dbxref', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['goldenbraid.Dbxref'])),
            ('vector', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['goldenbraid.Feature'], null=True)),
            ('prefix', self.gf('django.db.models.fields.CharField')(max_length=4)),
            ('suffix', self.gf('django.db.models.fields.CharField')(max_length=4)),
            ('genbank_file', self.gf('django.db.models.fields.files.FileField')(max_length=100)),
            ('timecreation', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('timelastmodified', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'goldenbraid', ['Feature'])

        # Adding model 'FeaturePerm'
        db.create_table(u'featureperm', (
            ('feature', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['goldenbraid.Feature'], unique=True, primary_key=True)),
            ('owner', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('is_public', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'goldenbraid', ['FeaturePerm'])

        # Adding model 'Featureprop'
        db.create_table(u'featureprop', (
            ('featureprop_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('feature', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['goldenbraid.Feature'])),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['goldenbraid.Cvterm'])),
            ('value', self.gf('django.db.models.fields.TextField')()),
            ('rank', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'goldenbraid', ['Featureprop'])

        # Adding model 'FeatureRelationship'
        db.create_table(u'feature_relationship', (
            ('featurerelationship_id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['goldenbraid.Cvterm'])),
            ('subject', self.gf('django.db.models.fields.related.ForeignKey')(related_name='subject', to=orm['goldenbraid.Feature'])),
            ('object', self.gf('django.db.models.fields.related.ForeignKey')(related_name='object', to=orm['goldenbraid.Feature'])),
        ))
        db.send_create_signal(u'goldenbraid', ['FeatureRelationship'])


    def backwards(self, orm):
        # Deleting model 'Db'
        db.delete_table(u'db')

        # Deleting model 'Cv'
        db.delete_table(u'cv')

        # Deleting model 'Dbxref'
        db.delete_table(u'dbxref')

        # Deleting model 'Cvterm'
        db.delete_table(u'cvterm')

        # Deleting model 'Count'
        db.delete_table(u'count')

        # Deleting model 'Contact'
        db.delete_table(u'contact')

        # Deleting model 'Stockcollection'
        db.delete_table(u'stockcollection')

        # Deleting model 'Stock'
        db.delete_table(u'stock')

        # Deleting model 'Feature'
        db.delete_table(u'feature')

        # Deleting model 'FeaturePerm'
        db.delete_table(u'featureperm')

        # Deleting model 'Featureprop'
        db.delete_table(u'featureprop')

        # Deleting model 'FeatureRelationship'
        db.delete_table(u'feature_relationship')


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