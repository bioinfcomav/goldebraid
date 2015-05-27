# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Contact',
            fields=[
                ('contact_id', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=255)),
                ('email', models.EmailField(unique=True, max_length=254)),
            ],
            options={
                'db_table': 'contact',
            },
        ),
        migrations.CreateModel(
            name='Count',
            fields=[
                ('count_id', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=255)),
                ('value', models.IntegerField()),
            ],
            options={
                'db_table': 'count',
            },
        ),
        migrations.CreateModel(
            name='Cv',
            fields=[
                ('cv_id', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=255)),
                ('definition', models.TextField()),
            ],
            options={
                'db_table': 'cv',
            },
        ),
        migrations.CreateModel(
            name='Cvterm',
            fields=[
                ('cvterm_id', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=1024)),
                ('definition', models.TextField()),
                ('cv', models.ForeignKey(to='goldenbraid.Cv')),
            ],
            options={
                'db_table': 'cvterm',
            },
        ),
        migrations.CreateModel(
            name='Db',
            fields=[
                ('db_id', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=255)),
                ('description', models.CharField(max_length=255)),
                ('urlprefix', models.CharField(max_length=255)),
                ('url', models.CharField(max_length=255)),
            ],
            options={
                'db_table': 'db',
            },
        ),
        migrations.CreateModel(
            name='Dbxref',
            fields=[
                ('dbxref_id', models.AutoField(serialize=False, primary_key=True)),
                ('accession', models.CharField(max_length=255)),
                ('db', models.ForeignKey(to='goldenbraid.Db')),
            ],
            options={
                'db_table': 'dbxref',
            },
        ),
        migrations.CreateModel(
            name='Experiment',
            fields=[
                ('experiment_id', models.AutoField(serialize=False, primary_key=True)),
                ('uniquename', models.CharField(unique=True, max_length=255)),
                ('chasis_1', models.CharField(max_length=255)),
                ('chasis_2', models.CharField(max_length=255)),
                ('description', models.TextField(max_length=255)),
                ('timecreation', models.DateTimeField(auto_now_add=True)),
            ],
            options={
                'db_table': 'experiment',
            },
        ),
        migrations.CreateModel(
            name='ExperimentFeature',
            fields=[
                ('experiment_feature_id', models.AutoField(serialize=False, primary_key=True)),
            ],
            options={
                'db_table': 'experimentfeature',
            },
        ),
        migrations.CreateModel(
            name='ExperimentPropExcel',
            fields=[
                ('experiment_prop_excel_id', models.AutoField(serialize=False, primary_key=True)),
                ('image', models.FileField(upload_to=b'result_files')),
                ('description', models.CharField(max_length=255)),
            ],
            options={
                'db_table': 'experimentpropexcel',
            },
        ),
        migrations.CreateModel(
            name='ExperimentPropImage',
            fields=[
                ('experiment_prop_image_id', models.AutoField(serialize=False, primary_key=True)),
                ('image', models.ImageField(upload_to=b'result_files')),
                ('description', models.CharField(max_length=255)),
            ],
            options={
                'db_table': 'experimentpropimage',
            },
        ),
        migrations.CreateModel(
            name='ExperimentPropNumeric',
            fields=[
                ('experiment_prop_numeric_id', models.AutoField(serialize=False, primary_key=True)),
                ('value', models.FloatField()),
            ],
            options={
                'db_table': 'experimentpropnumeric',
            },
        ),
        migrations.CreateModel(
            name='ExperimentPropText',
            fields=[
                ('experiment_prop_text_id', models.AutoField(serialize=False, primary_key=True)),
                ('title', models.CharField(max_length=255)),
                ('value', models.TextField(max_length=255)),
            ],
            options={
                'db_table': 'experimentproptext',
            },
        ),
        migrations.CreateModel(
            name='ExperimentSubFeature',
            fields=[
                ('experiment_key_subfeature_id', models.AutoField(serialize=False, primary_key=True)),
            ],
            options={
                'db_table': 'experimentkeysubfeature',
            },
        ),
        migrations.CreateModel(
            name='Feature',
            fields=[
                ('feature_id', models.AutoField(serialize=False, primary_key=True)),
                ('uniquename', models.CharField(unique=True, max_length=255)),
                ('name', models.CharField(max_length=255)),
                ('residues', models.TextField()),
                ('prefix', models.CharField(max_length=4)),
                ('suffix', models.CharField(max_length=4)),
                ('genbank_file', models.FileField(upload_to=b'genbank_files')),
                ('timecreation', models.DateTimeField(auto_now_add=True)),
                ('timelastmodified', models.DateTimeField(auto_now=True)),
            ],
            options={
                'db_table': 'feature',
            },
        ),
        migrations.CreateModel(
            name='Featureprop',
            fields=[
                ('featureprop_id', models.AutoField(serialize=False, primary_key=True)),
                ('value', models.TextField()),
                ('rank', models.IntegerField()),
            ],
            options={
                'db_table': 'featureprop',
            },
        ),
        migrations.CreateModel(
            name='FeatureRelationship',
            fields=[
                ('featurerelationship_id', models.AutoField(serialize=False, primary_key=True)),
            ],
            options={
                'db_table': 'feature_relationship',
            },
        ),
        migrations.CreateModel(
            name='Stock',
            fields=[
                ('stock_id', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=255)),
                ('uniquename', models.TextField()),
                ('description', models.TextField()),
            ],
            options={
                'db_table': 'stock',
            },
        ),
        migrations.CreateModel(
            name='Stockcollection',
            fields=[
                ('stockcollection_id', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=255)),
                ('uniquename', models.TextField()),
                ('contact', models.ForeignKey(to='goldenbraid.Contact')),
            ],
            options={
                'db_table': 'stockcollection',
            },
        ),
        migrations.CreateModel(
            name='ExperimentPerm',
            fields=[
                ('experiment', models.OneToOneField(primary_key=True, serialize=False, to='goldenbraid.Experiment')),
                ('is_public', models.BooleanField(default=False)),
                ('owner', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'db_table': 'experimentperm',
            },
        ),
        migrations.CreateModel(
            name='FeaturePerm',
            fields=[
                ('feature', models.OneToOneField(primary_key=True, serialize=False, to='goldenbraid.Feature')),
                ('is_public', models.BooleanField(default=False)),
                ('owner', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'db_table': 'featureperm',
            },
        ),
        migrations.AddField(
            model_name='stock',
            name='feature',
            field=models.ForeignKey(related_name='stocks', to='goldenbraid.Feature', null=True),
        ),
        migrations.AddField(
            model_name='stock',
            name='stockcollection',
            field=models.ForeignKey(to='goldenbraid.Stockcollection'),
        ),
        migrations.AddField(
            model_name='featurerelationship',
            name='object',
            field=models.ForeignKey(related_name='object', to='goldenbraid.Feature'),
        ),
        migrations.AddField(
            model_name='featurerelationship',
            name='subject',
            field=models.ForeignKey(related_name='subject', to='goldenbraid.Feature'),
        ),
        migrations.AddField(
            model_name='featurerelationship',
            name='type',
            field=models.ForeignKey(to='goldenbraid.Cvterm'),
        ),
        migrations.AddField(
            model_name='featureprop',
            name='feature',
            field=models.ForeignKey(to='goldenbraid.Feature'),
        ),
        migrations.AddField(
            model_name='featureprop',
            name='type',
            field=models.ForeignKey(to='goldenbraid.Cvterm'),
        ),
        migrations.AddField(
            model_name='feature',
            name='dbxref',
            field=models.ForeignKey(to='goldenbraid.Dbxref'),
        ),
        migrations.AddField(
            model_name='feature',
            name='type',
            field=models.ForeignKey(to='goldenbraid.Cvterm'),
        ),
        migrations.AddField(
            model_name='feature',
            name='vector',
            field=models.ForeignKey(to='goldenbraid.Feature', null=True),
        ),
        migrations.AddField(
            model_name='experimentsubfeature',
            name='experiment',
            field=models.ForeignKey(to='goldenbraid.Experiment'),
        ),
        migrations.AddField(
            model_name='experimentsubfeature',
            name='feature',
            field=models.ForeignKey(to='goldenbraid.Feature'),
        ),
        migrations.AddField(
            model_name='experimentproptext',
            name='experiment',
            field=models.ForeignKey(to='goldenbraid.Experiment'),
        ),
        migrations.AddField(
            model_name='experimentpropnumeric',
            name='experiment',
            field=models.ForeignKey(to='goldenbraid.Experiment'),
        ),
        migrations.AddField(
            model_name='experimentpropnumeric',
            name='type',
            field=models.ForeignKey(to='goldenbraid.Cvterm'),
        ),
        migrations.AddField(
            model_name='experimentpropimage',
            name='experiment',
            field=models.ForeignKey(to='goldenbraid.Experiment'),
        ),
        migrations.AddField(
            model_name='experimentpropexcel',
            name='experiment',
            field=models.ForeignKey(to='goldenbraid.Experiment'),
        ),
        migrations.AddField(
            model_name='experimentfeature',
            name='experiment',
            field=models.ForeignKey(to='goldenbraid.Experiment'),
        ),
        migrations.AddField(
            model_name='experimentfeature',
            name='feature',
            field=models.ForeignKey(to='goldenbraid.Feature'),
        ),
        migrations.AddField(
            model_name='experiment',
            name='dbxref',
            field=models.ForeignKey(to='goldenbraid.Dbxref'),
        ),
        migrations.AddField(
            model_name='experiment',
            name='type',
            field=models.ForeignKey(to='goldenbraid.Cvterm'),
        ),
        migrations.AlterUniqueTogether(
            name='experimentsubfeature',
            unique_together=set([('experiment', 'feature')]),
        ),
        migrations.AlterUniqueTogether(
            name='experimentfeature',
            unique_together=set([('experiment', 'feature')]),
        ),
    ]
