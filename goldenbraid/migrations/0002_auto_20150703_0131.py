# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('goldenbraid', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='ExperimentPropGenericFile',
            fields=[
                ('experiment_prop_generic_file_id', models.AutoField(serialize=False, primary_key=True)),
                ('file', models.FileField(upload_to=b'result_files')),
                ('description', models.CharField(max_length=255)),
                ('experiment', models.ForeignKey(to='goldenbraid.Experiment')),
            ],
            options={
                'db_table': 'experimentpropgenericfile',
            },
        ),
        migrations.AlterField(
            model_name='feature',
            name='genbank_file',
            field=models.FileField(upload_to=b'genbank_files'),
        ),
    ]
