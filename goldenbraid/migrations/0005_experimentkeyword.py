# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('goldenbraid', '0004_auto_20150917_0904'),
    ]

    operations = [
        migrations.CreateModel(
            name='ExperimentKeyword',
            fields=[
                ('experiment_keyword_id', models.AutoField(serialize=False, primary_key=True)),
                ('keyword', models.CharField(max_length=50)),
                ('experiment', models.ForeignKey(to='goldenbraid.Experiment')),
            ],
            options={
                'db_table': 'experimentkeyword',
            },
        ),
    ]
