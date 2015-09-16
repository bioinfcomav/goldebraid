# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('goldenbraid', '0002_auto_20150703_0131'),
    ]

    operations = [
        migrations.AlterField(
            model_name='experimentpropnumeric',
            name='value',
            field=models.FloatField(null=True),
        ),
    ]
