# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('goldenbraid', '0003_auto_20150915_0610'),
    ]

    operations = [
        migrations.AlterField(
            model_name='experiment',
            name='description',
            field=models.TextField(max_length=3000),
        ),
        migrations.AlterField(
            model_name='experimentproptext',
            name='value',
            field=models.TextField(max_length=3000),
        ),
    ]
