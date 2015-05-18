# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('goldenbraid', '0001_initial'),
    ]

    operations = [
        migrations.RenameField(
            model_name='experimentpropexcel',
            old_name='image',
            new_name='excel',
        ),
    ]
