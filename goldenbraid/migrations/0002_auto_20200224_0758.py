# Generated by Django 2.2.9 on 2020-02-24 13:58

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('goldenbraid', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='feature',
            name='field_description',
            field=models.CharField(default='No description', max_length=255),
        ),
        migrations.AddField(
            model_name='feature',
            name='field_owner',
            field=models.CharField(default='guest', max_length=255),
        ),
    ]