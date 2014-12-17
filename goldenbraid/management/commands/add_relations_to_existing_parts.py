'''
Created on 2014 abe 17

@author: peio
'''
import os
from django.core.management.base import BaseCommand, CommandError

from Bio import SeqIO
from goldenbraid.models import Feature
from goldenbraid.views.feature_views import add_relations


class Command(BaseCommand):
    help = 'Adds the given cvterm file to the pseudo_chado database'

    def handle(self, *args, **options):
        'Adds the given featuress to the pseudo_chado database'
        try:
            add_relations_to_existing_parts()
        except Exception as error:
            raise CommandError(str(error))


def add_relations_to_existing_parts():
    for feature in Feature.objects.all():
        if os.path.exists(feature.genbank_file.path):
            seq = SeqIO.read(feature.genbank_file.path, 'gb')
            add_relations(feature, seq)
