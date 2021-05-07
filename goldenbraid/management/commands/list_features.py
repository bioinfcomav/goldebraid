import os
import csv
import re
from django.core.management.base import BaseCommand, CommandError
from goldenbraid.models import Feature, FeatureRelationship
from gbdb import settings
from django.db.models import Q

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide


# GB_UC_5CA
# GB_UA_2131
# GB_UA_2132
# GB_UA_2133

class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('changes_fpath', type=str,
                            help="file with changes to do")

    def handle(self, *args, **kwargs):
        'Searches features using this genbank'
        if not kwargs:
            raise CommandError('No file provided')
        else:
            changes_fhand = open(kwargs['changes_fpath'], "r")
        # try:
        run_command(changes_fhand)
        # except Exception as error:
        #     raise CommandError(str(error))


def get_users(users):
    if users.upper() == "NO" or users.upper() == "ALL":
        return users.upper()
    else:
        users = users.split(",")
        return users


def run_command(fhand):

    report = {}
    features_to_change = []
    features_to_check = []
    total_features = []
    for genbank in csv.DictReader(fhand, delimiter="\t"):
        uniquename = genbank["genbank"]
        original_seq = genbank["original_sequence"]
        modified_seq = genbank["modified_sequence"]
        owner = genbank["owner"]
        users_to_propagate = get_users(genbank["propagate_to"])
        feature = Feature.objects.get(Q(uniquename__iexact=uniquename))
        features_to_change.append(feature)
        features_to_check.append(feature)

        while features_to_check:
            feature = features_to_check.pop(0)
            related_features = FeatureRelationship.objects.filter(Q(subject__uniquename=feature.uniquename))
            for related_feature in related_features:
            	total_features.append(related_feature)
            	feature = related_feature.object
            	feature_owner = str(feature.owner)
            	if users_to_propagate == "NO":
            		if feature_owner == owner:
            			features_to_change.append(feature)
            			features_to_check.append(feature)
            	elif users_to_propagate == "ALL":
            		features_to_change.append(feature)
            		features_to_check.append(feature)
            	elif feature_owner == owner or feature_owner in users_to_propagate:
            		features_to_change.append(feature)
            		features_to_check.append(feature)
            	else:
                	continue
    for feature in features_to_change:
        uniquename = feature.uniquename
    print("feature\tresidues\tgenbank")
    for feature in features_to_change:
        line = feature.uniquename
        print(line)
