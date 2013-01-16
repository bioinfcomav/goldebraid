'''
Created on 2013 urt 15

@author: peio
'''
from django.core.management.base import BaseCommand, CommandError

from goldenbraid.management.commands.add_cvterms import run_command
from goldenbraid.settings import DB
from goldenbraid.views import add_feature

MANDATORY_FIELDS = ('name', 'type', 'vector', 'genbank_file', 'properties')
FAIL_IF_EXISTS = False


class Command(BaseCommand):
    args = '<features_file>'
    help = 'Adds the given cvterm file to the pseudo_chado database'

    def handle(self, *args, **options):
        'Adds the given featuress to the pseudo_chado database'
        if not args:
            raise CommandError('No database and features file given')
        else:
            features_fpath = args[0]
        try:
            run_command(open(features_fpath), DB, load_features)
        except Exception as error:
            raise CommandError(str(error))


def load_features(database, reader):
    '''it loads the features in the database. The input file is a csv with
    these columns:
        name, type, vector, genbank_file, properties

        properties is a list separated by ';' with (key;value)
                key is the name of the cvterm(i must be in cv goldenbraid)
                value is the value of the property
     '''
    for line in reader:
        try:
            name = line['name']
            type_name = line['type']
            vector = line['vector']
            genbank_fpath = line['genbank_file']
            props = line['properties']
        except KeyError:
            msg = 'Malformed line: ' + str(line)
            raise RuntimeError(msg)
        props = dict([prop_pair.split('=') for prop_pair in props.split(';')])
        add_feature(database, name, type_name, vector, genbank_fpath,
                    props=props)



