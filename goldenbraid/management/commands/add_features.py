# Copyright 2013 Diego Orzaez, Univ.Politecnica Valencia, Consejo Superior de
# Investigaciones Cientificas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings as proj_settings
from goldenbraid.management.commands.add_cvterms import run_command
from goldenbraid import settings
from goldenbraid.views.feature import add_feature

MANDATORY_FIELDS = ('name', 'type', 'vector', 'genbank_file', 'properties',
                    'owner', 'is_public')
FAIL_IF_EXISTS = False


class Command(BaseCommand):
    help = 'Adds the given cvterm file to the pseudo_chado database'

    def add_arguments(self, parser):
        parser.add_argument("-f", "--feature_file", nargs="?", type=str)

    def handle(self, *args, **options):
        'Adds the given featuress to the pseudo_chado database'
        if not options:
            raise CommandError('No features file given')
        else:
            features_fpath = options["feature_file"]
        try:
            run_command(open(features_fpath), load_features, MANDATORY_FIELDS)
        except Exception as error:
            raise CommandError(str(error))


def load_features(reader):
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
            owner = line['owner']
            is_public = bool(line['is_public'])

        except KeyError:
            msg = 'Malformed line: ' + str(line)
            raise RuntimeError(msg)
        if props:
            prop_list = [prop_pair.strip().split('=') for prop_pair in props.split(';')]
            props = dict([(prop[0], prop[1].split(',')) for prop in prop_list])
        else:
            props = {}

        # remove the file from media root if exists for duplicity problems
        path_in_media = os.path.join(proj_settings.MEDIA_ROOT,
                                     settings.GENBANK_DIR,
                                     os.path.basename(genbank_fpath))
        if os.path.exists(path_in_media):
            os.remove(path_in_media)
        add_feature(name, type_name, vector, open(genbank_fpath, 'rb'),
                    props=props, owner=owner, is_public=is_public)
