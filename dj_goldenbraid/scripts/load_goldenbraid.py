#!/usr/bin/env python

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
from os.path import join
from subprocess import check_call

DJANGO_DIR = '/home/peio/devel/goldenbraid/dj_goldenbraid/'
GOLDENBRAID_DB_PATH = os.path.join(DJANGO_DIR, 'goldenbraid.db')

DB_LOADING = '/home/peio/devel/goldenbraid/db_loading'
CVTERMS_PATH = join(DB_LOADING, 'cvterms.csv')
VECTORS_PATH = join(DB_LOADING, 'vectors.csv')
PARTS_PATH = join(DB_LOADING, 'parts.csv')
TUS_PATH = join(DB_LOADING, 'tus.csv')


def create_database():
    'It creates the melonomics chado database'
    check_call(['rm', GOLDENBRAID_DB_PATH])
    _syncdb_django()


def _syncdb_django():
    'it runs the syncsdb of the goldenbraid database'
    old_cwd = os.getcwd()
    os.chdir(DJANGO_DIR)
    cmd = ['python', 'manage.py', 'syncdb']
    check_call(cmd)
    os.chdir(old_cwd)


def _add_feature(fpath):
    'It adds feature cvterm relations on the database'
    old_cwd = os.getcwd()
    os.chdir(DJANGO_DIR)
    cmd = ['python', 'manage.py', 'add_features', fpath]
    check_call(cmd)
    os.chdir(old_cwd)


def _add_cvterms(fpath):
    'It adds feature cvterm relations on the database'
    old_cwd = os.getcwd()
    os.chdir(DJANGO_DIR)
    cmd = ['python', 'manage.py', 'add_cvterms', fpath]
    check_call(cmd)
    os.chdir(old_cwd)


def load_initial_data():
    insert = "insert into db(name, description, urlprefix, url) values ('{0}',"
    insert += " '{1}', '{2}', '{3}');"
    insert = insert.format('goldenbraid', '', '/', '/')
    cmd = ['sqlite3', GOLDENBRAID_DB_PATH, insert]
    check_call(cmd)
    # load cvterms
    _add_cvterms(CVTERMS_PATH)
    # load vectors
    _add_feature(VECTORS_PATH)
    # load parts
    _add_feature(PARTS_PATH)
    # load uts
    _add_feature(TUS_PATH)


def main():

    actions = [create_database, load_initial_data]
    # actions = all_actions
    for action in actions:
          action()

if __name__ in '__main__':
    main()

