'''
Created on 2014 uzt 1

@author: peio
'''
from gb_genome_domestication.models import Db
from gb_genome_domestication.settings import DB


def get_db(name):
    'It gets the db object if exists else none'
    try:
        db = Db.objects.using(DB).get(name=name)
    except Db.DoesNotExist:
        db = None
    return db


def get_or_load_db(name, urlprefix, fail_if_exists=False):
    'It loads the db into the database'
    db = get_db(name)
    if db:
        if fail_if_exists:
            raise RuntimeError('db %s already in the database' % name)
    else:
        db = Db.objects.create(name=name, urlprefix=urlprefix)

    return db

