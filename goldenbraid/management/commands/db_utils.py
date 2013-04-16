'''
Created on 2013 urt 15

@author: peio
'''
from goldenbraid.models import Cv, Cvterm


def get_cv(name):
    'It gets the db object if exists else none'
    try:
        cv = Cv.objects.get(name=name)
    except Cv.DoesNotExist:
        cv = None
    return cv


def get_or_load_cv(name, definition=None, fail_if_exists=False):
    'It loads the db into the database'
    cv = get_cv(name)
    if cv:
        if fail_if_exists:
            raise RuntimeError('cv %s already in the database' % name)
    else:
        if definition is None:
            definition = name
        cv = Cv.objects.create(name=name, definition=definition)

    return cv


def get_cvterm(cv=None, name=None, dbxref=None):
    'It gets the cvterm object else None'
    if cv is None and name is None and dbxref is None:
        raise RuntimeError('it needs at least one argument')
    try:
        cvterm = Cvterm.objects.get(cv=cv, name=name)
    except Cvterm.DoesNotExist:
        cvterm = None
    return cvterm


def get_or_load_cvterm(cv, name, definition, fail_if_exists=False):
    '''It loads into the database the cvterm with the given data.
    all arguments are string.
    '''
    cvterm = get_cvterm(cv, name)
    if cvterm:
        if fail_if_exists:
            raise RuntimeError('cvterm already in database')
    else:
        cvterm = Cvterm.objects.create(cv=cv, name=name,
                                                       definition=definition)
    return cvterm

