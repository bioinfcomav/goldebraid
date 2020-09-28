'''
Created on 2014 uzt 1

@author: peio
'''
from gb_genome_domestication import settings

DB_GENOME_DOMESTICATION = settings.DB
APP_NAME = settings.APP_NAME


class GenomeDomesticationRouter(object):
    """
    A router to control all database operations on models in the
    auth application.
    """
    def db_for_read(self, model, **hints):
        """
        Attempts to read auth models go to auth_db.
        """
        if model._meta.app_label == APP_NAME:
            return DB_GENOME_DOMESTICATION
        return None

    def db_for_write(self, model, **hints):
        """
        Attempts to write auth models go to auth_db.
        """
        if model._meta.app_label == APP_NAME:
            return DB_GENOME_DOMESTICATION
        return None

    def allow_relation(self, obj1, obj2, **hints):
        """
        Allow relations if a model in the auth app is involved.
        """
        if obj1._meta.app_label == APP_NAME or \
           obj2._meta.app_label == APP_NAME:
            return True
        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        """
        Make sure the auth app only appears in the 'auth_db'
        database.
        """
        if db == DB_GENOME_DOMESTICATION:
            return app_label == APP_NAME
        elif app_label == APP_NAME:
            return False
        return None
