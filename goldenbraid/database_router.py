'''
Created on 2013 urt 16

@author: peio
'''
'''
Created on 27/01/2011

@author: jose
'''

from goldenbraid import settings

APP_LABEL = 'goldenbraid'
DB = settings.DB


class GoldenbraidRouter(object):
    """A router to control all database operations on models in
    the FeatureViewer application"""

    def db_for_read(self, model, **hints):
        "Point all operations on feature_viewer models to 'other'"
        if model._meta.app_label == APP_LABEL:
            return DB
        return None

    def db_for_write(self, model, **hints):
        "Point all operations on feature_viewer models to 'other'"
        if model._meta.app_label == APP_LABEL:
            return DB
        return None

    def allow_relation(self, obj1, obj2, **hints):
        "Allow any relation if a model in feature_viewer is involved"
        if (obj1._meta.app_label == APP_LABEL or
            obj2._meta.app_label == APP_LABEL):
            return True
        return None

    def allow_syncdb(self, db, model):
        "Make sure the feature_viewer app only appears on the 'other' db"
        if db == DB:
            return model._meta.app_label == APP_LABEL
        elif model._meta.app_label == APP_LABEL:
            return False
        return None
