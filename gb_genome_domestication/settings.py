'''
Created on 2014 uzt 1

@author: peio
'''
from django.conf import settings

APP_NAME = 'gb_genome_domestication'

ENZYMES_USED_IN_GOLDENBRAID = ('BsmBI', 'BsaI', 'BtgZI')

DB = getattr(settings, 'GB_GENOME_DOMESTICATION_DBNAME')

DB_URLPREFIX = '/genome_domestication/feature/'

RESTCMD_CLIENT_APP_NAME = 'restcmd_client'
