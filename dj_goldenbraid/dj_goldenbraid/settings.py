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

# Django settings for dj_goldenbraid project.
import os
from os.path import join, dirname, abspath
PROJECT_DIR = abspath(join(__file__, '..', '..'))
DEBUG = True
TEMPLATE_DEBUG = DEBUG

ADMINS = (
    # ('Your Name', 'your_email@example.com'),
)

MANAGERS = ADMINS

DATABASES = {}
SOUTH_TESTS_MIGRATE = False
REMOTE_DATABASE = {
         'ENGINE': 'django.db.backends.postgresql_psycopg2',  # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
         'NAME': 'goldenbraid',  # Or path to database file if using sqlite3.
         'USER': '',  # Not used with sqlite3.
         'PASSWORD': '',  # Not used with sqlite3.
         'HOST': '',  # Set to empty string for localhost. Not used with sqlite3.
         'PORT': '',  # Set to empty string for default. Not used with sqlite3.
      }
LOCAL_DATABASE = {
        'ENGINE': 'django.db.backends.sqlite3',  # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
         'NAME': 'goldenbraid.db',  # Or path to database file if using sqlite3.
         'USER': '',  # Not used with sqlite3.
         'PASSWORD': '',  # Not used with sqlite3.
        'HOST': '',  # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '',  # Set to empty string for default. Not used with sqlite3.
    }

GENOME_DOMETICATION_DB = {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': 'gb_genome_domest.db',
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    }


DATABASES['default'] = LOCAL_DATABASE
GB_GENOME_DOMESTICATION_DBNAME = 'gb_genome_domestication'
DATABASES[GB_GENOME_DOMESTICATION_DBNAME] = GENOME_DOMETICATION_DB

DATABASE_ROUTERS = ['gb_genome_domestication.db_router.GenomeDomesticationRouter']

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# In a Windows environment this must be set to your system time zone.
TIME_ZONE = 'America/Chicago'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale.
USE_L10N = True

# If you set this to False, Django will not use timezone-aware datetimes.
USE_TZ = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/home/media/media.lawrence.com/media/"
MEDIA_ROOT = os.path.join(PROJECT_DIR, 'media')

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://media.lawrence.com/media/", "http://example.com/media/"
MEDIA_URL = '/media/'

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"
STATIC_ROOT = os.path.join(PROJECT_DIR, 'static')

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
STATIC_URL = '/static/'

# Additional locations of static files
STATICFILES_DIRS = (
    # Put strings here, like "/home/html/static" or "C:/www/django/static".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
#    'django.contrib.staticfiles.finders.DefaultStorageFinder',
)

# Make this unique, and don't share it with anybody.
SECRET_KEY = '9a$4#47d&amp;a1zq2l--9b5-q0jnljkmsyujfk41x-he@p9g(rc1a'

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    # Uncomment the next line for simple clickjacking protection:
    # 'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'dj_goldenbraid.urls'
TEST_RUNNER = 'django.test.runner.DiscoverRunner'
# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'dj_goldenbraid.wsgi.application'

TEMPLATE_DIRS = (
    # Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
    os.path.join(PROJECT_DIR, 'templates'),
)
TEMPLATE_CONTEXT_PROCESSORS = (
                               'django.template.context_processors.request',
                               'django.contrib.auth.context_processors.auth')
INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    # Uncomment the next line to enable the admin:
    'django.contrib.admin',
    # Uncomment the next line to enable admin documentation:
    # 'django.contrib.admindocs',
    'goldenbraid',
    'gb_genome_domestication',
    'restcmd_client',
    'django_tables2',
)

# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error when DEBUG=False.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        }
    },
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'filters': ['require_debug_false'],
            'class': 'django.utils.log.AdminEmailHandler'
        }
    },
    'loggers': {
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': True,
        },
    }
}
GOLDENBRAID_SEARCH_MENU_TYPE_CHOICES = ('PROM+UTR+ATG', 'CDS', 'TER', 'TU',
                                        'Other', 'SP')

# restcmd
ARAB_DB = abspath(join(dirname(__file__), '..', '..',
                       'gb_genome_domestication', 'tests', 'data', 'at.fasta'))
RESTCMD_SERVER_URL = 'http://localhost:8001/cmd_server/'
RESTCMD_SERVER_USER = 'testclient'
RESTCMD_SERVER_PASS = 'test'
BLAST_VIEW_CONFIG = {os.path.basename(ARAB_DB):
                        {'subject_url': '/genome_domestication/feature/%s'}}

RESTCMD_TOOL_CONFIG = {'blastplus': {'database': {'choices':
                                                  [[ARAB_DB, 'arabidopsis'],
                                                   ['tair8_cdna', 'arabid8']]}}
                       }
