from django.conf import settings

DB = getattr(settings, 'GOLDENBRAID_DB', None)
if not DB:
    raise ValueError('GOLDENBRAID_DB is not defined in the settings')
