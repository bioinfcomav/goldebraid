import os

from django.conf import settings
import goldenbraid

DB = getattr(settings, 'GOLDENBRAID_DB', None)
if not DB:
    raise ValueError('GOLDENBRAID_DB is not defined in the settings')


GENBANK_DIR = getattr(settings, 'GOLDENBRAID_GENBANK_DIR', 'genbank_files')
REBASE = os.path.join(goldenbraid.__path__[0], 'rebase', 'withrefm.301')
REBASE_FILE = getattr(settings, 'GOLDENBRAID_REBASE_FILE', REBASE)
