import os
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict


from django.conf import settings

import goldenbraid
from goldenbraid.tags import MODULE_TYPE_NAME, TU_TYPE_NAME


GENBANK_DIR = getattr(settings, 'GOLDENBRAID_GENBANK_DIR', 'genbank_files')
REBASE = os.path.join(goldenbraid.__path__[0], 'rebase', 'withrefm.301')
REBASE_FILE = getattr(settings, 'GOLDENBRAID_REBASE_FILE', REBASE)
SEARCH_MENU_TYPE_CHOICES = getattr(settings,
                                   'GOLDENBRAID_SEARCH_MENU_TYPE_CHOICES',
                                   None)
DOMESTICATION_DEFAULT_MELTING_TEMP = 50
DOMESTICATION_MIN_OLIGO_LENGTH = 20

PARTS_TO_ASSEMBLE = {'basic': [('PROM+UTR+ATG', 'GGAG', 'AATG'),
                               ('CDS', 'AATG', 'GCTT'),
                               ('TER', 'GCTT', 'CGCT')],
                     'secreted': [('PROM+UTR+ATG', 'GGAG', 'AATG'),
                                 ('SP', 'AATG', 'AGCC'),
                                 ('CDS', 'AGCC', 'GCTT'),
                                 ('TER', 'GCTT', 'CGCT')],
                     'ct-fusion': [('PROM+UTR+ATG', 'GGAG', 'AATG'),
                                   ('CDS', 'AATG', 'GCAG'),
                                   ('CT', 'GCAG', 'GCTT'),
                                   ('TER', 'GCTT', 'CGCT')],
                     'nt-fusion': [('PROM+UTR', 'GGAG', 'CCAT'),
                                   ('NT', 'CCAT', 'AATG'),
                                   ('CDS', 'AATG', 'GCTT'),
                                   ('TER', 'GCTT', 'CGCT')],
                     'nt-ct-fusion': [('PROM+UTR', 'GGAG', 'CCAT'),
                                      ('NT', 'CCAT', 'AATG'),
                                      ('CDS', 'AATG', 'GCAG'),
                                      ('CT', 'GCAG', 'GCTT'),
                                      ('TER', 'GCTT', 'CGCT')],
                     'operated-promoter-a': [('OP', 'GGAG', 'TCCC'),
                                           ('MinPROM', 'TCCC', 'AATG'),
                                           ('CDS', 'AATG', 'GCTT'),
                                           ('TER', 'GCTT', 'CGCT')],
                     'operated-promoter-b': [('PROM', 'GGAG', 'TGAC'),
                                           ('OP', 'TGAC', 'TCCC'),
                                           ('MinPROM', 'TCCC', 'AATG'),
                                           ('CDS', 'AATG', 'GCTT'),
                                           ('TER', 'GCTT', 'CGCT')],
                     'protein-interaction': [('INTERACTION ADAPTOR', 'GGAG',
                                              'AATG'),
                                             ('CDS', 'AATG', 'GCTT'),
                                             ('TER', 'GCTT', 'CGCT')],
                     'amiRNA':  [('PROM+UTR', 'GGAG', 'CCAT'),
                                 ('''5'FS''', 'CCAT', 'GTAG'),
                                 ('Target', 'GTAG', 'TCTC'),
                                 ('''3'FS''', 'TCTC', 'GCTT'),
                                 ('TER', 'GCTT', 'CGCT')],
                     'hpRNA':  [('PROM+UTR', 'GGAG', 'CCAT'),
                                ('goi', 'CCAT', 'AGCC'),
                                ('int', 'AGCC', 'GCAG'),
                                ('iog', 'GCAG', 'GCTT'),
                                ('TER', 'GCTT', 'CGCT')],
                     'tasiRNA':  [('PROM+UTR+mir173', 'GGAG', 'CCAT'),
                                  ('goi', 'CCAT', 'GCTT'),
                                  ('TER', 'GCTT', 'CGCT')]
                     }
UT_PREFIX = PARTS_TO_ASSEMBLE['basic'][0][1]
UT_SUFFIX = PARTS_TO_ASSEMBLE['basic'][-1][2]

SITE_A = UT_PREFIX
SITE_B = UT_SUFFIX
SITE_C = 'GTCA'

BIPARTITE_ALLOWED_PARTS = (TU_TYPE_NAME, MODULE_TYPE_NAME)

CATEGORIES = OrderedDict()
CATEGORIES['01-02-03-11-12 (PROM+UTR+ATG)'] = ('PROM+UTR+ATG', 'GGAG', 'AATG')
CATEGORIES['01-02-03-11 (PROM+UTR)'] = ('PROM+UTR', 'GGAG', 'CCAT')
CATEGORIES['01-02 (OP)'] = ('OP', 'GGAG', 'TCCC')
CATEGORIES['03-11-12 (MinPROM)'] = ('MinPROM', 'TCCC', 'AATG')
CATEGORIES['01 (PROM)'] = ('PROM', 'GGAG', 'TGAC')
CATEGORIES['02 (OP)'] = ('OP', 'TGAC', 'TCCC')
CATEGORIES['01-02-03-11-12B (INTERACTION ADAPTOR)'] = \
                                        ('INTERACTION ADAPTOR', 'GGAG', 'AATG')
CATEGORIES['01-02-03-11-C (PROM+UTR+mir173)'] = \
                                        ('PROM+UTR+mir173', 'GGAG', 'CCAT')
CATEGORIES['12 (NT)'] = ('NT', 'CCAT', 'AATG')
CATEGORIES['13-14-15-16 (CDS)'] = ('CDS', 'AATG', 'GCTT')
CATEGORIES['13 (SP)'] = ('SP', 'AATG', 'AGCC')
CATEGORIES['14-15-16 (CDS)'] = ('CDS', 'AGCC', 'GCTT')
CATEGORIES['13-14-15 (CDS)'] = ('CDS', 'AATG', 'GCAG')
CATEGORIES['16 (CT)'] = ('CT', 'GCAG', 'GCTT')
CATEGORIES["12-13B (5'FS)"] = ('''5'FS''', 'CCAT', 'GTAG')
CATEGORIES['14B-15B (Target)'] = ('Target', 'GTAG', 'TCTC')
CATEGORIES["16B (3'FS)"] = ('''3'FS''', 'TCTC', 'GCTT')
CATEGORIES['12-13 (GOI)'] = ('goi', 'CCAT', 'AGCC')
CATEGORIES['14-15(INT)'] = ('int', 'AGCC', 'GCAG')
CATEGORIES['16 (IOG)'] = ('iog', 'GCAG', 'GCTT')
CATEGORIES['12-13-14-15-16 (GOI)'] = ('goi', 'CCAT', 'GCTT')
CATEGORIES['17-21 (TER)'] = ('TER', 'GCTT', 'CGCT')

ENZYMES_USED_IN_GOLDENBRAID = ('BsmBI', 'BsaI', 'BtgZI')

PUPD_PREFIX = 'CTCG'

# added to create all oligos in domesticator
OLIGO_UNIVERSAL = 'GCGCCGTCTCG'
ASSEMBLED_SEQ = 'GB_ASSEMB'
DOMESTICATED_SEQ = 'GB_DOMEST'
