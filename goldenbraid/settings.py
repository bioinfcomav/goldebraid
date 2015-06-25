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
from collections import OrderedDict

from django.conf import settings

import goldenbraid
from goldenbraid.tags import (MODULE_TYPE_NAME, TU_TYPE_NAME, TARGET_DICOT,
                              TARGET_MONOCOT, PROM_DICOT, PROM_MONOCOT,
                              TER_CRYSPER)


GENBANK_DIR = getattr(settings, 'GOLDENBRAID_GENBANK_DIR', 'genbank_files')
RESULTS_DIR = getattr(settings, 'GOLDENBRAID_RESULTS_DIR', 'result_files')

REBASE = os.path.join(goldenbraid.__path__[0], 'rebase', 'withrefm.301')
REBASE_FILE = getattr(settings, 'GOLDENBRAID_REBASE_FILE', REBASE)
SEARCH_MENU_TYPE_CHOICES = getattr(settings,
                                   'GOLDENBRAID_SEARCH_MENU_TYPE_CHOICES',
                                   None)
DOMESTICATION_DEFAULT_MELTING_TEMP = 50
DOMESTICATION_MIN_OLIGO_LENGTH = 20
MINIMUN_PCR_LENGTH = 50

BIPARTITE_ALLOWED_PARTS = (TU_TYPE_NAME, MODULE_TYPE_NAME)

CATEGORIES = OrderedDict()
CATEGORIES['PROM+5UTR+NTAG (A1-A2-A3-B1-B2)'] = ('PROM+UTR+ATG', 'GGAG',
                                                 'AATG')
CATEGORIES['PROM+5UTR (A1-A2-A3-B1)'] = ('PROM+UTR', 'GGAG', 'CCAT')
CATEGORIES['OP (A1-A2)'] = ('OP', 'GGAG', 'TCCC')
CATEGORIES['MinPROM (A3-B1-B2)'] = ('MinPROM', 'TCCC', 'AATG')
CATEGORIES['PROM (A1)'] = ('PROM', 'GGAG', 'TGAC')
CATEGORIES['OP (A2)'] = ('OP', 'TGAC', 'TCCC')
CATEGORIES['INTERACTION ADAPTOR (A1-A2-A3-B1-B2b)'] = ('INTERACTION ADAPTOR',
                                                       'GGAG', 'AATG')
CATEGORIES['PROM+5UTR+mir173 (A1-A2-A3-B1b)'] = ('PROM+UTR+mir173', 'GGAG',
                                                 'CCAT')
CATEGORIES['NTAG (B2)'] = ('NT', 'CCAT', 'AATG')
CATEGORIES['CDS (B3-B4-B5)'] = ('CDS', 'AATG', 'GCTT')
CATEGORIES['SP (B3)'] = ('SP', 'AATG', 'AGCC')
CATEGORIES['CDS (B4-B5)'] = ('CDS', 'AGCC', 'GCTT')
CATEGORIES['CDS (B3-B4)'] = ('CDS', 'AATG', 'TTCG')
CATEGORIES['CTAG (B5)'] = ('CT', 'TTCG', 'GCTT')
CATEGORIES["5'FS (B2-B3b)"] = ('''5'FS''', 'CCAT', 'GTAG')
CATEGORIES['Target (B4b)'] = ('Target', 'GTAG', 'TCTC')
CATEGORIES["3'FS (B5b)"] = ('''3'FS''', 'TCTC', 'GCTT')
CATEGORIES['goi (B2-B3)'] = ('goi', 'CCAT', 'AGCC')
CATEGORIES['int (B4)'] = ('int', 'AGCC', 'TTCG')
CATEGORIES['iog (B5)'] = ('iog', 'TTCG', 'GCTT')
CATEGORIES['goi (B2-B3-B4-B5)'] = ('goi', 'CCAT', 'GCTT')
CATEGORIES['3UTR+TERM (B6-C1)'] = ('TER', 'GCTT', 'CGCT')


CRYSPER_CATEGORIES = OrderedDict()
CRYSPER_CATEGORIES[PROM_DICOT] = (PROM_DICOT, 'GGAG', 'ATTG')
CRYSPER_CATEGORIES[PROM_MONOCOT] = (PROM_MONOCOT, 'GGAG', 'GGCA')
CRYSPER_CATEGORIES[TARGET_DICOT] = (TARGET_DICOT, 'ATTG', 'GTTT')
CRYSPER_CATEGORIES[TARGET_MONOCOT] = (TARGET_MONOCOT, 'GGCA', 'GTTT')
CRYSPER_CATEGORIES[TER_CRYSPER] = (TER_CRYSPER, 'GTTT', 'CGCT')
CRYSPER_TARGETS_TO_DOMESTICATE = (TARGET_MONOCOT, TARGET_DICOT)

MANDATORY_DOMEST_ENZYMES = ('BsmBI', 'BsaI')
OPTIONAL_DOMEST_ENZYMES = ('BtgZI', 'BpiI')


PARTS_TO_ASSEMBLE = {'basic': [CATEGORIES['PROM+5UTR+NTAG (A1-A2-A3-B1-B2)'],
                               CATEGORIES['CDS (B3-B4-B5)'] ,
                               CATEGORIES['3UTR+TERM (B6-C1)']],
                     'secreted': [CATEGORIES['PROM+5UTR+NTAG (A1-A2-A3-B1-B2)'],
                                  CATEGORIES['SP (B3)'],
                                  CATEGORIES['CDS (B4-B5)'],
                                  CATEGORIES['3UTR+TERM (B6-C1)']],
                     'ct-fusion': [CATEGORIES['PROM+5UTR+NTAG (A1-A2-A3-B1-B2)'],
                                   CATEGORIES['CDS (B3-B4)'],
                                   CATEGORIES['CTAG (B5)'],
                                   CATEGORIES['3UTR+TERM (B6-C1)']],
                     'nt-fusion': [CATEGORIES['PROM+5UTR (A1-A2-A3-B1)'],
                                   CATEGORIES['NTAG (B2)'],
                                   CATEGORIES['CDS (B3-B4-B5)'],
                                   CATEGORIES['3UTR+TERM (B6-C1)']],
                     'nt-ct-fusion': [CATEGORIES['PROM+5UTR (A1-A2-A3-B1)'],
                                      CATEGORIES['NTAG (B2)'],
                                      CATEGORIES['CDS (B3-B4)'],
                                      CATEGORIES['CTAG (B5)'],
                                      CATEGORIES['3UTR+TERM (B6-C1)']],
                     'operated-promoter-a': [CATEGORIES['OP (A1-A2)'],
                                             CATEGORIES['MinPROM (A3-B1-B2)'],
                                             CATEGORIES['CDS (B3-B4-B5)'],
                                             CATEGORIES['3UTR+TERM (B6-C1)']],
                     'operated-promoter-b': [CATEGORIES['PROM (A1)'],
                                             CATEGORIES['OP (A2)'],
                                             CATEGORIES['MinPROM (A3-B1-B2)'],
                                             CATEGORIES['CDS (B3-B4-B5)'],
                                             CATEGORIES['3UTR+TERM (B6-C1)']],
                     'protein-interaction': [CATEGORIES['INTERACTION ADAPTOR (A1-A2-A3-B1-B2b)'],
                                             CATEGORIES['CDS (B3-B4-B5)'] ,
                                             CATEGORIES['3UTR+TERM (B6-C1)']],
                     'amiRNA':  [CATEGORIES['PROM+5UTR (A1-A2-A3-B1)'],
                                 CATEGORIES["5'FS (B2-B3b)"],
                                 CATEGORIES['Target (B4b)'],
                                 CATEGORIES["3'FS (B5b)"],
                                 CATEGORIES['3UTR+TERM (B6-C1)']],
                     'hpRNA':  [CATEGORIES['PROM+5UTR (A1-A2-A3-B1)'],
                                CATEGORIES['goi (B2-B3)'],
                                CATEGORIES['int (B4)'],
                                CATEGORIES['iog (B5)'],
                                CATEGORIES['3UTR+TERM (B6-C1)']],
                     'tasiRNA':  [CATEGORIES['PROM+5UTR+mir173 (A1-A2-A3-B1b)'],
                                  CATEGORIES['goi (B2-B3-B4-B5)'],
                                  CATEGORIES['3UTR+TERM (B6-C1)']],

                     'gRNA_monocot': [CRYSPER_CATEGORIES[PROM_MONOCOT],
                                      CRYSPER_CATEGORIES[TARGET_MONOCOT],
                                      CRYSPER_CATEGORIES[TER_CRYSPER]],
                     'gRNA_dicot': [CRYSPER_CATEGORIES[PROM_DICOT],
                                    CRYSPER_CATEGORIES[TARGET_DICOT],
                                    CRYSPER_CATEGORIES[TER_CRYSPER]],
                     }
UT_PREFIX = PARTS_TO_ASSEMBLE['basic'][0][1]
UT_SUFFIX = PARTS_TO_ASSEMBLE['basic'][-1][2]

SITE_A = UT_PREFIX
SITE_B = UT_SUFFIX
SITE_C = 'GTCA'

DOMESTICATION_VECTORS_IN_GB = ('pUPD',)
DOMESTICATED_VECTOR = 'pUPD'

PUPD_PREFIX = 'CTCG'

# added to create all oligos in domesticator
OLIGO_UNIVERSAL = 'GCGCCGTCTCG'
ASSEMBLED_SEQ = 'GB_UA'
DOMESTICATED_SEQ = 'GB_UD'
CRYSPER_SEQ = 'GB_UC'
EXPERIMENT_ID_PREFIX = 'GB_EXP'
