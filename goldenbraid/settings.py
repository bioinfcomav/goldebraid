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
                              TER_CRYSPER, PROM_5UTR_NTAG, PROM_5UTR, PROX,
                              DIST_PROX, CORE_5UTR, DIST, INTERACTION_ADAPTOR,
                              PROM_5UTR_MIR173, NTAG, CDS, CDS2_CTAG,
                              CDS1_CDS2, CTAG, FS5, TARGET, FS3, GOI, INT, IOG,
                              FGOI, UTR3_TERM, CDS1)


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
CATEGORIES[PROM_5UTR_NTAG] = ('PROM+UTR+ATG', 'GGAG', 'AATG')
CATEGORIES[PROM_5UTR] = ('PROM+UTR', 'GGAG', 'CCAT')
CATEGORIES[DIST_PROX] = ('OP', 'GGAG', 'TCCC')
CATEGORIES[CORE_5UTR] = ('MinPROM', 'TCCC', 'AATG')
CATEGORIES[DIST] = ('PROM', 'GGAG', 'TGAC')
CATEGORIES[PROX] = ('OP', 'TGAC', 'TCCC')
CATEGORIES[INTERACTION_ADAPTOR] = ('INTERACTION ADAPTOR', 'GGAG', 'AATG')
CATEGORIES[PROM_5UTR_MIR173] = ('PROM+UTR+mir173', 'GGAG', 'CCAT')
CATEGORIES[NTAG] = ('NT', 'CCAT', 'AATG')
CATEGORIES[CDS] = ('CDS', 'AATG', 'GCTT')
CATEGORIES[CDS1] = ('SP', 'AATG', 'AGCC')
CATEGORIES[CDS2_CTAG] = ('CDS', 'AGCC', 'GCTT')
CATEGORIES[CDS1_CDS2] = ('CDS', 'AATG', 'TTCG')
CATEGORIES[CTAG] = ('CT', 'TTCG', 'GCTT')
CATEGORIES[FS5] = ('''5'FS''', 'CCAT', 'GTAG')
CATEGORIES[TARGET] = ('Target', 'GTAG', 'TCTC')
CATEGORIES[FS3] = ('''3'FS''', 'TCTC', 'GCTT')
CATEGORIES[GOI] = ('goi', 'CCAT', 'AGCC')
CATEGORIES[INT] = ('int', 'AGCC', 'TTCG')
CATEGORIES[IOG] = ('iog', 'TTCG', 'GCTT')
CATEGORIES[FGOI] = ('goi', 'CCAT', 'GCTT')
CATEGORIES[UTR3_TERM] = ('TER', 'GCTT', 'CGCT')


CRYSPER_CATEGORIES = OrderedDict()
CRYSPER_CATEGORIES[PROM_DICOT] = (PROM_DICOT, 'GGAG', 'ATTG')
CRYSPER_CATEGORIES[PROM_MONOCOT] = (PROM_MONOCOT, 'GGAG', 'GGCA')
CRYSPER_CATEGORIES[TARGET_DICOT] = (TARGET_DICOT, 'ATTG', 'GTTT')
CRYSPER_CATEGORIES[TARGET_MONOCOT] = (TARGET_MONOCOT, 'GGCA', 'GTTT')
CRYSPER_CATEGORIES[TER_CRYSPER] = (TER_CRYSPER, 'GTTT', 'CGCT')
CRYSPER_TARGETS_TO_DOMESTICATE = (TARGET_MONOCOT, TARGET_DICOT)

MANDATORY_DOMEST_ENZYMES = ('BsmBI', 'BsaI')
OPTIONAL_DOMEST_ENZYMES = ('BtgZI', 'BpiI')


PARTS_TO_ASSEMBLE = {'basic': [CATEGORIES[PROM_5UTR_NTAG],
                               CATEGORIES[CDS] ,
                               CATEGORIES[UTR3_TERM]],
                     'secreted': [CATEGORIES[PROM_5UTR_NTAG],
                                  CATEGORIES[CDS1],
                                  CATEGORIES[CDS1_CDS2],
                                  CATEGORIES[UTR3_TERM]],
                     'ct-fusion': [CATEGORIES[PROM_5UTR_NTAG],
                                   CATEGORIES[CDS1_CDS2],
                                   CATEGORIES[CTAG],
                                   CATEGORIES[UTR3_TERM]],
                     'nt-fusion': [CATEGORIES[PROM_5UTR],
                                   CATEGORIES[NTAG],
                                   CATEGORIES[CDS],
                                   CATEGORIES[UTR3_TERM]],
                     'nt-ct-fusion': [CATEGORIES[PROM_5UTR],
                                      CATEGORIES[NTAG],
                                      CATEGORIES[CDS1_CDS2],
                                      CATEGORIES[CTAG],
                                      CATEGORIES[UTR3_TERM]],
                     'operated-promoter-a': [CATEGORIES[DIST_PROX],
                                             CATEGORIES[CORE_5UTR],
                                             CATEGORIES[CDS],
                                             CATEGORIES[UTR3_TERM]],
                     'operated-promoter-b': [CATEGORIES[DIST],
                                             CATEGORIES[PROX],
                                             CATEGORIES[CORE_5UTR],
                                             CATEGORIES[CDS],
                                             CATEGORIES[UTR3_TERM]],
                     'protein-interaction': [CATEGORIES[INTERACTION_ADAPTOR],
                                             CATEGORIES[CDS],
                                             CATEGORIES[UTR3_TERM]],
                     'amiRNA':  [CATEGORIES[PROM_5UTR],
                                 CATEGORIES[FS5],
                                 CATEGORIES[TARGET],
                                 CATEGORIES[FS3],
                                 CATEGORIES[UTR3_TERM]],
                     'hpRNA':  [CATEGORIES[PROM_5UTR],
                                CATEGORIES[GOI],
                                CATEGORIES[INT],
                                CATEGORIES[IOG],
                                CATEGORIES[UTR3_TERM]],
                     'tasiRNA':  [CATEGORIES[PROM_5UTR_MIR173],
                                  CATEGORIES[FGOI],
                                  CATEGORIES[UTR3_TERM]],
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
