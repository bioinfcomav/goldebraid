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
                              FGOI, UTR3_TERM, CDS1, FORWARD, REVERSE,
                              OTHER_TYPE_NAME, CRISPR_MULTIPLEXING_TARGET,
                              EDIT_E1, EDIT_E2, EDIT_E3, EDIT_E4, EDIT_E3_N_MINUS_ONE,
                              EDIT_N_MINUS_1, EDIT_NTERM, EDIT_E1_N_MINUS_ONE,
                              EDIT_E2_N, EDIT_E1B, EDIT_E1B_N_MINUS_ONE,
                              REG_R1, REG_RN_MINUS_ONE, REG_NTERM, REG_ALL,
                              CRISPR_EDITING, CRISPR_REGULATION,
                              TRNA, SCAFFOLD, EDIT_E1_N_MINUS_ONE,
                              EDIT_E4_N_MINUS_ONE, EDIT_E1_NTERM,
                              TARGET_CAS12A, PROM_CAS12,
                              CLEAVAGE_SIGNAL,
                              CAS12_LEVEL_MINUS_ONE_2X,
                              CAS12_LEVEL_MINUS_ONE_3X,
                              CAS12_LEVEL_MINUS_ONE_4X,
                              CAS12_LEVEL_MINUS_ONE_5X,
                              CAS12_LEVEL_MINUS_ONE_6X,
                              FUNGAL_PROM_5UTR,
                              FUNGAL_CDS,
                              FUNGAL_UTR3_TERM,
                              REVERSE_MARKER,
                              FORWARD_MARKER,
                              FLANK_5UTR,
                              FLANK_3UTR,
                              FUNGAL_TU_TYPE_NAME,
                              FUNGAL_KNOCK_OUT,
                              )

GENBANK_DIR = getattr(settings, 'GOLDENBRAID_GENBANK_DIR', 'genbank')
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
FUNGAL_BIPARTITE_ALLOWED_PARTS = (FUNGAL_TU_TYPE_NAME, FUNGAL_KNOCK_OUT)

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
CATEGORIES[PROM_CAS12] = (PROM_CAS12, 'GGAG', 'AGAT')
CATEGORIES[CLEAVAGE_SIGNAL] = (CLEAVAGE_SIGNAL, 'GGCC', 'CGCT')


FUNGAL_CATEGORIES = OrderedDict()
FUNGAL_CATEGORIES[FUNGAL_PROM_5UTR] = (FUNGAL_PROM_5UTR, 'GGAG', 'AATG')
FUNGAL_CATEGORIES[FUNGAL_CDS] = (FUNGAL_CDS, 'AATG', 'GCTT')
FUNGAL_CATEGORIES[FUNGAL_UTR3_TERM] = (FUNGAL_UTR3_TERM, 'GCTT', 'CGCT')
FUNGAL_CATEGORIES[REVERSE_MARKER] = (REVERSE_MARKER, "GGAG", "TACT")
FUNGAL_CATEGORIES[FLANK_5UTR] = (FLANK_5UTR, "TACT", "AATG")
FUNGAL_CATEGORIES[FORWARD_MARKER] = (FORWARD_MARKER, "AATG", "GCTT")
FUNGAL_CATEGORIES[FLANK_3UTR] = (FLANK_3UTR, "GCTT", "CGCT")


CRYSPER_CATEGORIES = OrderedDict()
CRYSPER_CATEGORIES[PROM_DICOT] = (PROM_DICOT, 'GGAG', 'ATTG')
CRYSPER_CATEGORIES[PROM_CAS12] = (PROM_CAS12, 'GGAG', 'AGAT')
CRYSPER_CATEGORIES[CLEAVAGE_SIGNAL] = (CLEAVAGE_SIGNAL, 'GGCC', 'CGCT')
CRYSPER_CATEGORIES[PROM_MONOCOT] = (PROM_MONOCOT, 'GGAG', 'GGCA')
CRYSPER_CATEGORIES[TARGET_DICOT] = (TARGET_DICOT, 'ATTG', 'GTTT')
CRYSPER_CATEGORIES[TARGET_MONOCOT] = (TARGET_MONOCOT, 'GGCA', 'GTTT')
CRYSPER_CATEGORIES[TARGET_CAS12A] = (TARGET_CAS12A, 'AGAT', 'GGCC')
CRYSPER_CATEGORIES[TER_CRYSPER] = (TER_CRYSPER, 'GTTT', 'CGCT')
CRYSPER_CATEGORIES[CRISPR_MULTIPLEXING_TARGET] = (CRISPR_MULTIPLEXING_TARGET, 'GTGC', 'GTTT')
CRYSPER_TARGETS_TO_DOMESTICATE = (TARGET_MONOCOT, TARGET_DICOT,
                                  CRISPR_MULTIPLEXING_TARGET,
                                  TARGET_CAS12A)
CRISPR_CAS9_SINGLE_TO_DOMESTICATE = (TARGET_MONOCOT, TARGET_DICOT)

CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE = OrderedDict()
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E1] = (EDIT_E1, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E2] = (EDIT_E2, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E3] = (EDIT_E3, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E4] = (EDIT_E4, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_N_MINUS_1] = (EDIT_N_MINUS_1, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E3_N_MINUS_ONE] = (EDIT_E3_N_MINUS_ONE, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E4_N_MINUS_ONE] = (EDIT_E4_N_MINUS_ONE, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_NTERM] = (EDIT_NTERM, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E1_N_MINUS_ONE] = (EDIT_E1_N_MINUS_ONE, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E1_NTERM] = (EDIT_E1_NTERM, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E2_N] = (EDIT_E2_N, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E1B] = (EDIT_E1B, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E1B_N_MINUS_ONE] = (EDIT_E1B_N_MINUS_ONE, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[REG_R1] = (REG_R1, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[REG_RN_MINUS_ONE] = (REG_RN_MINUS_ONE, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[REG_NTERM] = (REG_NTERM, 'CTCG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[REG_ALL] = (REG_ALL, 'CTCG', 'GTGC')


CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE = OrderedDict()
CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[CAS12_LEVEL_MINUS_ONE_2X] = (CAS12_LEVEL_MINUS_ONE_2X, 'CTCG', 'TGAG')
CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[CAS12_LEVEL_MINUS_ONE_3X] = (CAS12_LEVEL_MINUS_ONE_3X, 'CTCG', 'TGAG')
CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[CAS12_LEVEL_MINUS_ONE_4X] = (CAS12_LEVEL_MINUS_ONE_4X, 'CTCG', 'TGAG')
CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[CAS12_LEVEL_MINUS_ONE_5X] = (CAS12_LEVEL_MINUS_ONE_5X, 'CTCG', 'TGAG')
CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[CAS12_LEVEL_MINUS_ONE_6X] = (CAS12_LEVEL_MINUS_ONE_6X, 'CTCG', 'TGAG')

CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_ZERO = OrderedDict()
CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_ZERO[CAS12_LEVEL_MINUS_ONE_2X] = (CAS12_LEVEL_MINUS_ONE_2X, 'AGAT', 'CGCT')
CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_ZERO[CAS12_LEVEL_MINUS_ONE_3X] = (CAS12_LEVEL_MINUS_ONE_3X, 'AGAT', 'CGCT')
CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_ZERO[CAS12_LEVEL_MINUS_ONE_4X] = (CAS12_LEVEL_MINUS_ONE_4X, 'AGAT', 'CGCT')
CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_ZERO[CAS12_LEVEL_MINUS_ONE_5X] = (CAS12_LEVEL_MINUS_ONE_5X, 'AGAT', 'CGCT')
CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_ZERO[CAS12_LEVEL_MINUS_ONE_6X] = (CAS12_LEVEL_MINUS_ONE_6X, 'AGAT', 'CGCT')

CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO = OrderedDict()
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1] = (EDIT_E1, 'ATTG', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2] = (EDIT_E2, 'GTGC', 'TCGG')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E3] = (EDIT_E3, 'TCGG', 'AGTC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E4] = (EDIT_E4, 'AGTC', 'CGAG')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E3_N_MINUS_ONE] = (EDIT_E3_N_MINUS_ONE, 'TCGG', 'GCAA')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E4_N_MINUS_ONE] = (EDIT_E4_N_MINUS_ONE, 'AGTC', 'GCAA')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_N_MINUS_1] = (EDIT_N_MINUS_1, 'CGAG', 'GCAA')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM] = (EDIT_NTERM, 'GCAA', 'CGCT')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1_N_MINUS_ONE] = (EDIT_E1_N_MINUS_ONE, 'ATTG', 'GCAA')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2_N] = (EDIT_E2_N, 'GTGC', 'GCAA')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1B] = (EDIT_E1B, 'GGCA', 'GTGC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1B_N_MINUS_ONE] = (EDIT_E1B_N_MINUS_ONE, 'GGCA', 'GCAA')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_R1] = (REG_R1, 'ATTG', 'TCCC')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_RN_MINUS_ONE] = (REG_RN_MINUS_ONE, 'TCCC', 'CCAA')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_NTERM] = (REG_NTERM, 'CCAA', 'CGCT')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1_NTERM] = (EDIT_E1_NTERM, 'ATTG', 'CGCT')
CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_ALL] = (REG_ALL, 'ATTG', 'CGCT')



CRYSPR_MULTIPLEX_MODES = [CRISPR_EDITING, CRISPR_REGULATION]

MONOCOT_EDIT_POS = [EDIT_E1B, EDIT_E2, EDIT_E3, EDIT_E4, EDIT_N_MINUS_1,
                    EDIT_NTERM, EDIT_E1B_N_MINUS_ONE, EDIT_E4_N_MINUS_ONE,
                    EDIT_E3_N_MINUS_ONE]
DICOT_EDIT_POS = [EDIT_E1, EDIT_E2, EDIT_E3, EDIT_E4, EDIT_N_MINUS_1,
                  EDIT_NTERM, EDIT_E1_N_MINUS_ONE, EDIT_E4_N_MINUS_ONE,
                  EDIT_E1_NTERM, EDIT_E3_N_MINUS_ONE]


CRYSPR_MULTIPLEX_EDITING_LEVEL_MINUS_ONE = (CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E1],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E2],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E3],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E4],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E3_N_MINUS_ONE],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E4_N_MINUS_ONE],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_N_MINUS_1],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E1_N_MINUS_ONE],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_NTERM],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E1_N_MINUS_ONE],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E2_N],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E1B],
                                            CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[EDIT_E1B_N_MINUS_ONE]
                                            )


CRYSPR_MULTIPLEX_REGULATION_LEVEL_MINUS_ONE = (CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[REG_R1],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[REG_RN_MINUS_ONE],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[REG_NTERM],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[REG_ALL])


CRYSPR_MULTIPLEX_EDITING_LEVEL_ZERO = (CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1],
                                       CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2],
                                       CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E3],
                                       CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E4],
                                       CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E3_N_MINUS_ONE],
                                       CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_N_MINUS_1],
                                       CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM],
                                       CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1_N_MINUS_ONE],
                                       CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2_N],
                                       CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1B],
                                       CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1B_N_MINUS_ONE]
                                       )


CRYSPR_MULTIPLEX_REGULATION_LEVEL_ZERO = (CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_R1],
                                          CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_RN_MINUS_ONE],
                                          CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_NTERM],
                                          CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_ALL])


MANDATORY_DOMEST_ENZYMES = ('BsmBI', 'BsaI')
OPTIONAL_DOMEST_ENZYMES = ('BtgZI', 'BpiI')

AUTO_CRISPR_SELECTION_CHOICES = OrderedDict([('crispr_edit_dicot_1', 'Dicot, 1 target'),
                                             ('crispr_edit_dicot_2', 'Dicot, 2 targets'),
                                             ('crispr_edit_dicot_3', 'Dicot, 3 targets'),
                                             ('crispr_edit_dicot_4', 'Dicot, 4 targets'),
                                             ('crispr_edit_dicot_5', 'Dicot, 5 targets'),
                                             ('crispr_edit_dicot_6', 'Dicot, 6 targets'),
                                             ('crispr_edit_monocot_2', 'Monocot, 2 targets'),
                                             ('crispr_edit_monocot_3', 'Monocot, 3 targets'),
                                             ('crispr_edit_monocot_4', 'Monocot, 4 targets'),
                                             ('crispr_edit_monocot_5', 'Monocot, 5 targets'),
                                             ('crispr_edit_monocot_6', 'Monocot, 6 targets')])

AUTO_CRISPR_SELECTION_POSITIONS_DOMEST = OrderedDict([(EDIT_E1, ('ATTG', 'GTTT')),
                                                      (EDIT_E1_N_MINUS_ONE, ('ATTG', 'GTTT')),
                                                      (EDIT_E1_NTERM, ('ATTG', 'GTTT')),
                                                      (EDIT_E1B, ('GGCA', 'GTTT')),
                                                      (EDIT_E1B_N_MINUS_ONE, ('GGCA', 'GTTT'))])

                                              
PARTS_TO_ASSEMBLE = {'basic': [CATEGORIES[PROM_5UTR_NTAG],
                               CATEGORIES[CDS],
                               CATEGORIES[UTR3_TERM]],
                     'secreted': [CATEGORIES[PROM_5UTR_NTAG],
                                  CATEGORIES[CDS1],
                                  CATEGORIES[CDS2_CTAG],
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
                     'gRNA_monocot_reg': [CRYSPER_CATEGORIES[PROM_MONOCOT],
                                          CRYSPER_CATEGORIES[TARGET_MONOCOT],
                                          CRYSPER_CATEGORIES[TER_CRYSPER]],
                     'gRNA_dicot_reg': [CRYSPER_CATEGORIES[PROM_DICOT],
                                        CRYSPER_CATEGORIES[TARGET_DICOT],
                                        CRYSPER_CATEGORIES[TER_CRYSPER]],

                     'grna_CAS12a_dicot': [CRYSPER_CATEGORIES[PROM_CAS12],
                                           CRYSPER_CATEGORIES[TARGET_CAS12A],
                                           CRYSPER_CATEGORIES[CLEAVAGE_SIGNAL]],
                     'crispr_multiplexing': [TRNA,
                                             CRISPR_MULTIPLEXING_TARGET,
                                             SCAFFOLD],
                     'crispr_edit_dicot_1': [CRYSPER_CATEGORIES[PROM_DICOT],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1_NTERM]],
                     'crispr_edit_monocot_6': [CRYSPER_CATEGORIES[PROM_MONOCOT],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1B],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E3],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E4],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_N_MINUS_1],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM]],
                     'crispr_edit_dicot_5': [CRYSPER_CATEGORIES[PROM_DICOT],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E3],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E4_N_MINUS_ONE],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM]],

                     'crispr_edit_dicot_4': [CRYSPER_CATEGORIES[PROM_DICOT],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E3_N_MINUS_ONE],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM]],

                     'crispr_edit_monocot_5': [CRYSPER_CATEGORIES[PROM_MONOCOT],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1B],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E3],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E4_N_MINUS_ONE],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM]],

                     'crispr_edit_monocot_4': [CRYSPER_CATEGORIES[PROM_MONOCOT],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1B],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E3_N_MINUS_ONE],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM]],

                     'crispr_edit_monocot_3': [CRYSPER_CATEGORIES[PROM_MONOCOT],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1B],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2_N],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM]],

                     'crispr_edit_monocot_2': [CRYSPER_CATEGORIES[PROM_MONOCOT],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1B_N_MINUS_ONE],
                                               CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM]],

                     'crispr_edit_dicot_6': [CRYSPER_CATEGORIES[PROM_DICOT],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E3],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E4],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_N_MINUS_1],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM]
                                             ],

                     'crispr_edit_dicot_3': [CRYSPER_CATEGORIES[PROM_DICOT],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E2_N],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM]],

                     'crispr_edit_dicot_2': [CRYSPER_CATEGORIES[PROM_DICOT],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_E1_N_MINUS_ONE],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[EDIT_NTERM]],

                     'crispr_regulation_3': [CRYSPER_CATEGORIES[PROM_DICOT],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_R1],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_RN_MINUS_ONE],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_NTERM]],

                     'crispr_regulation_1': [CRYSPER_CATEGORIES[PROM_DICOT],
                                             CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[REG_ALL]],


                     'gene_disruption': [FUNGAL_CATEGORIES[REVERSE_MARKER],
                                         FUNGAL_CATEGORIES[FLANK_5UTR],
                                         FUNGAL_CATEGORIES[FORWARD_MARKER],
                                         FUNGAL_CATEGORIES[FLANK_3UTR]],

                      'Fungal TU': [FUNGAL_CATEGORIES[FUNGAL_PROM_5UTR],
                                    FUNGAL_CATEGORIES[FUNGAL_CDS],
                                    FUNGAL_CATEGORIES[FUNGAL_UTR3_TERM],
                                    ],
                                    
                      'Fungi gene_knock-out': [FUNGAL_CATEGORIES[REVERSE_MARKER],
                                              FUNGAL_CATEGORIES[FLANK_5UTR],
                                              FUNGAL_CATEGORIES[FORWARD_MARKER],
                                              FUNGAL_CATEGORIES[FLANK_3UTR]
                                              ]

                     }
UT_PREFIX = PARTS_TO_ASSEMBLE['basic'][0][1]
UT_SUFFIX = PARTS_TO_ASSEMBLE['basic'][-1][2]

SITE_A = UT_PREFIX
SITE_B = UT_SUFFIX
SITE_C = 'GTCA'


DOMESTICATION_VECTORS_IN_GB = ('pUPD', 'pUPD2')
LEVEL_MINUS_1_VECTOR = "pVD1"
DOMESTICATED_VECTOR = 'pUPD2'

DOMEST_VECTOR_PREFIX = 'CTCG'
DOMEST_VECTOR_SUFFIX = 'CTCA'

MOCLO_INCOMPATIBLE_RESISTANCES = ('Ampicillin', 'AmpR')
# PUPD_PREFIX = 'CTCG'

# added to create all oligos in domesticator
OLIGO_UNIVERSAL = 'GCGCCGTCTCG'
ASSEMBLED_SEQ = 'GB_UA'
DOMESTICATED_SEQ = 'GB_UD'
CRYSPER_SEQ = 'GB_UC'
EXPERIMENT_ID_PREFIX = 'GB_EXP'

# CATEGORY SBOL IMAGE CORRESPONDENCE
SBOL_IMAGES = {FORWARD: {PROM_5UTR_NTAG: 'prom_5utr_ntag.png',
                         PROX: 'prox.png',
                         CORE_5UTR: 'core_5utr.png',
                         PROM_5UTR: 'prom_5utrf.png',
                         DIST: 'dist.png',
                         DIST_PROX: 'dist_prox.png',
                         INTERACTION_ADAPTOR: 'interaction_adaptor.png',
                         CDS: 'cds.png',
                         CDS1_CDS2: 'cds1_cds2.png',
                         CDS2_CTAG: 'cds2_ctag.png',
                         CDS1: 'cds1.png',
                         NTAG: 'ntag.png',
                         CTAG: 'ctag.png',
                         UTR3_TERM: '3utr_term.png',
                         GOI: 'goi.png',
                         INT: 'int.png',
                         IOG: 'iog.png',
                         PROM_5UTR_MIR173: 'prom_5utr_mir173.png',
                         FGOI: 'fgoi.png',
                         FS5: '5fs.png',
                         TARGET: 'target.png',
                         FS3: '3fs.png',
                         OTHER_TYPE_NAME: 'other.png',
                         TARGET_DICOT: 'target_dicot.png',
                         TARGET_MONOCOT: 'target_monocot.png',
                         PROM_DICOT: 'prom_dicot.png',
                         PROM_MONOCOT: 'prom_monocot.png',
                         TER_CRYSPER: 'sgrna.png',
                         CRISPR_MULTIPLEXING_TARGET: 'target_dicot.png',
                         CAS12_LEVEL_MINUS_ONE_2X: 'cas12a_2x.png',
                         CAS12_LEVEL_MINUS_ONE_3X: 'cas12a_3x.png',
                         CAS12_LEVEL_MINUS_ONE_4X: 'cas12a_4x.png',
                         CAS12_LEVEL_MINUS_ONE_5X: 'cas12a_5x.png',
                         CAS12_LEVEL_MINUS_ONE_6X: 'cas12a_6x.png',
                         CLEAVAGE_SIGNAL: 'cleavage_signal.png',
                         EDIT_E1: "edit_e1.png",
                         EDIT_E2: "edit_e2.png",
                         EDIT_E3: "edit_e3.png",
                         EDIT_E4: "edit_e4.png",
                         EDIT_E3_N_MINUS_ONE: "edit_e3_n_minus_one.png",
                         EDIT_E4_N_MINUS_ONE: "edit_e4_n_minus_one.png",
                         EDIT_N_MINUS_1: "edit_n_minus_one.png",
                         EDIT_NTERM: "edit_nterm.png",
                         EDIT_E1_N_MINUS_ONE: "edit_e1_n_minus_one.png",
                         EDIT_E1_NTERM: "edit_e1_nterm.png",
                         EDIT_E2_N: "edit_e2_n.png",
                         EDIT_E1B: "edit_e1b.png",
                         EDIT_E1B_N_MINUS_ONE: "edit_e1b_n_minus_one.png",
                         FLANK_3UTR: 'flank_3utr.png',
                         FLANK_5UTR: 'flank_5utr.png',
                         FORWARD_MARKER: 'forward_marker.png',
                         FUNGAL_UTR3_TERM: 'fungal_3utr_term.png',
                         FUNGAL_CDS: 'fungal_cds.png',
                         FUNGAL_PROM_5UTR: 'fungal_prom_5utr.png',
                         PROM_CAS12: 'prom_cas12.png',
                         REG_ALL: 'reg_all.png',
                         REG_NTERM: 'reg_nterm.png',
                         REG_R1: 'reg_r1.png',
                         REG_RN_MINUS_ONE: 'reg_rn_minus_one.png',
                         REVERSE_MARKER: 'reverse_marker.png',
                         TER_CRYSPER: 'sgrna.png',
                         TARGET: 'target.png',
                         TARGET_CAS12A: 'target_cas12a.png'},
               REVERSE: {PROM_5UTR_NTAG: 'prom_5utr_ntag.rev.png',
                         PROX: 'prox.rev.png',
                         CORE_5UTR: 'core_5utr.rev.png',
                         PROM_5UTR: 'prom_5utrf.rev.png',
                         DIST: 'dist.rev.png',
                         DIST_PROX: 'dist_prox.rev.png',
                         INTERACTION_ADAPTOR: 'interaction_adaptor.rev.png',
                         CDS: 'cds.rev.png',
                         CDS1_CDS2: 'cds1_cds2.rev.png',
                         CDS2_CTAG: 'cds2_ctag.rev.png',
                         CDS1: 'cds1.rev.png',
                         NTAG: 'ntag.rev.png',
                         CTAG: 'ctag.rev.png',
                         UTR3_TERM: '3utr_term.rev.png',
                         GOI: 'goi.rev.png',
                         INT: 'int.rev.png',
                         IOG: 'iog.rev.png',
                         PROM_5UTR_MIR173: 'prom_5utr_mir173.rev.png',
                         FGOI: 'fgoiv.png',
                         FS5: '5fs.rev.png',
                         TARGET: 'targetv.png',
                         FS3: '3fs.rev.png',
                         OTHER_TYPE_NAME: 'other.rev.png',
                         TARGET_DICOT: 'target_dicot.rev.png',
                         TARGET_MONOCOT: 'target_monocot.rev.png',
                         PROM_DICOT: 'prom_dicot.rev.png',
                         PROM_MONOCOT: 'prom_monocot.rev.png',
                         TER_CRYSPER: 'sgrna.rev.png',
                         CRISPR_MULTIPLEXING_TARGET: 'target_dicot.png',
                         CAS12_LEVEL_MINUS_ONE_2X: 'cas12a_2x.rev.png',
                         CAS12_LEVEL_MINUS_ONE_3X: 'cas12a_3x.rev.png',
                         CAS12_LEVEL_MINUS_ONE_4X: 'cas12a_4x.rev.png',
                         CAS12_LEVEL_MINUS_ONE_5X: 'cas12a_5x.rev.png',
                         CAS12_LEVEL_MINUS_ONE_6X: 'cas12a_6x.rev.png',
                         CLEAVAGE_SIGNAL: 'cleavage_signal.rev.png',
                         EDIT_E1: "edit_e1.rev.png",
                         EDIT_E2: "edit_e2.rev.png",
                         EDIT_E3: "edit_e3.rev.png",
                         EDIT_E4: "edit_e4.rev.png",
                         EDIT_E3_N_MINUS_ONE: "edit_e3_n_minus_one.rev.png",
                         EDIT_E4_N_MINUS_ONE: "edit_e4_n_minus_one.rev.png",
                         EDIT_N_MINUS_1: "edit_n_minus_one.rev.png",
                         EDIT_NTERM: "edit_nterm.rev.png",
                         EDIT_E1_N_MINUS_ONE: "edit_e1_n_minus_one.rev.png",
                         EDIT_E1_NTERM: "edit_e1_nterm.rev.png",
                         EDIT_E2_N: "edit_e2_n.rev.png",
                         EDIT_E1B: "edit_e1b.rev.png",
                         EDIT_E1B_N_MINUS_ONE: "edit_e1b_n_minus_one.rev.png",
                         FLANK_3UTR: 'flank_3utr.rev.png',
                         FLANK_5UTR: 'flank_5utr.rev.png',
                         FORWARD_MARKER: 'forward_marker.rev.png',
                         FUNGAL_UTR3_TERM: 'fungal_3utr_term.rev.png',
                         FUNGAL_CDS: 'fungal_cds.rev.png',
                         FUNGAL_PROM_5UTR: 'fungal_prom_5utr.rev.png',
                         PROM_CAS12: 'prom_cas12.rev.png',
                         REG_ALL: 'reg_all.rev.png',
                         REG_NTERM: 'reg_nterm.rev.png',
                         REG_R1: 'reg_r1.rev.png',
                         REG_RN_MINUS_ONE: 'reg_rn_minus_one.rev.png',
                         REVERSE_MARKER: 'reverse_marker.rev.png',
                         TER_CRYSPER: 'sgrna.rev.png',
                         TARGET: 'target.rev.png',
                         TARGET_CAS12A: 'target_cas12a.png'}
               }
