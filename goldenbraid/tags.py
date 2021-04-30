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

MONOCOT_TAXA = "Monocot"
DICOT_TAXA = "Dicot"
LEVEL_MINUS_1 = '-1'
VECTOR_TYPE_NAME = 'Vector'
ENZYME_TYPE_NAME = 'Enzyme'
TU_TYPE_NAME = 'TU'
OTHER_TYPE_NAME = 'Other'
MODULE_TYPE_NAME = 'Module'
DESCRIPTION_TYPE_NAME = 'Description'
ENZYME_IN_TYPE_NAME = 'Enzyme_in'
ENZYME_OUT_TYPE_NAME = 'Enzyme_out'
RESISTANCE_TYPE_NAME = 'Resistance'
REFERENCE_TYPE_NAME = 'Reference'
GOLDEN_DB = 'goldenbraid'
GOLDEN_CV = 'goldenbraid'
FORWARD = 'forward'
REVERSE = 'reverse'
DERIVES_FROM = 'derives_from'

# Cvs that a user has to save in DB
EXPERIMENT_TYPES = 'experiment_types'
NUMERIC_TYPES = 'numeric_types'

# CVterms defining roles for features in expeiment
MAIN_ROLE = 'main_role'
ACCESORY_ROLE = 'Accesory_role'

# Category Names
PROM_5UTR_NTAG = 'PROM+5UTR (A1-A2-A3-B1-B2)'
PROX = 'PROX (A2)'
CORE_5UTR = 'CORE+5UTR (A3-B1-B2)'
PROM_5UTR = 'PROM+5UTR(f) (A1-A2-A3-B1)'
DIST = 'DIST(A1)'
DIST_PROX = 'DIST+PROX (A1-A2)'
INTERACTION_ADAPTOR = 'INTERACTION ADAPTOR (A1-A2-A3-B1-B2b)'
CDS = 'CDS (B3-B4-B5)'
CDS1_CDS2 = 'CDS1+CDS2 (B3-B4)'
CDS2_CTAG = 'CDS2+CTAG (B4-B5)'
CDS1 = 'CDS1 (B3)'
CDS2 = 'CDS2 (B4)'
NTAG = 'NTAG (B2)'
CTAG = 'CTAG (B5)'
UTR3_TERM = '3UTR+TERM (B6-C1)'
GOI = 'goi (B2-B3)'
INT = 'int (B4)'
IOG = 'iog (B5)'
PROM_5UTR_MIR173 = 'PROM+5UTR+mir173 (A1-A2-A3-B1b)'
FGOI = 'fgoi (B2-B3-B4-B5)'
FS5 = "5'FS (B2-B3b)"
TARGET = 'Target (B4b)'
FS3 = "'3FS (B5b)"
TARGET_DICOT = 'D Target (B3c-B4-B5c)'
TARGET_MONOCOT = 'M Target (B3d-B4-B5d)'
TARGET_CAS12A = "D CAS12a target (B3e-B4-B5e)"
GRNA_CAS12_DICOT = 'grna_CAS12a_dicot'
PROM_CAS12 = "PROM DPolIII+DRCas12 (A1-A2-A3-B1-B2e)"
CLEAVAGE_SIGNAL = "3_prime processing (B6c-C1)"
PROM_DICOT = 'PROM DPolIII (A1-A2-A3-B1-B2c)'
PROM_MONOCOT = 'PROM MPolIII (A1-A2-A3-B1-B2d)'
TER_CRYSPER = 'sgRNA (B6b-C1)'
CRISPR_EDITING = "CRISPR Multiplexing Editing"
CRISPR_REGULATION = "CRISPR Multiplexing Regulation"
CRISPR_MULTIPLEXING_TARGET = "CRISPR Multiplexing Target"

CAS12_LEVEL_MINUS_ONE_2X = "CRISPR Multiplexing CAS12A_2X (B3e-B4-B5-B6-C1)"
CAS12_LEVEL_MINUS_ONE_3X = "CRISPR Multiplexing CAS12A_3X (B3e-B4-B5-B6-C1)"
CAS12_LEVEL_MINUS_ONE_4X = "CRISPR Multiplexing CAS12A_4X (B3e-B4-B5-B6-C1)"
CAS12_LEVEL_MINUS_ONE_5X = "CRISPR Multiplexing CAS12A_5X (B3e-B4-B5-B6-C1)"
CAS12_LEVEL_MINUS_ONE_6X = "CRISPR Multiplexing CAS12A_6X (B3e-B4-B5-B6-C1)"
EDIT_E1 = "Multiplexing Dicot Edit (E1)"
EDIT_E2 = "Multiplexing Edit (E2)"
EDIT_E3 = "Multiplexing Edit (E3)"
EDIT_E4 = "Multiplexing Edit (E4)"
EDIT_E3_N_MINUS_ONE = "Multiplexing Edit (E3-E4-En-1)"
EDIT_E4_N_MINUS_ONE = "Multiplexing Edit (E4-En-1)"
EDIT_N_MINUS_1 = "Multiplexing Edit (En-1)"
EDIT_NTERM = "Multiplexing Edit (EnC1)"
EDIT_E1_N_MINUS_ONE = "Multiplexing Dicot Edit (E1-E2-E3-E4-En-1)"
EDIT_E1_NTERM = "Multiplexing Dicot Edit (E1-E2-E3-E4-En-1-EnC1)"
EDIT_E1B_NTERM = "Multiplexing Monocot Edit (E1B-E2-E3-E4-En-1-EnC1)"
EDIT_E2_N = "Multiplexing Edit (E2-E3-E4-En-1)"
EDIT_E1B = "Multiplexing Monocot Edit (E1b)"
EDIT_E1B_N_MINUS_ONE = "Multiplexing Monocot Edit (E1b-E2-E3-E4-En-1)"


REG_R1 = "Multiplexing Regulation R1"
REG_RN_MINUS_ONE = "Multiplexing Regulation RN-1"
REG_NTERM = "Multiplexing Regulation RnC-1"
REG_ALL = "Multiplexing Regulation (R1-Rn-1-RnC-1)"
IREG_R1 = "Multiplexing Inducible Regulation IR1"
IREG_RN_MINUS_ONE = "Multiplexing Inducible Regulation IRN-1"
IREG_NTERM = "Multiplexing Inducible Regulation IRnB5"
IREG_ALL = "Multiplexing Inducible Regulation (IR1-IRn-1-IRnB5)"
SCAFFOLD = "scaffold"
TRNA = "tRNA"



#FungalBraid

FUNGAL_PROM_5UTR = 'fungal_PROM+5UTR (A1-A2-A3-B1-B2)'
FUNGAL_CDS = 'fungal_CDS (B3-B4-B5)'
FUNGAL_UTR3_TERM = 'fungal_3UTR+TERM (B6-C1)'
REVERSE_MARKER = "negative_marker (A1-A2-A3)"
FORWARD_MARKER = "positive_marker (B3-B4-B5)"
FLANK_5UTR = "5_prime_flank (B1-B2)"
FLANK_3UTR = "3_prime flank (B6-C1)"
FUNGAL_TU_TYPE_NAME = "Fungal TU"
FUNGAL_KNOCK_OUT = "Fungi gene_knock-out"
FUNGAL_MODULE_TYPE_NAME = "Fungal Module"
