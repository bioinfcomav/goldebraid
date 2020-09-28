from django import forms
from django.forms.widgets import Select

from goldenbraid.tags import (MONOCOT_TAXA, DICOT_TAXA,
                              REG_R1, REG_RN_MINUS_ONE, REG_NTERM,
                              REG_ALL, EDIT_E2_N, VECTOR_TYPE_NAME,
                              PROM_DICOT, PROM_MONOCOT)
from goldenbraid.utils import has_rec_sites
from goldenbraid.settings import (MONOCOT_EDIT_POS, DICOT_EDIT_POS,
                                  CRYSPR_MULTIPLEX_MODES, AUTO_CRISPR_SELECTION_CHOICES,
                                  PARTS_TO_ASSEMBLE, CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE)
from goldenbraid.models import Feature
from goldenbraid.tags import (CRISPR_MULTIPLEXING_TARGET)
from django.db.models import Q
from django.core.exceptions import ValidationError

REGULATING_CATEGORIES = [REG_R1, REG_RN_MINUS_ONE, REG_NTERM, REG_ALL]


def _seq_is_dna(string):
    len_sum = sum([string.count(l.upper()) + string.count(l.lower()) for l
                                                                  in ('ATCG')])
    return False if len_sum != len(string) else True

def populate_location_form(form, taxa):
    if taxa == MONOCOT_TAXA:
        position_choices = [('', '')]
        for item in MONOCOT_EDIT_POS:
            position_choices.append((item, item))
    elif taxa == DICOT_TAXA:
        position_choices = [('', '')]
        for item in DICOT_EDIT_POS:
            position_choices.append((item, item))
    position_choices.append((EDIT_E2_N, EDIT_E2_N))
    form.fields['position'].widget.choices = position_choices
    form.fields['taxa'].widget.choices = [(taxa, taxa)]


def populate_regulation_location_form(form):
    position_choices = [('', '')]
    for item in REGULATING_CATEGORIES:
        position_choices.append((item, item))
    form.fields['position'].widget.choices = position_choices


class LocationForm(forms.Form):
    taxa = forms.CharField(max_length=100, required=True, widget=Select(),
                           label="Taxa selected")
    position = forms.CharField(max_length=100,
                               label="Choose a position",
                               required=True,
                               widget=Select())


class AutoMultiTUForm(forms.Form):
    tu_choices = [('', ''), 
                  ("GB3264", "pEGB3a1_nptII-Cas9"),
                  ("GB2234", "nptII-Cas9-DsRed_en_al")]
    transcription_unit = forms.CharField(max_length=100,
                                         label="Choose a transcriptional unit",
                                         required=True,
                                         widget=Select(choices=tu_choices))
    vector_choices = [('', '')]
    available_vectors = Feature.objects.filter(Q(type__name="Vector") & 
                                               Q(name__contains="omega1") & 
                                               Q(featureperm__owner__username="admin"))
    for vector in available_vectors:
        uniquename = vector.uniquename
        vector_name = vector.name
        vector_choices.append((uniquename, vector_name))
    vector = forms.CharField(max_length=100,
                             label="Choose a vector",
                             required=True,
                             widget=Select(choices=vector_choices))


class RegulationLocationForm(forms.Form):
    position = forms.CharField(max_length=100,
                               label="Choose a position",
                               required=True,
                               widget=Select())


def populate_sequence_form(form, taxa, position, user):
    form.fields['taxa'].widget.choices = [(taxa, taxa)]
    form.fields['position'].widget.choices = [(position, position)]
    targets = Feature.objects.filter(type__name=CRISPR_MULTIPLEXING_TARGET)
    if user.is_authenticated:
        targets = targets.filter(Q(featureperm__owner__username=user) |
                                 Q(featureperm__is_public=True))
    target_choices = [('', '')]
    for target in targets:
        target_name = "{} - {}".format(target.uniquename, target.name)
        target_choices.append((target.uniquename, target_name))
    form.fields['target'].widget.choices = target_choices


class AutoMultiTypeForm(forms.Form):
    multiplexing_type = forms.CharField(max_length=100, required=True, widget=Select(),
                                        label="Multiplexing type")


def populate_auto_type_form(form):
    choices = [('', '')]
    for category, name in AUTO_CRISPR_SELECTION_CHOICES.items():
        choices.append((category, name))
    form.fields['multiplexing_type'].widget.choices = choices


class SelectSequenceForm(forms.Form):

    taxa = forms.CharField(max_length=100, required=True, widget=Select(),
                           label="Taxa selected")
    position = forms.CharField(max_length=100, required=True, widget=Select(),
                               label="Position selected")
    target = forms.CharField(max_length=100, label="Choose a target Sequence",
                             widget=Select(), required=True)


def populate_regulation_sequence_form(form, position, user):
    form.fields['position'].widget.choices = [(position, position)]
    targets = Feature.objects.filter(type__name=CRISPR_MULTIPLEXING_TARGET)
    if user.is_authenticated:
        targets = targets.filter(Q(featureperm__owner__username=user) |
                                 Q(featureperm__is_public=True))
    target_choices = [('', '')]
    for target in targets:
        target_name = "{} - {}".format(target.uniquename, target.name)
        target_choices.append((target.uniquename, target_name))
    form.fields['target'].widget.choices = target_choices


class RegulationSequenceForm(forms.Form):

    position = forms.CharField(max_length=100, required=True, widget=Select(),
                               label="Position selected")
    target = forms.CharField(max_length=100, label="Choose a target Sequence",
                             widget=Select(), required=True)



class AutoMultiTargetsForm(forms.Form):

    def __init__(self, *args, **kwargs):
        category = kwargs.pop('category')
        super().__init__(*args, **kwargs)
        category_elements = PARTS_TO_ASSEMBLE[category]
        for element in category_elements:
            if element[0] in CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE:
                self.fields[element[0]] = forms.CharField(required=True)

    def clean(self):
        cleaned_data = self.cleaned_data
        for category, seq in cleaned_data.items():
            if len(seq) < 20:
                raise ValidationError('{} CRISPR target seq length must be 20 nucleotides'.format(category))

            if len(seq) != 20:
                msg = '{} CRISPR target seq length must be 20 nucleotides'.format(category)
                raise ValidationError(msg)

            if not _seq_is_dna(seq):
                raise ValidationError('{} Seq must only contain ACTG'.format(category))

            if has_rec_sites(seq):
                msg = '{} This sequence can not be domesticated.'
                msg += ' It has internal restriction sites'
                raise ValidationError(msg.format(category))

        return cleaned_data



class TargetTypeForm(forms.Form):
    choices = [('', '')]
    for item in CRYSPR_MULTIPLEX_MODES:
        choices.append((item, item))
    target_type = forms.CharField(max_length=100,
                                  label='Choose a target type',
                                  widget=Select(choices=choices),
                                  required=True)


class TaxaChoiceForm(forms.Form):
    choices = [('', '')]
    for item in (MONOCOT_TAXA, DICOT_TAXA):
        choices.append((item, item))
    taxa_choice = forms.CharField(max_length=100,
                                  label='Choose between monocot and dicot',
                                  widget=Select(choices=choices),
                                  required=True)

def populate_vector_select_form(form, user):
    if user.is_authenticated:
        vectors = Feature.objects.filter(type__name=VECTOR_TYPE_NAME)
    vector_choices = [('', '')]
    for vector in vectors:
        if vector.resistance[0] == 'Kanamycin':
            vector_name = "{} - {}".format(vector.uniquename, vector.name)
            vector_choices.append((vector, vector_name))
    form.fields['vector'].widget.choices = vector_choices


class VectorSelectForm(forms.Form):

    vector = forms.CharField(max_length=100, label='Select a Vector',
                             widget=Select(),
                             required=True)


def populate_promoter_choice_form(form, vector):
    form.fields['vector'].widget.choices = [(vector, vector)]
    choices = [('', '')]
    dicot_prom = Feature.objects.filter(Q(suffix='ATTG') & Q(prefix='GGAG')
                                        & Q(featureperm__owner__username='admin'))
    monocot_prom = Feature.objects.filter(Q(suffix='GGCA') & Q(prefix='GGAG')
                                          & Q(featureperm__owner__username='admin'))
    if len(dicot_prom) != 1 or len(monocot_prom) != 1:
        msg = "The number of availables promoters"
        msg += " for each kind is different than one"
        msg += " Dicot: {}, Monocot: {}".format(len(dicot_prom),
                                                len(monocot_prom))
        raise RuntimeError(msg)
    for choice in (dicot_prom[0], monocot_prom[0]):
        choices.append((choice.type, choice.type))
    form.fields['promoter'].widget.choices = choices


class PromoterChoiceForm(forms.Form):
    vector = forms.CharField(max_length=100, label='Select a Vector',
                             widget=Select(), required=True)
    promoter = forms.CharField(max_length=100, label='Select a Promoter',
                               widget=Select(), required=True)


def populate_configuration_choice_form(form, vector, promoter):
    form.fields['vector'].widget.choices = [(vector, vector)]
    form.fields['promoter'].widget.choices = [(promoter, promoter)]


class ConfigurationChoiceForm(forms.Form):
    vector = forms.CharField(max_length=100, label='Select a Vector',
                             widget=Select(), required=True)
    promoter = forms.CharField(max_length=100, label='Select a Promoter',
                               widget=Select(), required=True)
    configuration = forms.CharField(max_length=100, label='Select a Configuration',
                                    widget=Select(), required=True)
