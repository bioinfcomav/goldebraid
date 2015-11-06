from django import forms
from django.core.exceptions import ValidationError
from django.forms.models import ModelForm


from goldenbraid.models import (Cvterm, Feature, Experiment,
                                ExperimentPropNumeric, ExperimentPropText, Cv,
                                ExperimentPropImage)
from goldenbraid.tags import (EXPERIMENT_TYPES, NUMERIC_TYPES)
from goldenbraid.forms.widgets import (AutocompleteTextInput,
                                       DinamicSelectMultiple)
from django.forms.widgets import Select
from goldenbraid.excel import parse_xlsx
from zipfile import BadZipfile
from django.forms.formsets import BaseFormSet


class ExperimentForm(ModelForm):
    def __init__(self, *args, **kwargs):
        super(ExperimentForm, self).__init__(*args, **kwargs)
        cv = Cv.objects.get(name=EXPERIMENT_TYPES)
        exp_type_choices = [('', '')]
        for cvterm in Cvterm.objects.filter(cv=cv):
            exp_type_choices.append((cvterm.cvterm_id, cvterm.name))
        self.fields['type'] = forms.ChoiceField(choices=exp_type_choices)

    class Meta:
        model = Experiment
        fields = ['chasis_1', 'chasis_2', 'description', 'type']

    def clean_type(self):
        cvterm_id = self.cleaned_data['type']
        try:
            cvterm = Cvterm.objects.get(cvterm_id=cvterm_id)
        except Cvterm.DoesNotExist:
            raise ValidationError('This experiment type does not exist!')
        if cvterm.cv.name != EXPERIMENT_TYPES:
            raise ValidationError('This type is not an experiment type')
        return cvterm


class BaseExperimentNumFormset(BaseFormSet):
    def clean(self):
        if any(self.errors):
            # Don't bother validating the formset unless each form is valid on its own
            return
        values = []
        for form in self.forms:
            values.append(form.cleaned_data['value'])
        if all([True if v is None else False for v in values]):
            raise ValidationError('At least one quantitative value is needed')


class ExperimentNumForm(ModelForm):
    def __init__(self, *args, **kwargs):
        super(ExperimentNumForm, self).__init__(*args, **kwargs)
        cv = Cv.objects.get(name=NUMERIC_TYPES)
        exp_type_choices = [('', '')]
        for cvterm in Cvterm.objects.filter(cv=cv):
            exp_type_choices.append((cvterm.cvterm_id, cvterm.name))
        self.fields['type'] = forms.ChoiceField(choices=exp_type_choices)
        self.fields['value'].required = False

    class Meta:
        model = ExperimentPropNumeric
        exclude = ['experiment']

    def clean_type(self):
        cvterm_id = self.cleaned_data['type']
        try:
            cvterm = Cvterm.objects.get(cvterm_id=cvterm_id)
        except Cvterm.DoesNotExist:
            raise ValidationError('This experiment type does not exist!')
        if cvterm.cv.name != NUMERIC_TYPES:
            raise ValidationError('This type is not an experiment type')
        return cvterm


class FeatureField(forms.CharField):
    '''A specialized Field that validates the feature type given'''

    def to_python(self, value):
        'It transforms the value into a cvterm'

        if not value:
            return ''
        elif value.isdigit():
            return self._search_item_id(value, 'id')
        else:
            return self._search_item_id(value, 'uniquename')

    def _search_item_id(self, value, kind):
        'It returns the featureSrc given the name or id'
        try:
            if kind == 'id':
                feature = Feature.objects.get(feature_id=value)
            else:
                feature = Feature.objects.get(uniquename=value)
        except Feature.DoesNotExist:
            raise ValidationError('feature does not exists: {}'.format(value))
        return str(feature.feature_id)


class ExperimentFeatureForm(forms.Form):
    feature = FeatureField(max_length=100,
                           widget=AutocompleteTextInput(source='/api/feature_uniquenames/',
                                                        min_length=1))


class DynamicMultipleChoiceField(forms.MultipleChoiceField):
    def valid_value(self, value):
        return True


class ExperimentSubFeatureForm(forms.Form):
    _widget = DinamicSelectMultiple(source='/api/features_key_elements/',
                                    parent_class='ui-autocomplete-input')
    features = DynamicMultipleChoiceField(widget=_widget,
                                          label='GB Elements')

    def clean_features(self):
        feature_uniquenames = self.cleaned_data['features']
        feats = []
        for feature_uniquename in feature_uniquenames:
            try:
                feat = Feature.objects.get(uniquename=feature_uniquename)
                feats.append(feat)
            except Feature.DoesNotExist:
                msg = 'feature does not exists: {}'.format(feature_uniquename)
                raise ValidationError(msg)
        return feats


class ExperimentTextForm(ModelForm):
    class Meta:
        model = ExperimentPropText
        exclude = ['experiment']


class ExperimentImageForm(ModelForm):
    class Meta:
        model = ExperimentPropImage
        exclude = ['experiment']


class ExperimentGenericFileForm(forms.Form):
    description = forms.CharField(max_length=100)
    file = forms.FileField(required=False)

    def clean(self):
        cleaned_data = super(ExperimentGenericFileForm, self).clean()
        description = cleaned_data.get('description', None)
        file_ = cleaned_data.get('file', None)

        if description == 'protocol' and file_ is None:
            self.cleaned_data['description'] = None


class ExperimentExcelForm(forms.Form):
    description = forms.CharField(max_length=100)
    excel = forms.FileField()

    def clean_excel(self):
        excel = self.cleaned_data['excel']
        try:
            parse_xlsx(excel)
        except RuntimeError:
            raise ValidationError('excel file is malformed')
        except BadZipfile:
            raise ValidationError('excel file is not an excel(xlsx) file')
        return excel


class ExperimentKeywordForm(forms.Form):
    keyword = forms.CharField(max_length=50,
                              widget=AutocompleteTextInput(source='/api/exp_keywords/',
                                                           min_length=1,
                                                           force_check=False))


class ExperimentProtocolForm(forms.Form):
    protocol = forms.FileField(label='Upload protocol')


def _get_numeric_choices():
    choices = []
    cv = Cv.objects.get(name=NUMERIC_TYPES)
    for num_types in Cvterm.objects.filter(cv=cv):
        name = num_types.name
        choices.append((name, name))
    return choices


class ExperimentSearchForm(forms.Form):
    help_name = 'Accession or description'
    name_or_description = forms.CharField(max_length=100, required=False,
                                          label=help_name)
    exp_types = Cvterm.objects.filter(cv__name=EXPERIMENT_TYPES).order_by('name')
    _choices = [('', '')] + [(cvterm.name, cvterm.name) for cvterm in exp_types]
    help_kind = 'Choose the type of the experiment'
    experiment_type = forms.CharField(max_length=200, label=help_kind,
                                      required=False,
                                      widget=Select(choices=_choices))
    feature = FeatureField(max_length=100, required=False,
                           widget=AutocompleteTextInput(source='/api/feature_uniquenames/',
                                                        min_length=1))
    num_choices = _get_numeric_choices()
    numeric_types = DynamicMultipleChoiceField(required=False,
              widget=forms.widgets.CheckboxSelectMultiple(choices=num_choices))
    ge = forms.FloatField(required=False)
    le = forms.FloatField(required=False)

    def clean_numeric_types(self):
        num_types = self.cleaned_data['numeric_types']
        return num_types

    def clean(self):
        cleaned_data = super(ExperimentSearchForm, self).clean()
        numeric_types = cleaned_data.get('numeric_types', [])
        ge = cleaned_data.get('ge', None)
        le = cleaned_data.get('le', None)
        if not numeric_types and (ge is not None or le is not None):
            raise ValidationError('You need the numeric type')


class ExperimentManagementForm(forms.Form):
    experiment = forms.CharField(max_length=30, widget=forms.HiddenInput())
    action = forms.CharField(max_length=30, widget=forms.HiddenInput())

    def clean_action(self):
        action = self.cleaned_data['action']
        if action in ('delete', 'make_public', 'make_private'):
            return action
        raise ValidationError('action must be delete or make_public')

#     def clean_experiment(self):
#
#         return create_feature_validator('feature')(self)
