from django import forms
from django.core.exceptions import ValidationError
from django.forms.models import ModelForm


from goldenbraid.models import (Cvterm, Feature, Experiment,
                                ExperimentPropNumeric, ExperimentPropText, Cv,
                                ExperimentPropImage, ExperimentPropExcel)
from goldenbraid.tags import (EXPERIMENT_TYPES, NUMERIC_TYPES)
from goldenbraid.forms.widgets import (AutocompleteTextInput,
                                       DinamicSelectMultiple)
from django.forms.widgets import Select
from goldenbraid.excel import parse_xlsx
from zipfile import BadZipfile


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


class ExperimentNumForm(ModelForm):
    def __init__(self, *args, **kwargs):
        super(ExperimentNumForm, self).__init__(*args, **kwargs)
        cv = Cv.objects.get(name=NUMERIC_TYPES)
        exp_type_choices = [('', '')]
        for cvterm in Cvterm.objects.filter(cv=cv):
            exp_type_choices.append((cvterm.cvterm_id, cvterm.name))
        self.fields['type'] = forms.ChoiceField(choices=exp_type_choices)

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
    features = DynamicMultipleChoiceField(widget=_widget)

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


class ExperimentExcelForm_old(ModelForm):
    class Meta:
        model = ExperimentPropExcel
        exclude = ['experiment']

    def clean_excel(self):
        excel = self.cleaned_data['excel']
        print type(excel)


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


class ExperimentSearchForm(forms.Form):
    help_name = 'Accession or description'
    name_or_description = forms.CharField(max_length=100, required=False,
                                          label=help_name)
    chasis_1 = forms.CharField(max_length=100, required=False)
    chassis_2 = forms.CharField(max_length=100, required=False)

    _choices = [('', '')] + [(cvterm.name, cvterm.name) for cvterm in Cvterm.objects.filter(cv__name=EXPERIMENT_TYPES)]
    help_kind = 'Choose the type of the experiment'
    experiment_type = forms.CharField(max_length=200, label=help_kind,
                                      required=False,
                                      widget=Select(choices=_choices))
    feature = FeatureField(max_length=100, required=False,
                           widget=AutocompleteTextInput(source='/api/feature_uniquenames/',
                                                        min_length=1))
