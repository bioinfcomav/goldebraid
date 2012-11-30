from django.template.context import RequestContext
from django.core.context_processors import csrf
from django import forms
from django.shortcuts import render_to_response
from django.core.exceptions import ValidationError
from Bio import SeqIO

from goldenbraid.models import Cvterm, Feature, Db, Dbxref, Featureprop
from goldenbraid.settings import DB
from goldenbraid.tags import (GOLDEN_DB, VECTOR_TYPE_NAME,
                              DESCRIPTION_TYPE_NAME, ENZYME_IN_TYPE_NAME,
                              ENZYME_OUT_TYPE_NAME, ENZYME_TYPE_NAME)


class FeatureForm(forms.Form):
    'Form to add features to db'
    name = forms.CharField(max_length=255, required=False)
    description = forms.CharField(max_length=255, required=False)
    type = forms.CharField()
    vector = forms.CharField(required=False)
    enzyme_in = forms.CharField(required=False)
    # TODO we have to change the widgte to allow multiple values
    # Now we do spliting the text with a comma
    enzyme_out = forms.CharField(required=False)

    gbfile_label = 'Select a GenBank-formatted local file on your computer'
    gbfile = forms.FileField(label=gbfile_label, required=True)

#    def clean_uniquename(self):
#        'It checks that the unique name is unique in database'
#        uniquename = self.cleaned_data['uniquename']
#        try:
#            Feature.objects.using(DB).get(uniquename=uniquename)
#        except Feature.DoesNotExist:
#            return uniquename
#        raise ValidationError('There is already a feature with this uniquename')

    def clean_type(self):
        'It validates the type field'
        type_str = self.cleaned_data['type']
        try:
            Cvterm.objects.using(DB).get(name=type_str)
        except Cvterm.DoesNotExist:
            raise ValidationError('This type does not exist in the database')
        return type_str

    def clean_vector(self):
        '''It validates the vector.

        If the feature is a vector is doe not validate anything
        if featuer is not a vector if validates that the vector
        is in the database'''
        vector = self.cleaned_data['vector']
        error_in_type = self._errors.get('type', False)
        if error_in_type:
            return vector
        type_str = self.cleaned_data['type']

        if type_str != VECTOR_TYPE_NAME:
            try:
                vector_type = Cvterm.objects.using(DB).get(name=VECTOR_TYPE_NAME)
                Feature.objects.using(DB).get(uniquename=vector,
                                              type=vector_type)
            except Feature.DoesNotExist:
                raise ValidationError('The given vector does not exist')

        else:
            if vector:
                raise ValidationError('A vector does not have a vector')

        return vector

    def _validate_enzyme(self, kind):
        '''It validates the vector.

        If the feature is a vector is doe not validate anything
        if featuer is not a vector if validates that the vector
        is in the database'''
        # TODO. change how we deal with two enzymes for the same field
        enzymes = self.cleaned_data['enzyme_{}'.format(kind)].split(',')
        error_in_type = self._errors.get('type', False)
        error_in_vector = self._errors.get('vector', False)

        if error_in_type or error_in_vector:
            return enzymes

        type_ = self.cleaned_data['type']
        if type_ == VECTOR_TYPE_NAME:
            if enzymes[0] == u'':
                err = 'A vector must have a enzyme {}'.format(kind)
                raise ValidationError(err)

            enzyme_type = Cvterm.objects.using(DB).get(name=ENZYME_TYPE_NAME)
            for enzyme in enzymes:
                enzyme = enzyme.strip()
                try:
                    Feature.objects.using(DB).get(uniquename=enzyme,
                                                  type=enzyme_type)
                except Feature.DoesNotExist:
                    msg = 'The given enzyme {} does not exist'.format(enzyme)
                    raise ValidationError(msg)

        else:
            if enzymes:
                err = 'Only vectors have enzyme {}'.format(kind)
                raise ValidationError(err)

        return enzymes

    def clean_enzyme_in(self):
        return self._validate_enzyme('in')

    def clean_enzyme_out(self):
        return self._validate_enzyme('out')


def add_feature(form_data):
    'With this function we add a feature to the database'
    seq = SeqIO.read(form_data['gbfile'], 'gb')
    residues = str(seq.seq)
    name = form_data['name']
    uniquename = seq.id
    type_ = Cvterm.objects.using(DB).get(name=form_data['type'])
    db = Db.objects.using(DB).get(name=GOLDEN_DB)
    dbxref = Dbxref.objects.using(DB).create(db=db, accession=uniquename)
    vector = form_data['vector']
    vector_type = Cvterm.objects.using(DB).get(name=VECTOR_TYPE_NAME)
    if vector:
        vector = Feature.objects.using(DB).get(uniquename=vector,
                                               type=vector_type)
    else:
        vector = None

    feature = Feature.objects.using(DB).create(uniquename=uniquename,
                                               name=name, type=type_,
                                               residues=residues,
                                               dbxref=dbxref, vector=vector)

    props = {}
    if form_data['description']:
        props[DESCRIPTION_TYPE_NAME] = [form_data['description']]
    if type_ == vector_type:
        props[ENZYME_IN_TYPE_NAME] = form_data['enzyme_in']
        props[ENZYME_OUT_TYPE_NAME] = form_data['enzyme_out']
    for type_name, values in props.items():
        type_ = Cvterm.objects.using(DB).get(name=type_name)
        rank = 0
        for value in values:
            Featureprop.objects.using(DB).create(feature=feature, type=type_,
                                             value=value, rank=rank)
            rank += 1
    return feature


def add_feature_view(request):
    'The add feature view'
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    else:
        request_data = None

    if request_data:
        form = FeatureForm(request_data, request.FILES)
        if form.is_valid():
            feat_form_data = form.cleaned_data
            try:
                feature = add_feature(feat_form_data)
            except Exception as error:
                print 'error', error
                feature = None
            if feature:
                return render_to_response('feature_template.html',
                                          {'feature': feature},
                                          context_instance=RequestContext(request))

    else:
        form = FeatureForm()
    context['form'] = form
    template = 'feature_add_template.html'
    return render_to_response(template, context)


def feature_view(request, uniquename):
    'The feature view'
    try:
        feature = Feature.objects.using(DB).get(uniquename=uniquename)
    except Feature.DoesNotExist:
        feature = None
    return render_to_response('feature_template.html', {'feature': feature},
                              context_instance=RequestContext(request))

