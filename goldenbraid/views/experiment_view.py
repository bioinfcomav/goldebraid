'''
Created on 2015 mar. 13

@author: peio
'''

from django.shortcuts import render_to_response, redirect
from django.template.context import RequestContext
from django.forms.formsets import formset_factory
from django.core.context_processors import csrf
from django.db import transaction
from django.db.utils import IntegrityError
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required


from goldenbraid.forms import (ExperimentForm, ExperimentNumForm,
                               ExperimentTextForm, ExperimentFeatureForm)
from goldenbraid.models import (Experiment, Count, Db, Dbxref, ExperimentPerm,
                                ExperimentPropNumeric, ExperimentPropText,
    Cvterm, Cv, ExperimentPropImage, ExperimentFeature, Feature)
from goldenbraid.settings import EXPERIMENT_ID_PREFIX
from goldenbraid.tags import GOLDEN_DB
from django.forms.models import modelformset_factory
from django.contrib.admin.filters import ChoicesFieldListFilter


def experiment_view(request, uniquename):
    'The feature view'
    try:
        experiment = Experiment.objects.get(uniquename=uniquename)
    except Experiment.DoesNotExist:
        experiment = None

    if experiment is None:
        return render_to_response('goldenbraid_info.html',
                                  {'title': 'Experiment not exist',
                                   'info': 'This experiment ({0}) does not exist in the database'.format(uniquename)},
                                  context_instance=RequestContext(request))
    else:
        if (experiment.is_public or
           (request.user.is_staff or request.user == experiment.owner)):
            return render_to_response('experiment_template.html',
                                      {'experiment': experiment},
                                      context_instance=RequestContext(request))
        else:
            return render_to_response('goldenbraid_info.html',
                                      {'title': 'Not Allowed',
                                       'info': 'You are not allowed to view this experiment'},
                                      context_instance=RequestContext(request))


def _add_experiment(form, numeric_formset, text_formset, image_formset,
                    feat_formset, user, is_public=False):
    try:
        with transaction.atomic():
            experiment = form.save(commit=False)
            try:
                count = Count.objects.get(name=EXPERIMENT_ID_PREFIX)
            except Count.DoesNotExist:
                count = Count.objects.create(name=EXPERIMENT_ID_PREFIX,
                                             value=1)
            next_value = count.next
            uniquename = EXPERIMENT_ID_PREFIX + '_' + next_value
            experiment.uniquename = uniquename
            db = Db.objects.get(name=GOLDEN_DB)
            try:
                dbxref = Dbxref.objects.create(db=db, accession=uniquename)
            except IntegrityError as error:
                raise IntegrityError('feature already in db' + str(error))
            experiment.dbxref = dbxref
            experiment.save()

            try:
                user = User.objects.get(username=user)
            except User.DoesNotExist:
                raise RuntimeError('the given user does not exist')
            ExperimentPerm.objects.create(experiment=experiment, owner=user,
                                          is_public=is_public)
            for feat_id_in_exp in feat_formset.cleaned_data:
                try:
                    feat_id = feat_id_in_exp['feature']
                    feat = Feature.objects.get(feature_id=feat_id)
                except Feature.DoesNotExist:
                    feat = None
                if feat:
                    ExperimentFeature.objects.create(experiment=experiment,
                                                     feature=feat)

            for numeric_prop in numeric_formset.cleaned_data:
                if not numeric_prop:
                    continue
                type_ = numeric_prop['type']
                value = numeric_prop['value']
                ExperimentPropNumeric.objects.create(experiment=experiment,
                                                     type=type_, value=value)

            text_props = text_formset.save(commit=False)
            for text_prop in text_props:
                text_prop.experiment = experiment
                text_prop.save()

            image_props = image_formset.save(commit=False)
            for image_prop in image_props:
                image_prop.experiment = experiment
                image_prop.save()

    except (IntegrityError, RuntimeError):
        transaction.rollback()
        raise
    return experiment


@login_required
def add_experiment_view(request):
    'The add feature view'
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST if request.method == 'POST' else None
    FeatFormset = formset_factory(ExperimentFeatureForm)
    NumericFormset = formset_factory(ExperimentNumForm)
    TextFormset = modelformset_factory(ExperimentPropText,
                                       exclude=('experiment',))
    ImageFormset = modelformset_factory(ExperimentPropImage,
                                        exclude=('experiment',))
    if request_data:
        form = ExperimentForm(request_data, instance=Experiment())
        feat_formset = FeatFormset(request_data, prefix='feature')
        numeric_formset = NumericFormset(request_data, prefix='numeric')
        text_formset = TextFormset(request_data, prefix='text')
        image_formset = ImageFormset(request_data, request.FILES,
                                     prefix='image')
        print request_data
        if (form.is_valid() and numeric_formset.is_valid()
           and text_formset.is_valid() and feat_formset.is_valid()):
            print "valid"
            try:
                experiment = _add_experiment(form, numeric_formset,
                                             text_formset, image_formset,
                                             feat_formset,
                                             user=request.user)
            except IntegrityError as error:
                print error
                raise
            return redirect(experiment.url)
        else:
            print "no valid"
    else:
        form = ExperimentForm(instance=Experiment())
        feat_formset = FeatFormset(prefix='feature')
        numeric_formset = NumericFormset(prefix='numeric')
        text_formset = TextFormset(prefix='text',
                                   queryset=ExperimentPropText.objects.none())
        image_formset = ImageFormset(prefix='image',
                                     queryset=ExperimentPropImage.objects.none())

    context['form'] = form
    context['feature_formset'] = feat_formset
    context['numeric_formset'] = numeric_formset
    context['text_formset'] = text_formset
    context['image_formset'] = image_formset
    template = 'experiment_add_template.html'
    return render_to_response(template, context)
