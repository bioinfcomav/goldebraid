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


from django.shortcuts import render_to_response, redirect
from django.template.context import RequestContext
from django.forms.formsets import formset_factory
from django.core.context_processors import csrf
from django.db import transaction
from django.db.utils import IntegrityError
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required
from django import forms
from django.db.models import Q
from django.core.paginator import Paginator, EmptyPage, InvalidPage
from django.forms.models import modelformset_factory
from django.http.response import HttpResponseForbidden, HttpResponseBadRequest

from goldenbraid.forms.experiment import (ExperimentForm, ExperimentNumForm,
                                          ExperimentFeatureForm,
                                          ExperimentSubFeatureForm,
                                          ExperimentSearchForm,
                                          ExperimentExcelForm,
                                          ExperimentManagementForm,
    ExperimentGenericFileForm)
from goldenbraid.models import (Experiment, Count, Db, Dbxref, ExperimentPerm,
                                ExperimentPropNumeric, ExperimentPropText,
                                Feature, ExperimentFeature,
                                ExperimentPropImage, ExperimentSubFeature,
                                ExperimentPropExcel, ExperimentPropGenericFile,
                                Cvterm)
from goldenbraid.settings import EXPERIMENT_ID_PREFIX
from goldenbraid.tags import GOLDEN_DB, EXPERIMENT_TYPES, NUMERIC_TYPES
from sphinx.addnodes import desc


def experiment_view(request, uniquename):
    'The feature view'
    try:
        experiment = Experiment.objects.get(uniquename=uniquename)
    except Experiment.DoesNotExist:
        experiment = None

    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if experiment is None:
        return render_to_response('goldenbraid_info.html',
                                  {'title': 'Experiment not exist',
                                   'info': 'This experiment ({0}) does not exist in the database'.format(uniquename)},
                                  context_instance=RequestContext(request))
    if not request_data:
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
    else:
        form = ExperimentManagementForm(request_data)
        if form.is_valid():
            form_data = form.cleaned_data
            action = form_data['action']
            experiment = Experiment.objects.get(uniquename=form_data['experiment'])
            req_context = RequestContext(request)
            if action == 'delete':
                if request.user.is_staff or request.user == experiment.owner:
                    experiment.delete()
                    return render_to_response('goldenbraid_info.html',
                                              {'title': 'Experiment deleted',
                                               'info': 'Experiment Deleted'},
                                              context_instance=req_context)
                else:
                    info_txt = 'You are not allowed to delete this experiment'
                    return render_to_response('goldenbraid_info.html',
                                              {'title': 'Not Allowed',
                                               'info': info_txt},
                                              context_instance=req_context)

            elif action in 'make_public' or 'make_private':
                if request.user.is_staff:
                    if action == 'make_public' and not experiment.is_public:
                        experiment.is_public = True
                    elif action == 'make_private' and experiment.is_public:
                        experiment.is_public = False
                    else:
                        raise RuntimeError('bad conbinations of input request')
                    return render_to_response('experiment_template.html',
                                              {'experiment': experiment,
                                               'info': 'Feature modified'},
                                              context_instance=req_context)
                else:
                    info_text = 'You are not allowed to modify this feature'
                    return render_to_response('Goldenbraid_info.html',
                                              {'title': 'Not Allowed',
                                               'info': info_text},
                                              context_instance=req_context)

            else:
                return HttpResponseBadRequest()


def _add_experiment(form, user, feat_formset, subfeat_form,
                    numeric_formset=None, text_formset=None,
                    image_formset=None, excel_formset=None,
                    generic_file_formset=None, is_public=False):
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
                if not feat_id_in_exp:
                    continue
                try:
                    feat_id = feat_id_in_exp['feature']
                    feat = Feature.objects.get(feature_id=feat_id)
                except Feature.DoesNotExist:
                    feat = None
                if feat:
                    if not user.is_staff and (feat.owner != user and
                                              not feat.is_public):
                        msg = 'You can not add an experiment with features in '
                        msg = 'which you are not the owner'
                        raise RuntimeError(msg)
                    ExperimentFeature.objects.create(experiment=experiment,
                                                     feature=feat)

            for subfeat_uniquename in subfeat_form.cleaned_data['features']:
                try:
                    feat = Feature.objects.get(uniquename=subfeat_uniquename)
                except Feature.DoesNotExist:
                    feat = None
                if feat:
                    ExperimentSubFeature.objects.create(experiment=experiment,
                                                        feature=feat)
            if numeric_formset is not None:
                for numeric_prop in numeric_formset.cleaned_data:
                    if not numeric_prop:
                        continue
                    type_ = numeric_prop['type']
                    value = numeric_prop['value']
                    ExperimentPropNumeric.objects.create(experiment=experiment,
                                                         type=type_,
                                                         value=value)
            if text_formset is not None:

                text_props = text_formset.save(commit=False)
                print text_props
                for text_prop in text_props:
                    text_prop.experiment = experiment
                    text_prop.save()
            if image_formset is not None:
                image_props = image_formset.save(commit=False)
                for image_prop in image_props:
                    image_prop.experiment = experiment
                    image_prop.save()

            if generic_file_formset is not None:
                for generic_file_props in generic_file_formset.cleaned_data:
                    desc = generic_file_props['description']
                    file_ = generic_file_props['file']
                    if file_ is None and desc is None:
                        continue
                    ExperimentPropGenericFile.objects.create(experiment=experiment,
                                                             description=desc,
                                                             file=file_)

            if excel_formset is not None:
                for excel_formdata in excel_formset.cleaned_data:
                    if not excel_formdata:
                        continue
                    description = excel_formdata['description']
                    excel_file = excel_formdata['excel']
                    ExperimentPropExcel.objects.create(experiment=experiment,
                                                       description=description,
                                                       excel=excel_file)

    except (IntegrityError, RuntimeError):
        transaction.rollback()
        raise
    return experiment


@login_required
def add_experiment_view(request, exp_type=None):
    if exp_type is None:
        return _add_experiment_free(request)
    elif exp_type == 'SE_001':
        return _add_experiment_SE_001(request)


@login_required
def _add_experiment_SE_001(request):
    exp_type_name = 'SE_001'
    exp_type = Cvterm.objects.get(cv__name=EXPERIMENT_TYPES,
                                  name=exp_type_name)
    context = RequestContext(request)
    context.update(csrf(request))

    request_data = request.POST if request.method == 'POST' else None
    FeatFormset = formset_factory(ExperimentFeatureForm)
    NumericFormset = formset_factory(ExperimentNumForm, extra=8)
    ImageFormset = modelformset_factory(ExperimentPropImage,
                                        exclude=('experiment',))
    GenericFileFormset = formset_factory(ExperimentGenericFileForm)
    TextFormset = modelformset_factory(ExperimentPropText,
                                       exclude=('experiment',))
    ExcelFormset = formset_factory(ExperimentExcelForm)

    if request_data:
        form = ExperimentForm(request_data, instance=Experiment())
        feat_formset = FeatFormset(request_data, prefix='feature')
        subfeat_form = ExperimentSubFeatureForm(request_data)
        numeric_formset = NumericFormset(request_data, prefix='numeric')
        generic_file_formset = GenericFileFormset(request_data, request.FILES,
                                                  prefix='generic_file')
        image_formset = ImageFormset(request_data, request.FILES,
                                     prefix='image')
        text_formset = TextFormset(request_data, prefix='text')
        excel_formset = ExcelFormset(request_data, request.FILES,
                                     prefix='excel')
        if (form.is_valid() and feat_formset.is_valid() and
            subfeat_form.is_valid() and numeric_formset.is_valid() and
            image_formset.is_valid() and text_formset.is_valid() and
                generic_file_formset.is_valid() and excel_formset.is_valid()):
            print text_formset
            try:
                experiment = _add_experiment(form=form,
                                             feat_formset=feat_formset,
                                             subfeat_form=subfeat_form,
                                             user=request.user,
                                             numeric_formset=numeric_formset,
                                             image_formset=image_formset,
                                             text_formset=text_formset,
                                             generic_file_formset=generic_file_formset,
                                             excel_formset=excel_formset
                                             )
                print "valid"
            except IntegrityError as error:
                print error
                raise
            except RuntimeError as error:
                msg = 'This user is not entitled to add this experiment'
                msg += ' as it is'
                return HttpResponseForbidden(msg)
            return redirect(experiment.url)
        else:
            print numeric_formset.errors
            print "no valid"

    else:
        form = ExperimentForm(instance=Experiment())
        get_request_data = request.GET if request.method == 'GET' else None
        if get_request_data:
            initial = [{'feature': get_request_data.get('feature')}]
        else:
            initial = None
        feat_formset = FeatFormset(initial=initial, prefix='feature')
        subfeat_form = ExperimentSubFeatureForm()
        numeric_formset = NumericFormset(prefix='numeric')
        generic_file_formset = GenericFileFormset(prefix='generic_file')
        none_image_query = ExperimentPropImage.objects.none()
        image_formset = ImageFormset(prefix='image', queryset=none_image_query)
        text_formset = TextFormset(prefix='text',
                                   queryset=ExperimentPropText.objects.none())
        excel_formset = ExcelFormset(prefix='excel')

    context['form'] = form
    context['feature_formset'] = feat_formset
    context['subfeat_form'] = subfeat_form
    context['numeric_formset'] = numeric_formset
    context['generic_file_formset'] = generic_file_formset
    context['image_formset'] = image_formset
    context['text_formset'] = text_formset
    context['excel_formset'] = excel_formset
    context['exp_cv_type'] = exp_type
    context['plant_species'] = 'Nicotiana bentamiana'
    context['chassis'] = "Agroinfiltrated leaves"
    context['quantitative_outputs'] = Cvterm.objects.filter(cv__name=NUMERIC_TYPES,
                                                            name__icontains='Normalized Fluc/Rluc')
    template = 'experiment_add_{}.html'.format(exp_type_name)
    return render_to_response(template, context)


@login_required
def _add_experiment_free(request):
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
    GeneriFileFormset = modelformset_factory(ExperimentPropGenericFile,
                                             exclude=('experiment',))
    ExcelFormset = formset_factory(ExperimentExcelForm)
    if request_data:
        form = ExperimentForm(request_data, instance=Experiment())
        feat_formset = FeatFormset(request_data, prefix='feature')
        subfeat_form = ExperimentSubFeatureForm(request_data)
        numeric_formset = NumericFormset(request_data, prefix='numeric')
        text_formset = TextFormset(request_data, prefix='text')
        image_formset = ImageFormset(request_data, request.FILES,
                                     prefix='image')
        generic_file_formset = GeneriFileFormset(request_data, request.FILES,
                                                 prefix='generic_file')
        excel_formset = ExcelFormset(request_data, request.FILES,
                                     prefix='excel')
        # print request_data, request.FILES
        if (form.is_valid() and numeric_formset.is_valid() and
            text_formset.is_valid() and feat_formset.is_valid() and
                subfeat_form.is_valid() and excel_formset.is_valid() and
                generic_file_formset.is_valid()):
            print "valid"
            try:
                experiment = _add_experiment(form=form,
                                             numeric_formset=numeric_formset,
                                             text_formset=text_formset,
                                             image_formset=image_formset,
                                             excel_formset=excel_formset,
                                             feat_formset=feat_formset,
                                             subfeat_form=subfeat_form,
                                             generic_file_formset=generic_file_formset,
                                             user=request.user)
            except IntegrityError as error:
                print error
                raise
            except RuntimeError as error:
                msg = 'This user is not entitled to add this experiment'
                msg += ' as it is'
                return HttpResponseForbidden(msg)
            return redirect(experiment.url)
        else:
            print "no valid"
    else:
        form = ExperimentForm(instance=Experiment())
        get_request_data = request.GET if request.method == 'GET' else None
        if get_request_data:
            initial = [{'feature': get_request_data.get('feature')}]
        else:
            initial = None

        feat_formset = FeatFormset(initial=initial, prefix='feature')
#         feat_formset = FeatFormset(prefix='feature')
        subfeat_form = ExperimentSubFeatureForm()
        numeric_formset = NumericFormset(prefix='numeric')
        text_formset = TextFormset(prefix='text',
                                   queryset=ExperimentPropText.objects.none())
        none_image_query = ExperimentPropImage.objects.none()
        image_formset = ImageFormset(prefix='image', queryset=none_image_query)
        none_genericf_query = ExperimentPropGenericFile.objects.none()
        generic_file_formset = GeneriFileFormset(prefix='generic_file',
                                                 queryset=none_genericf_query)
        excel_formset = ExcelFormset(prefix='excel')
    context['form'] = form
    context['feature_formset'] = feat_formset
    context['subfeat_form'] = subfeat_form
    context['numeric_formset'] = numeric_formset
    context['text_formset'] = text_formset
    context['image_formset'] = image_formset
    context['excel_formset'] = excel_formset
    context['generic_file_formset'] = generic_file_formset

    template = 'experiment_add_template.html'
    return render_to_response(template, context)


def _build_experiment_query(criteria, user=None):
    query = Experiment.objects
    if 'name_or_description' in criteria and criteria['name_or_description']:
        text = criteria['name_or_description']
        name_criteria = (Q(uniquename__icontains=text) |
                         Q(description__icontains=text))
        query = query.filter(name_criteria)
    if 'chasis_1' in criteria and criteria['chasis_1']:
        query = query.filter(chasis_1__icontains=criteria['chasis_1'])
    if 'chasis_2' in criteria and criteria['chasis_2']:
        query = query.filter(chasis_2__icontains=criteria['chasis_2'])
    if 'experiment_type' in criteria and criteria['experiment_type']:
        query = query.filter(type_cvterm_id=criteria['experiment_type'])
    if 'feature' in criteria and criteria['feature']:
        query = query.filter(experimentfeature__feature__feature_id=criteria['feature'])

    if user.is_staff:
        if 'only_user' in criteria and criteria['only_user']:
            query = query.filter(experimentperm__owner=user)

    else:
        if 'only_user' in criteria and criteria['only_user']:
            query = query.filter(experimentperm__owner=user)
        else:
            query = query.filter(Q(experimentperm__is_public=True) |
                                 Q(experimentperm__owner=user))

    query = query.distinct()
    return query


def search_experiment(request):
    'The feature search view'

    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    template = 'search_experiment.html'
    content_type = None  # default
    if request_data:
        form = ExperimentSearchForm(request_data)
        if request.user.is_authenticated():
            _label = "Search only in my parts?"
            form.fields['only_user'] = forms.BooleanField(label=_label,
                                                          initial=False,
                                                          required=False)
        if form.is_valid():
            search_criteria = form.cleaned_data
            search_criteria = dict([(key, value) for key, value in
                                    search_criteria.items() if value])
            context['search_criteria'] = search_criteria
            experiment_queryset = _build_experiment_query(search_criteria,
                                                          user=request.user)
            if experiment_queryset:
                if experiment_queryset.count() == 1:
                    experiment_uniquename = experiment_queryset[0].uniquename
                    return redirect(experiment_view,
                                    uniquename=experiment_uniquename)
                paginator = Paginator(experiment_queryset, 25)
                # Make sure page request is an int. If not, deliver first page.
                try:
                    page_number = int(request.POST.get('page', '1'))
                except ValueError:
                    page_number = 1
                # If page request (9999) is out of range, deliver last page of
                # results.
                try:
                    page_object = paginator.page(page_number)
                except (EmptyPage, InvalidPage):
                    page_object = paginator.page(paginator.num_pages)
                context['experiments_page'] = page_object
            else:
                context['experiments_page'] = None

    else:
        form = ExperimentSearchForm()
        if request.user.is_authenticated():
            _label = "Search only in my parts?"
            form.fields['only_user'] = forms.BooleanField(label=_label,
                                                          initial=False,
                                                          required=False)

    context['form'] = form

    return render_to_response(template, context, content_type=content_type)
