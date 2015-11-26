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

import operator

from django.shortcuts import render_to_response, redirect
from django.template.context import RequestContext
from django.forms.formsets import formset_factory
from django.core.context_processors import csrf
from django.db import transaction, connection
from django.db.utils import IntegrityError
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required
from django import forms
from django.db.models import Q
from django.forms.models import modelformset_factory
from django.http.response import HttpResponseForbidden, HttpResponseBadRequest
from django.utils.html import escape
from django.utils.safestring import mark_safe

import django_tables2 as tables
from django_tables2 import RequestConfig
from django_tables2.utils import A

from goldenbraid.forms.experiment import (ExperimentForm, ExperimentNumForm,
                                          ExperimentFeatureForm,
                                          ExperimentSubFeatureForm,
                                          ExperimentSearchForm,
                                          ExperimentExcelForm,
                                          ExperimentManagementForm,
                                          ExperimentGenericFileForm,
                                          BaseExperimentNumFormset,
                                          ExperimentKeywordForm,
                                          ExperimentProtocolForm)
from goldenbraid.models import (Experiment, Count, Db, Dbxref, ExperimentPerm,
                                ExperimentPropNumeric, ExperimentPropText,
                                Feature, ExperimentFeature,
                                ExperimentPropImage, ExperimentSubFeature,
                                ExperimentPropExcel, ExperimentPropGenericFile,
                                Cvterm, ExperimentKeyword)
from goldenbraid.settings import EXPERIMENT_ID_PREFIX
from goldenbraid.tags import GOLDEN_DB, EXPERIMENT_TYPES, NUMERIC_TYPES


EXP_SETTINGS = {'SE_001': {'plant_species': 'Nicotiana benthamiana',
                           'chassis': "Agroinfiltrated leaves",
                           'excel_mandatory': True},
                'SE_002': {'plant_species': 'Nicotiana benthamiana',
                           'chassis': "Agroinfiltrated leaves",
                           'quantitative_outputs_def': 'SE_001',
                           'excel_mandatory': True},
                'SE_003': {'excel_mandatory': False},
                'SE_004': {'plant_species': 'Nicotiana benthamiana',
                           'chassis': "Agroinfiltrated leaves",
                           'excel_mandatory': False},
                'SE_005': {'excel_mandatory': False}
                }


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
                if request.user.is_staff or request.user == experiment.owner:
                    if action == 'make_public' and not experiment.is_public:
                        experiment.is_public = True
                    elif action == 'make_private' and experiment.is_public:
                        experiment.is_public = False
                    else:
                        raise RuntimeError('bad conbinations of input request')
                    return render_to_response('experiment_template.html',
                                              {'experiment': experiment,
                                               'info': 'Experiment modified'},
                                              context_instance=req_context)
                else:
                    info_text = 'You are not allowed to modify this experiment'
                    return render_to_response('Goldenbraid_info.html',
                                              {'title': 'Not Allowed',
                                               'info': info_text},
                                              context_instance=req_context)

            else:
                return HttpResponseBadRequest()


def _add_experiment(form, user, feat_formset, subfeat_form,
                    numeric_formset=None, text_formset=None,
                    image_formset=None, excel_formset=None,
                    generic_file_formset=None, keyword_formset=None,
                    is_public=False, protocol_form=None):
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
            for keyword in keyword_formset.cleaned_data:
                if not keyword:
                    continue
                ExperimentKeyword.objects.create(experiment=experiment,
                                                 keyword=keyword['keyword'])

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
                    if not generic_file_props:
                        continue
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
            if protocol_form is not None:
  	        if protocol_form.cleaned_data['protocol'] is not None:
                    file_ = protocol_form.cleaned_data['protocol']
                    ExperimentPropGenericFile.objects.create(experiment=experiment,
                                                             description='Protocol',
                                                             file=file_)
                pass

    except (IntegrityError, RuntimeError):
        transaction.rollback()
        raise
    return experiment


@login_required
def add_experiment_view(request, exp_type):
    if exp_type == 'NS_000':
        return _add_experiment_free(request, exp_type)
    else:
        return _add_experiment_SE(request, exp_type)


@login_required
def _add_experiment_SE(request, exp_type_name):
    settings = EXP_SETTINGS[exp_type_name]
    exp_type = Cvterm.objects.get(cv__name=EXPERIMENT_TYPES,
                                  name=exp_type_name)
    context = RequestContext(request)
    context.update(csrf(request))
    quantitative_exp_def = settings.get('quantitative_outputs_def', exp_type)
    quantitative = Cvterm.objects.filter(cv__name=NUMERIC_TYPES,
                                         definition=quantitative_exp_def)
    request_data = request.POST if request.method == 'POST' else None
    FeatFormset = formset_factory(ExperimentFeatureForm)
    NumericFormset = formset_factory(ExperimentNumForm,
                                     extra=len(quantitative),
                                     formset=BaseExperimentNumFormset)
    ImageFormset = modelformset_factory(ExperimentPropImage,
                                        exclude=('experiment',))
    GenericFileFormset = formset_factory(ExperimentGenericFileForm)
    TextFormset = modelformset_factory(ExperimentPropText,
                                       exclude=('experiment',))
    ExcelFormset = formset_factory(ExperimentExcelForm)
    KeywordFormset = formset_factory(ExperimentKeywordForm)

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
        keyword_formset = KeywordFormset(request_data, prefix='keyword')
        protocol_form = ExperimentProtocolForm(request_data, request.FILES)

        if (form.is_valid() and feat_formset.is_valid() and
            subfeat_form.is_valid() and numeric_formset.is_valid() and
            image_formset.is_valid() and text_formset.is_valid() and
                generic_file_formset.is_valid() and excel_formset.is_valid()
                and keyword_formset.is_valid() and protocol_form.is_valid()):
            try:
                experiment = _add_experiment(form=form,
                                             feat_formset=feat_formset,
                                             subfeat_form=subfeat_form,
                                             user=request.user,
                                             numeric_formset=numeric_formset,
                                             image_formset=image_formset,
                                             text_formset=text_formset,
                                             generic_file_formset=generic_file_formset,
                                             excel_formset=excel_formset,
                                             keyword_formset=keyword_formset,
                                             protocol_form=protocol_form)
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
        keyword_formset = KeywordFormset(prefix='keyword')
        protocol_form = ExperimentProtocolForm()

    context['form'] = form
    context['feature_formset'] = feat_formset
    context['subfeat_form'] = subfeat_form
    context['numeric_formset'] = numeric_formset
    context['generic_file_formset'] = generic_file_formset
    context['image_formset'] = image_formset
    context['text_formset'] = text_formset
    context['excel_formset'] = excel_formset
    context['keyword_formset'] = keyword_formset
    context['protocol_form'] = protocol_form
    context['exp_cv_type'] = exp_type
    context['plant_species'] = settings.get('plant_species', None)
    context['chassis'] = settings.get('chassis', None)

    context['quantitative_outputs'] = quantitative.order_by('name')
    context['excel_mandatory'] = settings['excel_mandatory']
    return render_to_response('experiment_add_standard.html', context)


@login_required
def _add_experiment_free(request, exp_type_name):
    'The add feature view'
    exp_type_name = 'NS_000'
    context = RequestContext(request)
    context.update(csrf(request))

    exp_type = Cvterm.objects.get(cv__name=EXPERIMENT_TYPES,
                                  name=exp_type_name)

    request_data = request.POST if request.method == 'POST' else None
    FeatFormset = formset_factory(ExperimentFeatureForm)
    TextFormset = modelformset_factory(ExperimentPropText,
                                       exclude=('experiment',))
    ImageFormset = modelformset_factory(ExperimentPropImage,
                                        exclude=('experiment',))
    GeneriFileFormset = modelformset_factory(ExperimentPropGenericFile,
                                             exclude=('experiment',))
    ExcelFormset = formset_factory(ExperimentExcelForm)
    KeywordFormset = formset_factory(ExperimentKeywordForm)

    if request_data:
        form = ExperimentForm(request_data, instance=Experiment())
        feat_formset = FeatFormset(request_data, prefix='feature')
        subfeat_form = ExperimentSubFeatureForm(request_data)
        text_formset = TextFormset(request_data, prefix='text')
        image_formset = ImageFormset(request_data, request.FILES,
                                     prefix='image')
        generic_file_formset = GeneriFileFormset(request_data, request.FILES,
                                                 prefix='generic_file')
        excel_formset = ExcelFormset(request_data, request.FILES,
                                     prefix='excel')
        keyword_formset = KeywordFormset(request_data, prefix='keyword')
        if (form.is_valid() and text_formset.is_valid() and
            feat_formset.is_valid() and subfeat_form.is_valid() and
                excel_formset.is_valid() and generic_file_formset.is_valid()):
            try:
                experiment = _add_experiment(form=form,
                                             text_formset=text_formset,
                                             image_formset=image_formset,
                                             excel_formset=excel_formset,
                                             feat_formset=feat_formset,
                                             subfeat_form=subfeat_form,
                                             generic_file_formset=generic_file_formset,
                                             keyword_formset=keyword_formset,
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
        text_formset = TextFormset(prefix='text',
                                   queryset=ExperimentPropText.objects.none())
        none_image_query = ExperimentPropImage.objects.none()
        image_formset = ImageFormset(prefix='image', queryset=none_image_query)
        none_genericf_query = ExperimentPropGenericFile.objects.none()
        generic_file_formset = GeneriFileFormset(prefix='generic_file',
                                                 queryset=none_genericf_query)
        excel_formset = ExcelFormset(prefix='excel')
        keyword_formset = KeywordFormset(prefix='keyword')
    context['form'] = form
    context['feature_formset'] = feat_formset
    context['subfeat_form'] = subfeat_form
    context['text_formset'] = text_formset
    context['image_formset'] = image_formset
    context['excel_formset'] = excel_formset
    context['generic_file_formset'] = generic_file_formset
    context['keyword_formset'] = keyword_formset
    context['exp_cv_type'] = exp_type
    template = 'experiment_add_free.html'
    return render_to_response(template, context)


def _get_experiments_for_feature(feature_id):
    # It returns the experiments of the parents of the given feature
    sql_recursive = '''WITH RECURSIVE subparts (subpart, part) AS(
        SELECT subject_id as subpart, object_id as part
        FROM feature_relationship
        WHERE type_id=34
   UNION
        SELECT subparts.subpart, feature_relationship.object_id
        FROM feature_relationship, subparts
        WHERE subparts.part = feature_relationship.subject_id
)
SELECT experiment.experiment_id
FROM subparts, experimentfeature, experiment
where (experimentfeature.feature_id= subparts.part AND subparts.subpart=%s
AND experimentfeature.experiment_id=experiment.experiment_id);
        '''
    cursor = connection.cursor()
    cursor.execute(sql_recursive, [feature_id])
    return [e[0] for e in cursor.fetchall()]


def _build_experiment_query(criteria, user=None):
    query = Experiment.objects
    if 'feature' in criteria and criteria['feature']:
        # With this sql we get the experiment_ids where
        exp_ids = _get_experiments_for_feature(criteria['feature'])
        query = query.filter(Q(experiment_id__in=exp_ids) |
                             Q(experimentfeature__feature__feature_id=criteria['feature']))

    if 'name_or_description' in criteria and criteria['name_or_description']:
        text = criteria['name_or_description']
        name_criteria = (Q(uniquename__icontains=text) |
                         Q(description__icontains=text) |
                         Q(experimentkeyword__keyword__icontains=text))
        query = query.filter(name_criteria)
    if 'experiment_type' in criteria and criteria['experiment_type']:
        query = query.filter(type__name=criteria['experiment_type'])

    if 'numeric_types' in criteria and criteria['numeric_types']:
        ge = criteria['ge'] if 'ge' in criteria and criteria['ge'] else None
        le = criteria['le'] if 'le' in criteria and criteria['le'] else None
        num_criteria = []
        num_type_names = criteria['numeric_types']
        for num_type_name in num_type_names:
            num_type = Cvterm.objects.get(cv__name=NUMERIC_TYPES,
                                          name=num_type_name)
            each_criteria = [Q(experimentpropnumeric__type=num_type)]
            if ge:
                each_criteria.append(Q(experimentpropnumeric__value__gte=ge))
            if le:
                each_criteria.append(Q(experimentpropnumeric__value__lte=le))
            num_criteria.append(reduce(operator.and_, each_criteria))
        query = query.filter(reduce(operator.or_, num_criteria))

    if user.is_staff:
        if 'only_user' in criteria and criteria['only_user']:
            query = query.filter(experimentperm__owner=user.id)

    else:
        if 'only_user' in criteria and criteria['only_user']:
            query = query.filter(experimentperm__owner=user.id)
        else:
            query = query.filter(Q(experimentperm__is_public=True) |
                                 Q(experimentperm__owner=user.id))

    query = query.distinct()
    # print query.query
    return query


class ExperimentTable(tables.Table):

    uniquename = tables.LinkColumn('experiment_view', args=[A('uniquename')],
                                   verbose_name='Uniquename')
    type = tables.Column(verbose_name='Type', accessor='type.name')
    keywords = tables.Column(verbose_name='Keywords', orderable=False)
    features_used_in_experiment = tables.Column(verbose_name='GBelements',
                                                orderable=False)
    owner = tables.Column(verbose_name='Owner',
                          accessor='experimentperm.owner')
    timecreation = tables.DateColumn(verbose_name='Creation Time', short=False)

    def render_keywords(self, value):
        return ", ".join(value)

    def render_features_used_in_experiment(self, value):
        return mark_safe(", ".join([self._make_url(v.uniquename)
                                    for v in value]))

    def _make_url(self, value):
        return mark_safe("<a href='/feature/{0}'>{0}</a>".format(escape(value)))

    class Meta:
        attrs = {"class": "searchresult"}


def search_experiment(request):
    'The feature search view'

    context = RequestContext(request)
    context.update(csrf(request))
    getdata = False
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
        getdata = True
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

            download_search = request.GET.get('download_search', False)
            if download_search:
                context['experiments'] = experiment_queryset
                template = 'search_experiment_download.txt'
                content_type = 'text/plain'

            elif experiment_queryset:
                if experiment_queryset.count() == 1:
                    experiment_uniquename = experiment_queryset[0].uniquename
                    return redirect(experiment_view,
                                    uniquename=experiment_uniquename)
                else:
                    experiment_table = ExperimentTable(experiment_queryset,
                                                       template='table.html')
                    RequestConfig(request).configure(experiment_table)
                    experiment_table.paginate(page=request.GET.get('page', 1),
                                              per_page=25)
                    context['experiments'] = experiment_table
                    if not getdata:
                        context['criteria'] = "".join([';{}={}'.format(k, v)
                                          for k, v in search_criteria.items()])
            else:
                context['experiments'] = None
        else:
            print form.errors
            print "no_valid"

    else:
        form = ExperimentSearchForm()
        if request.user.is_authenticated():
            _label = "Search only in my parts?"
            form.fields['only_user'] = forms.BooleanField(label=_label,
                                                          initial=False,
                                                          required=False)

    context['form'] = form

    return render_to_response(template, context, content_type=content_type)
