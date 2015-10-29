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

import re
import os

from django.template.context import RequestContext
from django.core.context_processors import csrf
from django.shortcuts import render_to_response, redirect
from django.db.utils import IntegrityError
from django.http import HttpResponseServerError
from django.core.files import File
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required
from django.db import transaction
from django.http.response import HttpResponseBadRequest
from django.contrib.admin.views.decorators import staff_member_required
from django.core.exceptions import MultipleObjectsReturned
from django import forms
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.db.models import Q

from Bio import SeqIO

import django_tables2 as tables
from django_tables2 import RequestConfig
from django_tables2.utils import A

from goldenbraid.models import (Cvterm, Feature, Db, Dbxref, Featureprop,
                                FeaturePerm, FeatureRelationship)
from goldenbraid.tags import (GOLDEN_DB, VECTOR_TYPE_NAME,
                              DESCRIPTION_TYPE_NAME, ENZYME_IN_TYPE_NAME,
                              REFERENCE_TYPE_NAME, ENZYME_OUT_TYPE_NAME,
                              RESISTANCE_TYPE_NAME, DERIVES_FROM,
                              OTHER_TYPE_NAME, MODULE_TYPE_NAME, TU_TYPE_NAME)

from goldenbraid.forms.feature import (FeatureForm, FeatureManagementForm,
                                       get_all_vectors_as_choices, VectorForm,
                                       SearchFeatureForm,
                                       SPECIAL_SEARCH_CATEGORIES)
from goldenbraid.utils import get_prefix_and_suffix_index
from goldenbraid.settings import CATEGORIES, CRYSPER_CATEGORIES


def get_prefix_and_suffix(seq, enzyme):
    'it gets the prefix and the suffix of the feature seq'
    try:
        result = get_prefix_and_suffix_index(seq, enzyme)
        prefix_index, suffix_index, prefix_size = result
    except RuntimeError as error:
        if str(error) == "No rec_site":
            return None, None
        else:
            raise
    return _get_pref_suff_from_index(seq, prefix_index, suffix_index,
                                     prefix_size)


def _get_pref_suff_from_index(seq, prefix_index, suffix_index, prefix_size):

    prefix = seq[prefix_index:prefix_index + prefix_size]
    suffix = seq[suffix_index:suffix_index + prefix_size]

    if len(suffix) < prefix_size:
        remaining = prefix_size - len(suffix)
        suffix += seq[0:remaining]
    if len(prefix) < prefix_size:
        remaining = prefix_size - len(prefix)
        prefix += seq[0:remaining]

    return str(prefix), str(suffix)


def get_or_create_feature_relationship(object_, subject):
    derives_from = Cvterm.objects.get(name=DERIVES_FROM)
    try:
        FeatureRelationship.objects.get(subject=subject,
                                        type=derives_from, object=object_)
    except FeatureRelationship.DoesNotExist:
        FeatureRelationship.objects.create(subject=subject,
                                           type=derives_from, object=object_)


def add_relations(feature, seq):
    children = _parse_children_relations_from_gb(seq)
    if not children:
        return

    for child in children:
        child = Feature.objects.get(uniquename=child)
        get_or_create_feature_relationship(object_=feature, subject=child)


def _parse_children_relations_from_gb(seq):
    definition = seq.description
    if '(' in definition and ')' in definition:
        match = re.match('\((.+)\)', definition)
        return match.group(1).split(',')
    else:
        return None


def add_feature(name, type_name, vector, genbank, props, owner,
                is_public=False, prefix=None, suffix=None):
    'it adds a feature to the database'
    feature = None
    # transaction.set_autocommit(False)
    try:
        with transaction.atomic():
            seq = SeqIO.read(genbank, 'gb')
            residues = str(seq.seq)
            name = name
            uniquename = seq.id
            type_ = Cvterm.objects.get(name=type_name)
            db = Db.objects.get(name=GOLDEN_DB)
            genbank_file = File(genbank)
            try:
                dbxref = Dbxref.objects.create(db=db, accession=uniquename)
            except IntegrityError as error:
                raise IntegrityError('feature already in db' + str(error))
            vector_type = Cvterm.objects.get(name=VECTOR_TYPE_NAME)
            if vector and type_ == vector_type:
                # already checked in form validation
                raise RuntimeError("a vector feature can't  have a vector")
            if vector:
                vector = Feature.objects.get(uniquename=vector,
                                             type=vector_type)
            else:
                vector = None
            try:
                user = User.objects.get(username=owner)
            except User.DoesNotExist:
                raise RuntimeError('the given user does not exist')
            try:
                feature = Feature.objects.create(uniquename=uniquename,
                                                 name=name, type=type_,
                                                 residues=residues,
                                                 dbxref=dbxref, vector=vector,
                                                 genbank_file=genbank_file)

            except IntegrityError as error:
                raise IntegrityError('feature already in db: ' + str(error))

            FeaturePerm.objects.create(feature=feature, owner=user,
                                       is_public=is_public)

            feature.add_relations(seq)

            for type_name, values in props.items():
                try:
                    prop_type = Cvterm.objects.get(name=type_name)
                except Cvterm.DoesNotExist:
                    msg = 'Trying to add a property which cvterm does '
                    msg += 'not exist: {0}'
                    msg = msg.format(type_name)
                    raise RuntimeError(msg)
                except MultipleObjectsReturned:
                    for p in Cvterm.objects.filter(name=type_name):
                        print p.name
                        print p.cvterm_id
                        print p.definition
                    print "type_name", type_name
                    print "feature", feature.uniquename
                    raise
                rank = 0
                for value in values:
                    Featureprop.objects.create(feature=feature, type=prop_type,
                                               value=value, rank=rank)
                    rank += 1
            if suffix is None or prefix is None:
                if type_ == vector_type:
                    enzyme = props[ENZYME_IN_TYPE_NAME][0]
                else:
                    enzyme = vector.enzyme_out[0]
                prefix, suffix = get_prefix_and_suffix(residues, enzyme)
                if prefix is None or suffix is None:
                    raise RuntimeError('The given vector is not compatible with this part')

            if not _check_category(type_.name, prefix, suffix, vector):
                msg = 'It looks like yout construct does not match any of the '
                msg += 'standar GBCloning categories. You should uso Other '
                msg += 'category for this piece'
                raise RuntimeError(msg)
            feature.prefix = prefix.upper()
            feature.suffix = suffix.upper()
            feature.save()

    except (IntegrityError, RuntimeError):
        if feature:
            os.remove(feature.genbank_file.path)
        # transaction.rollback()
        raise
    return feature


def _check_category(type_name, prefix, suffix, vector):
    if type_name in (OTHER_TYPE_NAME, VECTOR_TYPE_NAME):
        return True

    for values in CATEGORIES.values():
        if values == (type_name, prefix, suffix):
            return True
    for values in CRYSPER_CATEGORIES.values():
        if values == (type_name, prefix, suffix):
            return True
    if (type_name, prefix, suffix) in ((TU_TYPE_NAME, 'GGAG', 'GTCA'),
                                       (TU_TYPE_NAME, 'GTCA', 'CGCT'),
                                       (MODULE_TYPE_NAME, 'GGAG', 'GTCA'),
                                       (MODULE_TYPE_NAME, 'GTCA', 'CGCT')):
        return True
    return False


def add_vector_from_form(form_data, user):
    'With this function we add a feature to the database'
    props = {}
    props[DESCRIPTION_TYPE_NAME] = [form_data['description']]
    props[ENZYME_IN_TYPE_NAME] = [form_data['enzyme_in']]
    props[ENZYME_OUT_TYPE_NAME] = [form_data['enzyme_out']]
    props[RESISTANCE_TYPE_NAME] = [form_data['resistance']]

    if form_data['reference']:
        props[REFERENCE_TYPE_NAME] = [form_data['reference']]

    vector = add_feature(name=form_data['name'], type_name=VECTOR_TYPE_NAME,
                         vector=None, genbank=form_data['gbfile'],
                         props=props, owner=user)

    return vector


def add_feature_from_form(form_data, user):
    'With this function we add a feature to the database'
    props = {}
    feature_type_name = form_data['type'].split(',')[0]
    if form_data['description']:
        props[DESCRIPTION_TYPE_NAME] = [form_data['description']]
    if form_data['reference']:
        props[REFERENCE_TYPE_NAME] = [form_data['reference']]

    feature = add_feature(name=form_data['name'], type_name=feature_type_name,
                          vector=form_data['vector'],
                          genbank=form_data['gbfile'],
                          props=props, owner=user)

    return feature


@staff_member_required
def add_vector_view(request):
    'The add feature view'
    context = RequestContext(request)
    context.update(csrf(request))
    request_data_post = request.POST if request.method == 'POST' else None
#    request_data_get = request.GET if request.method == 'GET' else None

    if request_data_post:
        form = VectorForm(request_data_post, request.FILES)
        if form.is_valid():
            vector_form_data = form.cleaned_data
            try:
                vector = add_vector_from_form(vector_form_data, request.user)
            except IntegrityError as error:
                if 'feature already in db' in str(error):
                    # TODO choose a template
                    req_context = RequestContext(request)
                    return render_to_response('feature_exists.html', {},
                                              context_instance=req_context)
                else:
                    return HttpResponseServerError(str(error))

            except Exception as error:
                req_context = RequestContext(request)
                return render_to_response('goldenbraid_info.html',
                                          {'title': 'Error',
                                           'info': str(error)},
                                          context_instance=req_context)
            # if everithing os fine we show the just added feature
            return redirect(vector.url)

    else:
        form = VectorForm()

    context['form'] = form
    template = 'vector_add_template.html'
    return render_to_response(template, context)


@login_required
def add_feature_view(request):
    'The add feature view'
    context = RequestContext(request)
    context.update(csrf(request))
    request_data_post = request.POST if request.method == 'POST' else None
#    request_data_get = request.GET if request.method == 'GET' else None

    if request_data_post:
        form = FeatureForm(request_data_post, request.FILES)
        if form.is_valid():
            feat_form_data = form.cleaned_data
            try:
                feature = add_feature_from_form(feat_form_data, request.user)
            except IntegrityError as error:
                if 'feature already in db' in str(error):
                    # TODO choose a template
                    req_context = RequestContext(request)
                    return render_to_response('feature_exists.html', {},
                                              context_instance=req_context)
                else:
                    return HttpResponseServerError(str(error))

            except Exception as error:
                req_context = RequestContext(request)
                return render_to_response('goldenbraid_info.html',
                                          {'title': 'Error',
                                           'info': str(error)},
                                          context_instance=req_context)
            # if everithing os fine we show the just added feature
            return redirect(feature.url)

    else:
        form = FeatureForm()
        _choices = get_all_vectors_as_choices(request.user)
        form.fields['vector'].widget.choices = _choices

    context['form'] = form
    template = 'feature_add_template.html'
    return render_to_response(template, context)


def feature_view(request, uniquename):
    'The feature view'
    try:
        feature = Feature.objects.get(uniquename=uniquename)
    except Feature.DoesNotExist:
        feature = None

    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if feature is None:
        info_text = 'This feature ({0}) does not exist in the database'
        return render_to_response('goldenbraid_info.html',
                                  {'title': 'Feature not exist',
                                   'info': info_text.format(uniquename)},
                                  context_instance=RequestContext(request))

    if not request_data:
        if (feature.is_public or
           (request.user.is_staff or request.user == feature.owner)):
            return render_to_response('feature_template.html',
                                      {'feature': feature},
                                      context_instance=RequestContext(request))
        else:
            info_text = 'You are not allowed to view this feature'
            return render_to_response('goldenbraid_info.html',
                                      {'title': 'Not Allowed',
                                       'info': info_text},
                                      context_instance=RequestContext(request))

    else:
        form = FeatureManagementForm(request_data)
        if form.is_valid():
            form_data = form.cleaned_data
            action = form_data['action']
            feature = Feature.objects.get(uniquename=form_data['feature'])
            req_context = RequestContext(request)
            if action == 'delete':
                if request.user.is_staff or request.user == feature.owner:
                    file_path = feature.genbank_file.path
                    feature.delete()
                    os.remove(file_path)
                    return render_to_response('goldenbraid_info.html',
                                              {'title': 'Feature deleted',
                                               'info': 'Feature Deleted'},
                                              context_instance=req_context)
                else:
                    info_txt = 'You are not allowed to delete this feature'
                    return render_to_response('goldenbraid_info.html',
                                              {'title': 'Not Allowed',
                                               'info': info_txt},
                                              context_instance=req_context)

            elif action in 'make_public' or 'make_private':

                if request.user.is_staff:
                    if action == 'make_public' and not feature.is_public:
                        featperm = FeaturePerm.objects.get(feature=feature)
                        featperm.is_public = True
                        featperm.save()
                    elif action == 'make_private' and feature.is_public:
                        featperm = FeaturePerm.objects.get(feature=feature)
                        featperm.is_public = False
                        featperm.save()
                    else:
                        raise RuntimeError('bad conbinations of input request')
                    return render_to_response('feature_template.html',
                                              {'feature': feature,
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


# search
def _build_name_or_prop_query(query, text, exact):
    'It looks in the name or in the feature property tables'
    if exact == 'True':
        name_criteria = Q(name__iexact=text) | Q(uniquename__iexact=text)
    else:
        name_criteria = (Q(name__icontains=text) |
                         Q(uniquename__icontains=text) |
                         Q(Q(featureprop__type__name=DESCRIPTION_TYPE_NAME) &
                           Q(featureprop__value__icontains=text)))

    query = query.filter(name_criteria)
    return query


def _build_feature_query(search_criteria, user):
    'Given a search criteria dict it returns a feature queryset'
    criteria = search_criteria
    query = Feature.objects
    if 'name_or_description' in criteria and criteria['name_or_description']:
        if 'name_exact' not in criteria:
            criteria['name_exact'] = 'False'
        query = _build_name_or_prop_query(query,
                                          criteria['name_or_description'],
                                          criteria['name_exact'])
    if 'category' in criteria and criteria['category']:
        kind, prefix, suffix = criteria['category'].split(',')
        query = query.filter(type__name=kind)
        if kind not in SPECIAL_SEARCH_CATEGORIES:
            query = query.filter(prefix=prefix, suffix=suffix)
    if user.is_staff:
        if 'only_user' in criteria and criteria['only_user']:
            query = query.filter(featureperm__owner__username=user)

    else:
        if 'only_user' in criteria and criteria['only_user']:
            query = query.filter(featureperm__owner__username=user)
        else:
            query = query.filter(Q(featureperm__is_public=True) |
                                 Q(featureperm__owner__username=user))

    query = query.distinct()
    return query


class FeatureTable(tables.Table):

    uniquename = tables.LinkColumn('feature_view', args=[A('uniquename')],
                                   verbose_name='Uniquename')
    gb_category_name = tables.Column(verbose_name='Type', order_by='type.name')
    description = tables.Column(verbose_name='Description', orderable=False)
    owner = tables.Column(verbose_name='Owner',
                          accessor='featureperm.owner')
    timecreation = tables.DateColumn(verbose_name='Creation Time', short=False)
    genbank_file = tables.FileColumn(verbose_name='Associated_genbank',
                                     orderable=False)

    class Meta:
        # model = Experiment
        attrs = {"class": "searchresult"}


def search_features_view(request):
    'The feature search view'

    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    template = 'search_feature.html'
    content_type = None  # default
    if request_data:
        form = SearchFeatureForm(request_data)
        if request.user.is_authenticated():
            usr_only_label = "Search only in my parts?"
            form.fields['only_user'] = forms.BooleanField(label=usr_only_label,
                                                          initial=False,
                                                          required=False)
        # _update_form_init_values(form, database)
        if form.is_valid():
            search_criteria = form.cleaned_data
            search_criteria = dict([(key, value)
                                    for key, value in search_criteria.items()
                                    if value])
            context['search_criteria'] = search_criteria
            feature_queryset = _build_feature_query(search_criteria,
                                                    user=request.user)
            download_search = search_criteria.get('download_search', False)
            if feature_queryset and download_search:
                context['features'] = feature_queryset
                template = 'search_feature_download.txt'
                content_type = 'text/plain'
            elif feature_queryset and not download_search:
                if feature_queryset.count() == 1:
                    feature_uniquename = feature_queryset[0].uniquename
                    return redirect(feature_view,
                                    uniquename=feature_uniquename)

                feature_table = FeatureTable(feature_queryset)
                RequestConfig(request).configure(feature_table)
                context['features'] = feature_table

            else:
                context['features'] = None
    else:
        form = SearchFeatureForm()
        if request.user.is_authenticated():
            only_usr_label = "Search only in my parts?"
            form.fields['only_user'] = forms.BooleanField(initial=False,
                                                          required=False,
                                                          label=only_usr_label)

    context['form'] = form

    return render_to_response(template, context, content_type=content_type)
