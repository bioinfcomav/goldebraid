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

from django import forms
from django.forms.widgets import Select
from django.core.exceptions import ValidationError
from django.core.context_processors import csrf
from django.template.context import RequestContext
from django.shortcuts import render_to_response, redirect
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.db.models import Q

from goldenbraid.models import Cvterm, Feature
from goldenbraid.views.feature import feature_view
from goldenbraid import settings
from goldenbraid.tags import (DESCRIPTION_TYPE_NAME, VECTOR_TYPE_NAME,
                              TU_TYPE_NAME, MODULE_TYPE_NAME, OTHER_TYPE_NAME)
from goldenbraid.settings import CATEGORIES

SPECIAL_SEARCH_CATEGORIES = (VECTOR_TYPE_NAME, TU_TYPE_NAME, MODULE_TYPE_NAME,
                             OTHER_TYPE_NAME)


def _get_category_name(category):
    if category[0] in SPECIAL_SEARCH_CATEGORIES:
        return category[0]
    for name, category_def in CATEGORIES.items():
        if category == category_def:
            return name
    return '{0}: {1}'.format(category[0], ','.join(category))
    raise ValueError('The given category not in the CATEGORY dictionary')


def _prepare_feature_kind():
    'It prepares the feature kind select choices to put in the type widget'
    if settings.SEARCH_MENU_TYPE_CHOICES:
        feature_categories = [(kind, kind) for kind in settings.SEARCH_MENU_TYPE_CHOICES]
    else:
        categories = Feature.objects.distinct('type', 'suffix', 'prefix').values('type__name', 'prefix', 'suffix')
        # VECTOR, other is special. manually added
        feature_categories = []
        for special_category in SPECIAL_SEARCH_CATEGORIES:
            feature_categories.append(('{0},None,None'.format(special_category),
                                                              special_category))
        for dict_category in categories:
            if dict_category['type__name'] in SPECIAL_SEARCH_CATEGORIES:
                continue
            category = (dict_category['type__name'], dict_category['prefix'],
                        dict_category['suffix'])
            category_name = _get_category_name(category)
            feature_categories.append((','.join(category), category_name))

    feature_categories.insert(0, ('', ''))  # no kind
    return feature_categories


class SearchFeatureForm(forms.Form):

    help_name = 'Accession or name or description'
    name_or_description = forms.CharField(max_length=100, required=False,
                                          label=help_name)
    choices = _prepare_feature_kind()
    help_kind = 'Type of feature'
    category = forms.CharField(max_length=200, label=help_kind, required=False,
                               widget=Select(choices=choices))

    def xclean_kind(self):
        type_str = self.cleaned_data['kind']
        if not type_str:
            return type_str
        try:
            Cvterm.objects.get(name=type_str)
        except Cvterm.DoesNotExist:
            raise ValidationError('This type does not exist in the database')
        return type_str


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
            form.fields['only_user'] = forms.BooleanField(label="Search only in my parts?", initial=False,
                                                         required=False)
        # _update_form_init_values(form, database)
        if form.is_valid():
            search_criteria = form.cleaned_data
            search_criteria = dict([(key, value) for key, value in search_criteria.items() if value])
            context['search_criteria'] = search_criteria
            feature_queryset = _build_feature_query(search_criteria, user=request.user)
            download_search = search_criteria.get('download_search', False)
            if feature_queryset and download_search:
                context['queryset'] = feature_queryset
                template = 'search_feature_download.txt'
                content_type = 'text/plain'
            elif feature_queryset and not download_search:
                if feature_queryset.count() == 1:
                    feature_uniquename = feature_queryset[0].uniquename
                    return redirect(feature_view,
                                    uniquename=feature_uniquename)

                paginator = Paginator(feature_queryset, 25)
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

                context['features_page'] = page_object
            else:
                context['features_page'] = None
    else:
        form = SearchFeatureForm()
        if request.user.is_authenticated():
            form.fields['only_user'] = forms.BooleanField(initial=False,
                                                         required=False,
                                              label="Search only in my parts?")

    context['form'] = form

    return render_to_response(template, context, content_type=content_type)
