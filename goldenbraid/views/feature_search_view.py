from django import forms
from django.forms.widgets import Select
from django.core.exceptions import ValidationError
from django.core.context_processors import csrf
from django.template.context import RequestContext
from django.shortcuts import render_to_response, redirect
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.db.models import Q

from goldenbraid.models import Cvterm, Feature
from goldenbraid.settings import DB
from goldenbraid.views.feature_views import feature_view
from goldenbraid import settings


def _prepare_feature_kind(database):
    'It prepares the feature kind select choices to put in the type widget'
    if settings.SEARCH_MENU_TYPE_CHOICES:
        feature_kinds = [(kind, kind) for kind in settings.SEARCH_MENU_TYPE_CHOICES]
    else:
        kinds = Feature.objects.using(database).distinct('type').values('type__name')
        kinds = [kind['type__name'] for kind in kinds]
        feature_kinds = [(kind, kind.replace('_', ' ')) for kind in kinds]

    feature_kinds.insert(0, ('', ''))  # no kind
    return feature_kinds


class SearchFeatureForm(forms.Form):

    help_name = 'Accession or name or description'
    name_or_description = forms.CharField(max_length=100, required=False,
                                          label=help_name)
    choices = _prepare_feature_kind(DB)
    help_kind = 'Type of feature'
    kind = forms.CharField(max_length=200, label=help_kind, required=False,
                           widget=Select(choices=choices))

    def clean_kind(self):
        type_str = self.cleaned_data['kind']
        if not type_str:
            return type_str
        try:
            Cvterm.objects.using(DB).get(name=type_str)
        except Cvterm.DoesNotExist:
            raise ValidationError('This type does not exist in the database')
        return type_str


def _build_name_or_prop_query(query, text, exact):
    'It looks in the name or in the feature property tables'
    if exact == 'True':
        name_criteria = Q(name__iexact=text)
    else:
        name_criteria = Q(name__icontains=text)

    # prop_criteria_value = Q(featureprop__value__icontains=text)
    query = query.filter(name_criteria)  # | prop_criteria_value)
    return query


def _build_feature_query(search_criteria):
    'Given a search criteria dict it returns a feature queryset'
    criteria = search_criteria
    query = Feature.objects.using(DB)
    if 'name_or_description' in criteria and criteria['name_or_description']:
        if 'name_exact' not in criteria:
            criteria['name_exact'] = 'False'
        query = _build_name_or_prop_query(query,
                                          criteria['name_or_description'],
                                          criteria['name_exact'])
    if 'kind' in criteria and criteria['kind']:
        query = query.filter(type__name=criteria['kind'])
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

    template = 'search_feature_template.html'
    mimetype = None  # default
    if request_data:
        form = SearchFeatureForm(request_data)
        # _update_form_init_values(form, database)
        if form.is_valid():
            search_criteria = form.cleaned_data
            search_criteria = dict([(key, value) for key, value in search_criteria.items() if value])
            context['search_criteria'] = search_criteria
            feature_queryset = _build_feature_query(search_criteria)
            download_search = search_criteria.get('download_search', False)
            if feature_queryset and download_search:
                context['queryset'] = feature_queryset
                template = 'search_feature_download_template.txt'
                mimetype = 'text/plain'
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

    context['form'] = form

    return render_to_response(template, context, mimetype=mimetype)
