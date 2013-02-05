'''
Created on 2013 ots 5

@author: peio
'''
from Bio.Seq import Seq
from django.template.context import RequestContext
from django.core.context_processors import csrf
from django.shortcuts import render_to_response
from django import forms
from django.forms.widgets import Select
from django.db.models import Q

from goldenbraid.models import Feature
from goldenbraid.settings import DB
from goldenbraid.views.multipartite_views import (create_feature_validator,
                                                  vector_by_direction_choice,
                                                  assemble_parts,
                                                  PARTS_TO_ASSEMBLE)
from goldenbraid.tags import VECTOR_TYPE_NAME, ENZYME_IN_TYPE_NAME
from django.http import HttpResponse, HttpResponseBadRequest


BIPARTITE_ALLOWED_PARTS = ('TU', 'Phrase')

SITE_A = 'GGAG'
SITE_B = 'CGCT'
SITE_C = 'GTCA'


def _select_part_choices(parts):
    parts_forw = parts.filter(vector__prefix=SITE_B, vector__suffix=SITE_A)
    parts_rev = parts.filter(vector__prefix=Seq(SITE_A).reverse_complement(),
                             vector__suffix=Seq(SITE_B).reverse_complement())
    part_forw_choices = []
    for part_forw in parts_forw:
        part_forw_choices.append((part_forw.uniquename,
                                   part_forw.uniquename))
    part_rev_choices = []
    for part_rev in parts_rev:
        part_rev_choices.append((part_rev.uniquename,
                                  part_rev.uniquename))
    part_choices = (('', ''),
                      ('Forward parts', part_forw_choices),
                      ('Reverse parts', part_rev_choices))
    return part_choices


class BipartiteForm1(forms.Form):
    # vector
    _bi_parts = Feature.objects.using(DB).filter(
                                        type__name__in=BIPARTITE_ALLOWED_PARTS)
    _parts = _bi_parts.filter(prefix=SITE_A, suffix=SITE_C)
    _part_choices = _select_part_choices(_parts)

    part_1 = forms.CharField(max_length=100,
                                         widget=Select(choices=_part_choices))

    def clean_part_1(self):
        return create_feature_validator('part_1')(self)


class BipartiteForm2(forms.Form):
    part_1 = forms.CharField(max_length=100)
    part_1.widget.attrs['readonly'] = True

    part_2 = forms.CharField(max_length=100, widget=Select(choices=[]))

    def clean_part_2(self):
        return create_feature_validator('part_2')(self)


class BipartiteForm3(forms.Form):
    part_1 = forms.CharField(max_length=100)
    part_1.widget.attrs['readonly'] = True

    part_2 = forms.CharField(max_length=100)
    part_2.widget.attrs['readonly'] = True

    Vector = forms.CharField(max_length=100, widget=Select(choices=[]))

    def clean_part_1(self):
        return create_feature_validator('part_1')(self)

    def clean_part_2(self):
        return create_feature_validator('part_2')(self)

    def clean_Vector(self):
        return create_feature_validator('Vector')(self)


def _get_part2_choices(part1_uniquename):
    part1 = Feature.objects.using(DB).get(uniquename=part1_uniquename)
    part1_enzyme_out = part1.enzyme_out
    bi_parts = Feature.objects.using(DB).filter(
                                        type__name__in=BIPARTITE_ALLOWED_PARTS)
    parts = bi_parts.filter(prefix=SITE_C, suffix=SITE_B)

    parts_forw = parts.filter(vector__prefix=SITE_B, vector__suffix=SITE_A)
    parts_rev = parts.filter(vector__prefix=Seq(SITE_A).reverse_complement(),
                             vector__suffix=Seq(SITE_B).reverse_complement())
    part_forw_choices = []
    for part_forw in parts_forw:
        if part_forw.enzyme_out == part1_enzyme_out:
            part_forw_choices.append((part_forw.uniquename,
                                      part_forw.uniquename))
    part_rev_choices = []
    for part_rev in parts_rev:
        if part_rev.enzyme_out == part1_enzyme_out:
            part_rev_choices.append((part_rev.uniquename,
                                     part_rev.uniquename))
    part_choices = (('', ''),
                      ('Forward parts', part_forw_choices),
                      ('Reverse parts', part_rev_choices))
    return part_choices


def _get_bipart_vector_choices(part_uniquename):
    part = Feature.objects.using(DB).get(uniquename=part_uniquename)
    part_enzyme_out = part.enzyme_out

    vectors = Feature.objects.using(DB).filter(type__name=VECTOR_TYPE_NAME)
    vectors = vectors.filter(featureprop__type__name=ENZYME_IN_TYPE_NAME,
                             featureprop__value=part_enzyme_out)
    print vectors.query
    for vector in vectors:
        print vector.enzyme_in, part_enzyme_out
    part_def = PARTS_TO_ASSEMBLE.values()[0]

    vector_suffix = part_def[0][1]
    vector_prefix = part_def[-1][2]
    return vector_by_direction_choice(vectors, vector_prefix, vector_suffix)


def bipartite_view(request, form_num):
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    form = None
    if form_num is None:
        form = BipartiteForm1()
        context['form_num'] = '1'
    elif form_num == '1':
        if request_data:
            form = BipartiteForm1(request_data)
            if form.is_valid():
                form1_data = form.cleaned_data
                form = BipartiteForm2()
                # form.fields['Vector'].initial = form1_data['Vector']
                form.fields['part_1'].initial = form1_data['part_1']
                choices_part2 = _get_part2_choices(form1_data['part_1'])
                form.fields['part_2'].widget.choices = choices_part2
                context['form_num'] = '2'
            else:
                context['form_num'] = '1'
    elif form_num == '2':
        if request_data:
            form = BipartiteForm2(request_data)
            if form.is_valid():
                form2_data = form.cleaned_data
                form = BipartiteForm3()
                form.fields['part_1'].initial = form2_data['part_1']
                form.fields['part_2'].initial = form2_data['part_2']
                choices_vector = _get_bipart_vector_choices(form2_data['part_1'])
                form.fields['Vector'].widget.choices = choices_vector
                context['form_num'] = '3'
    elif form_num == '3':
        if request_data:
            form = BipartiteForm3(request_data)
            if form.is_valid():
                return render_to_response('bipartite_result.html',
                                          {'used_parts': form.cleaned_data},
                                    context_instance=RequestContext(request))
            else:
                context['form_num'] = '3'
    if form is None:
        form = BipartiteForm1()
        context['form_num'] = '1'
    context['form'] = form
    template = 'bipartite_template.html'
    mimetype = None
    return render_to_response(template, context, mimetype=mimetype)


def bipartite_view_genbank(request):
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        form = BipartiteForm3(request_data)
        if form.is_valid():
            seq = assemble_parts(form.cleaned_data, ['part_1', 'part_2'])
            return  HttpResponse(seq.format('genbank'),
                             mimetype='text/plain')
    return HttpResponseBadRequest()
