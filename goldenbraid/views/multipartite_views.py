'''
Created on 2013 urt 17

@author: peio
'''
from django.core.exceptions import ValidationError
from django.forms.widgets import Select
from django.core.context_processors import csrf
from django.template.context import RequestContext
from django.shortcuts import render_to_response
from django.http import Http404
from django import forms
from Bio.Seq import Seq

from goldenbraid.models import Feature
from goldenbraid.settings import DB

from goldenbraid.tags import VECTOR_TYPE_NAME
from goldenbraid.views.feature_views import get_prefix_and_suffix_index


PARTS_TO_ASSEMBLE = {'basic': [('PROM+UTR+ATG', 'GGAG', 'AATG'),
                               ('CDS', 'AATG', 'GCTT'),
                               ('TER', 'GCTT', 'CGCT')],
                     'secreted': [('PROM+UTR+ATG', 'GGAG', 'AATG'),
                                 ('SP', 'AATG', 'AGCC'),
                                 ('CDS', 'AGCC', 'GCTT'),
                                 ('TER', 'GCTT', 'CGCT')]
                     }


def create_field_validator(field_name):

    def validator(self):
        uniquename_str = self.cleaned_data[field_name]
        try:
            Feature.objects.using(DB).get(uniquename=uniquename_str)
        except Feature.DoesNotExist:
            raise ValidationError('This feature does not exist in the database')
        return uniquename_str

    return validator


def _get_multipartite_form(multi_type):
    'It returns a form for the given multipartite'
    form_fields = {}

    part_defs = PARTS_TO_ASSEMBLE[multi_type]

    # first we need to add the vector to the form
    vector_choices = [('', '')]
    vector_suffix = part_defs[0][1]
    vector_prefix = part_defs[-1][2]
    vectors = Feature.objects.using(DB).filter(type__name=VECTOR_TYPE_NAME)
    vectors = vectors.filter(prefix=vector_prefix, suffix=vector_suffix)
    for vector in vectors:
        vector_choices.append((vector.uniquename, vector.uniquename))

    form_fields[VECTOR_TYPE_NAME] = forms.CharField(max_length=100,
                                        widget=Select(choices=vector_choices))

    for parts in part_defs:
        choices = [('', '')]
        for feat in Feature.objects.using(DB).filter(type__name=parts[0],
                                                     prefix=parts[1],
                                                     suffix=parts[2]):
            choices.append((feat.uniquename, feat.uniquename))

        name = parts[0]
        form_fields[name] = forms.CharField(max_length=100,
                                            widget=Select(choices=choices))

    form = type('MultiPartiteForm', (forms.BaseForm,),
                {'base_fields': form_fields})
    for field_name in form_fields.keys():
        setattr(form, 'clean_{}'.format(field_name),
                create_field_validator(field_name))
    return form


def _assemble_parts(parts, multi_type):
    'We build the parts using the form data'
    part_types = [p[0] for p in PARTS_TO_ASSEMBLE[multi_type]]
    part_types.append(VECTOR_TYPE_NAME)
    joined_seq = Seq('')
    for part_type in part_types:
        part_uniquename = parts[part_type]
        part = Feature.objects.using(DB).get(uniquename=part_uniquename)
        seq = Seq(part.residues)
        if part.type.name == VECTOR_TYPE_NAME:
            enzyme = part.enzyme_in[0]
        else:
            enzyme = part.enzyme_out[0]
        pref_idx, suf_idx = get_prefix_and_suffix_index(seq, enzyme)[:2]
        if suf_idx >= pref_idx:
            part_sub_seq = seq[pref_idx:suf_idx]
        else:
            part_sub_seq = seq[pref_idx:]
            part_sub_seq += seq[:suf_idx]
        joined_seq += part_sub_seq
    return joined_seq


def multipartite_view(request, multi_type=None):
    'view of the multipartite tool'
    if multi_type is None:
        return render_to_response('multipartite_initial.html', {},
                                  context_instance=RequestContext(request))
    elif multi_type not in PARTS_TO_ASSEMBLE.keys():
        return Http404

    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    form_class = _get_multipartite_form(multi_type)
    if request_data:
        form = form_class(request_data)
        if form.is_valid():
            multi_form_data = form.cleaned_data
            assembled_seq = _assemble_parts(multi_form_data, multi_type)
            return render_to_response('multipartite_result_template.html',
                                          {'assembled_seq': assembled_seq,
                                           'used_parts': request_data},
                                    context_instance=RequestContext(request))

    else:
        form = form_class()

    context['form'] = form

    template = 'multipartite_template.html'
    mimetype = None
    return render_to_response(template, context, mimetype=mimetype)

