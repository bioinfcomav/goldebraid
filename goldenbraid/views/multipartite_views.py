# -*- coding: utf-8 -*-
'''
Created on 2013 urt 17

@author: peio
'''
from django.core.exceptions import ValidationError
from django.forms.widgets import Select
from django.core.context_processors import csrf
from django.template.context import RequestContext
from django.shortcuts import render_to_response
from django.http import Http404, HttpResponseBadRequest, HttpResponse
from django import forms
from Bio.Seq import Seq

from goldenbraid.models import Feature
from goldenbraid.settings import DB

from goldenbraid.tags import VECTOR_TYPE_NAME
from goldenbraid.views.feature_views import get_prefix_and_suffix_index
from collections import OrderedDict


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
            used_parts = OrderedDict()
            for part_type in [p[0] for p in PARTS_TO_ASSEMBLE[multi_type]]:
                used_parts[part_type] = request_data[part_type]
            used_parts[VECTOR_TYPE_NAME] = request_data[VECTOR_TYPE_NAME]
            assembled_seq = _assemble_parts(multi_form_data, multi_type)
            return render_to_response('multipartite_result_template.html',
                                          {'assembled_seq': assembled_seq,
                                           'used_parts': used_parts,
                                           'multi_type': multi_type},
                                    context_instance=RequestContext(request))

    else:
        form = form_class()

    context['form'] = form

    template = 'multipartite_template.html'
    mimetype = None
    return render_to_response(template, context, mimetype=mimetype)


def write_protocol(protocol_data):
    "it writes the protocol in a variable"
    protocol = []
    protocol.append("Multipartite Assembly Protocol")
    protocol.append("")

    part_types = [p[0] for p in PARTS_TO_ASSEMBLE[protocol_data['multi_type']]]
    fragments = []
    for part_type in part_types:
        part_name = protocol_data[part_type]
        fragments.append(part_name)
    part_str = "({}){}".format(":".join(fragments), protocol_data[VECTOR_TYPE_NAME])

    protocol.append("Entities to assemble: {}".format(part_str))
    protocol.append("Reactions should be peformed as follows:")

    part_types.append(VECTOR_TYPE_NAME)
    for part_type in part_types:
        part_name = protocol_data[part_type]
        protocol.append("\t75 ng of {}".format(part_name))

    for enzyme in get_enzymes_for_protocol(protocol_data):
        protocol.append("\t3u of {}".format(enzyme))
    protocol.append("")
    protocol.append(u"\t1 microlitre Ligase Buffer")
    protocol.append("")
    protocol.append(u"Final volume: 10 microlitre")
    protocol.append("")
    long_line1 = "Set your reaction in a thermocycler: 25 cycles x "
    long_line1 += "(37C 2', 16C 5')."
    protocol.append(long_line1)

    lline2 = "One microlitre of the reaction is enough to be transform E.coli "
    lline2 += "electrocompetent cells. Positive clones are selected in {}"
    lline2 += " (50 micrograme ml-1), IPTG (0.5mM) and Xgal (40 micrograme ml-1) plates"
    lline2 += " You will distinguish between colonies carrying intact vectors "
    lline2 += "(blue) and those transformed with your construction (white)."
    vector = Feature.objects.using(DB).get(uniquename=protocol_data[VECTOR_TYPE_NAME])
    protocol.append(lline2.format(vector.resistance[0]))

    protocol = "\n".join(protocol)

    return protocol


def get_enzymes_for_protocol(protocol_data):
    'it gets the necesary enzymes'
    enzymes = set()
    vector = Feature.objects.using(DB).get(uniquename=protocol_data[VECTOR_TYPE_NAME])
    vec_enzyme_in = vector.enzyme_in[0]
    vec_enzyme_out = vector.enzyme_out

    enzymes.add(vec_enzyme_in)

    part_types = [p[0] for p in PARTS_TO_ASSEMBLE[protocol_data['multi_type']]]
    for part_type in part_types:
        part_name = protocol_data[part_type]
        part = Feature.objects.using(DB).get(uniquename=part_name)
        enzyme_outs = part.enzyme_out
        if vec_enzyme_in not in enzyme_outs:
            for enzyme_out in enzyme_outs:
                if enzyme_out != vec_enzyme_out:
                    enzymes.add(enzyme_out)
                    break
    return list(enzymes)


def multipartite_protocol_view(request):
    "it returns the protocol "
    if not request.GET:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    protocol = write_protocol(request.GET)
    return HttpResponse(protocol, mimetype='text/plain')


