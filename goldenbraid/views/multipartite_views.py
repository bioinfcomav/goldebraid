# -*- coding: utf-8 -*-
'''
Created on 2013 urt 17

@author: peio
'''
import os
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from django.core.exceptions import ValidationError
from django.forms.widgets import Select
from django.core.context_processors import csrf
from django.template.context import RequestContext
from django.shortcuts import render_to_response
from django.http import Http404, HttpResponseBadRequest, HttpResponse
from django import forms
from django.conf import settings as proj_settings

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO

from goldenbraid.models import Feature
from goldenbraid.settings import DB
from goldenbraid.tags import VECTOR_TYPE_NAME, FORWARD, REVERSE
from goldenbraid.views.feature_views import get_prefix_and_suffix_index


PARTS_TO_ASSEMBLE = {'basic': [('PROM+UTR+ATG', 'GGAG', 'AATG'),
                               ('CDS', 'AATG', 'GCTT'),
                               ('TER', 'GCTT', 'CGCT')],
                     'secreted': [('PROM+UTR+ATG', 'GGAG', 'AATG'),
                                 ('SP', 'AATG', 'AGCC'),
                                 ('CDS', 'AGCC', 'GCTT'),
                                 ('TER', 'GCTT', 'CGCT')],
                     'ct-fusion': [('PROM+UTR+ATG', 'GGAG', 'AATG'),
                                   ('CDS', 'AATG', 'GCAG'),
                                   ('CT', 'GCAG', 'GCTT'),
                                   ('TER', 'GCTT', 'CGCT')],
                     'nt-fusion': [('PROM+UTR', 'GGAG', 'CCAT'),
                                   ('NT', 'CCAT', 'AATG'),
                                   ('CDS', 'AATG', 'GCTT'),
                                   ('TER', 'GCTT', 'CGCT')],
                     'nt-ct-fusion': [('PROM+UTR', 'GGAG', 'CCAT'),
                                      ('NT', 'CCAT', 'AATG'),
                                      ('CDS', 'AATG', 'GCAG'),
                                      ('CT', 'GCAG', 'GCTT'),
                                      ('TER', 'GCTT', 'CGCT')],
                     'operated-promoter-a': [('OP', 'GGAG', 'TCCC'),
                                           ('MinPROM', 'TCCC', 'AATG'),
                                           ('CDS', 'AATG', 'GCTT'),
                                           ('TER', 'GCTT', 'CGCT')],
                     'operated-promoter-b': [('PROM', 'GGAG', 'TGAC'),
                                           ('OP', 'TGAC', 'TCCC'),
                                           ('MinPROM', 'TCCC', 'AATG'),
                                           ('CDS', 'AATG', 'GCTT'),
                                           ('TER', 'GCTT', 'CGCT')],
                     'protein-interaction': [('InteractionADAPTOR', 'GGAG',
                                              'AATG'),
                                             ('CDS', 'AATG', 'GCTT'),
                                             ('TER', 'GCTT', 'CGCT')],
                     'amiRNA':  [('PROM+UTR', 'GGAG', 'CCAT'),
                                 ('5FS', 'CCAT', 'GTGA'),
                                 ('Target', 'GTGA', 'TCTC'),
                                 ('3FS', 'TCTC', 'GCTT'),
                                 ('TER', 'GCTT', 'CGCT')],
                     'hpRNA':  [('PROM+UTR', 'GGAG', 'CCAT'),
                                ('goi', 'CCAT', 'GTGA'),
                                ('int', 'GTGA', 'TCTC'),
                                ('iog', 'TCTC', 'GCTT'),
                                ('TER', 'GCTT', 'CGCT')],
                     'tasiRNA':  [('PROM+UTR+mir173', 'GGAG', 'CCAT'),
                                  ('goi', 'CCAT', 'GCTT'),
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

    for_vectors = vectors.filter(prefix=vector_prefix, suffix=vector_suffix)
    rev_vectors = vectors.filter(prefix=Seq(vector_suffix).reverse_complement(),
                                     suffix=Seq(vector_prefix).reverse_complement())

    vector_choices.append(('', 'Forward vectors'))
    for vector in for_vectors:
        vector_choices.append((vector.uniquename, vector.uniquename))
    vector_choices.append(('', 'Reverse vectors'))
    for vector in rev_vectors:
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
        setattr(form, 'clean_{0}'.format(field_name),
                create_field_validator(field_name))
    return form


def _assemble_parts(parts, multi_type):
    'We build the parts using the form data'
    part_types = [p[0] for p in PARTS_TO_ASSEMBLE[multi_type]]
    part_types.append(VECTOR_TYPE_NAME)
    joined_seq = SeqRecord(Seq('', alphabet=generic_dna))
    for part_type in part_types:
        part_uniquename = parts[part_type]
        part = Feature.objects.using(DB).get(uniquename=part_uniquename)
        gb_path = os.path.join(proj_settings.MEDIA_ROOT,
                               part.genbank_file.name)
        part_record = SeqIO.read(gb_path, 'gb')
        seq = Seq(part.residues)
        if part.type.name == VECTOR_TYPE_NAME:
            enzyme = part.enzyme_in[0]
        else:
            enzyme = part.enzyme_out[0]
        pref_idx, suf_idx = get_prefix_and_suffix_index(seq, enzyme)[:2]
        if suf_idx >= pref_idx:
            part_sub_seq = part_record[pref_idx:suf_idx]
        else:
            part_sub_seq = part_record[pref_idx:]
            part_sub_seq += part_record[:suf_idx]

        # VECTOR must be always the last part to add
        if part.type.name == VECTOR_TYPE_NAME and part.direction == REVERSE:
            joined_seq = joined_seq.reverse_complement() + part_sub_seq
        else:
            joined_seq += part_sub_seq

    joined_seq.id = 'assembled_parts'
    joined_seq.name = joined_seq.id

    return joined_seq


def multipartite_view_genbank(request, multi_type=None):
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

            return  HttpResponse(assembled_seq.format('genbank'),
                                 mimetype='text/plain')
    return HttpResponseBadRequest()


def multipartite_view(request, multi_type=None):
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
            used_parts = OrderedDict()
            multi_form_data = form.cleaned_data
            for part_type in [p[0] for p in PARTS_TO_ASSEMBLE[multi_type]]:
                used_parts[part_type] = multi_form_data[part_type]
                used_parts[VECTOR_TYPE_NAME] = multi_form_data[VECTOR_TYPE_NAME]
            posted_data = multi_form_data
            return render_to_response('multipartite_result_template.html',
                                      {'used_parts': used_parts,
                                       'multi_type': multi_type,
                                       'posted_data': posted_data},
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
    part_str = "({0}){1}".format(":".join(fragments), protocol_data[VECTOR_TYPE_NAME])

    protocol.append("Entities to assemble: {0}".format(part_str))
    protocol.append("Reaction should be performed as follows:")

    part_types.append(VECTOR_TYPE_NAME)
    for part_type in part_types:
        part_name = protocol_data[part_type]
        protocol.append("\t75 ng of {0}".format(part_name))

    for enzyme in get_enzymes_for_protocol(protocol_data):
        protocol.append("\t3u of {0}".format(enzyme))
    protocol.append("")
    protocol.append(u"\t1 microlitre Ligase Buffer")
    protocol.append("")
    protocol.append(u"Final volume: 10 microlitre")
    protocol.append("")
    long_line1 = "Set your reaction in a thermocycler: 25 cycles x "
    long_line1 += "(37C 2', 16C 5')."
    protocol.append(long_line1)

    lline2 = "One microlitre of the reaction is enough to be transform E.coli "
    lline2 += "electrocompetent cells. Positive clones are selected in {0}"
    lline2 += " (50 microgram ml-1), IPTG (0.5mM) and Xgal (40 microgram ml-1) plates"
    lline2 += " You will distinguish between colonies carrying intact vectors "
    lline2 += "(blue) and those transformed with your construction (white)."
    vector = Feature.objects.using(DB).get(uniquename=protocol_data[VECTOR_TYPE_NAME])
    protocol.append(lline2.format(vector.resistance[0]))

    protocol = "\n".join(protocol)

    return protocol


def get_enzymes_for_protocol(protocol_data):
    'it gets the necessary enzymes'
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
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    protocol = write_protocol(request.POST)
    return HttpResponse(protocol, mimetype='text/plain')


