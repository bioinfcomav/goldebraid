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
from goldenbraid.tags import (VECTOR_TYPE_NAME, REVERSE, TU_TYPE_NAME,
                              PHRASE_TYPE_NAME, MODULE_TYPE_NAME)
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
UT_PREFIX = PARTS_TO_ASSEMBLE['basic'][0][1]
UT_SUFFIX = PARTS_TO_ASSEMBLE['basic'][-1][2]


def create_feature_validator(field_name):

    def validator(self):
        uniquename_str = self.cleaned_data[field_name]
        try:
            Feature.objects.using(DB).get(uniquename=uniquename_str)
        except Feature.DoesNotExist:
            raise ValidationError('This feature does not exist in the database')
        return uniquename_str

    return validator


def vectors_to_choice(vectors):
    "it returns the given vectors but prepared to use as choices in a select"
    for_vectors = vectors.filter(prefix=UT_SUFFIX, suffix=UT_PREFIX)
    rev_vectors = vectors.filter(prefix=Seq(UT_PREFIX).reverse_complement(),
                                 suffix=Seq(UT_SUFFIX).reverse_complement())
    for_vector_choices = []
    for vector in for_vectors:
        for_vector_choices.append((vector.uniquename, vector.uniquename))
    rev_vector_choices = []
    for vector in rev_vectors:
        rev_vector_choices.append((vector.uniquename, vector.uniquename))

    vector_choices = (('', ''),
                      ('Forward vectors', for_vector_choices),
                      ('Reverse vectors', rev_vector_choices))

    return vector_choices


def _get_multipartite_form(multi_type):
    'It returns a form for the given multipartite'
    form_fields = OrderedDict()

    part_defs = PARTS_TO_ASSEMBLE[multi_type]
    for parts in part_defs:
        choices = [('', '')]
        for feat in Feature.objects.using(DB).filter(type__name=parts[0],
                                                     prefix=parts[1],
                                                     suffix=parts[2]):
            choices.append((feat.uniquename, feat.uniquename))

        name = parts[0]
        form_fields[name] = forms.CharField(max_length=100,
                                            widget=Select(choices=choices))

    # last we need to add the vector to the form
    vectors = Feature.objects.using(DB).filter(type__name=VECTOR_TYPE_NAME)
    vector_choices = vectors_to_choice(vectors)
    form_fields[VECTOR_TYPE_NAME] = forms.CharField(max_length=100,
                                        widget=Select(choices=vector_choices))

    form = type('MultiPartiteForm', (forms.BaseForm,),
                {'base_fields': form_fields})
    for field_name in form_fields.keys():
        setattr(form, 'clean_{0}'.format(field_name),
                create_feature_validator(field_name))
    return form


def assemble_parts(parts, part_types):
    'We build the parts using the form data'
    part_types.append(VECTOR_TYPE_NAME)
    joined_seq = SeqRecord(Seq('', alphabet=generic_dna))
    names = {'parts': []}
    for part_type in part_types:
        part_uniquename = parts[part_type]
        part = Feature.objects.using(DB).get(uniquename=part_uniquename)
        gb_path = os.path.join(proj_settings.MEDIA_ROOT,
                               part.genbank_file.name)
        part_record = SeqIO.read(gb_path, 'gb')
        seq = Seq(part.residues)
        if part.type.name == VECTOR_TYPE_NAME:
            enzyme = part.enzyme_in[0]
            names['vector'] = part.uniquename
        else:
            names['parts'].append(part.uniquename)
            enzyme = part.enzyme_out[0]

        pref_idx, suf_idx = get_prefix_and_suffix_index(seq, enzyme)[:2]
        # VECTOR must be always the last part to add
        if part.type.name == VECTOR_TYPE_NAME and part.direction == REVERSE:
            if suf_idx >= pref_idx and suf_idx + 4 < len(seq):
                part_sub_seq = part_record[pref_idx + 4:suf_idx + 4]
            else:
                part_sub_seq = part_record[pref_idx + 4:]
                part_sub_seq += part_record[:suf_idx + 4]

            joined_seq = joined_seq.reverse_complement() + part_sub_seq
        else:
            if suf_idx >= pref_idx:
                part_sub_seq = part_record[pref_idx:suf_idx]
            else:
                part_sub_seq = part_record[pref_idx:]
                part_sub_seq += part_record[:suf_idx]
            joined_seq += part_sub_seq

    joined_seq.id = 'assembled_seq'
    joined_seq.name = joined_seq.id
    joined_seq.description = "({0}){1}".format(','.join(names['parts']),
                                             names['vector'])

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
            part_types = [p[0] for p in PARTS_TO_ASSEMBLE[multi_type]]
            assembled_seq = assemble_parts(multi_form_data, part_types)

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
            return render_to_response('multipartite_result.html',
                                      {'used_parts': used_parts,
                                       'multi_type': multi_type},
                                context_instance=RequestContext(request))
    else:
        form = form_class()

    context['form'] = form

    template = 'multipartite_template.html'
    mimetype = None
    return render_to_response(template, context, mimetype=mimetype)


def write_protocol(protocol_data, assembly_type, part_order):
    "it writes the protocol in a variable"
    protocol = []
    protocol.append("{0} Assembly Protocol".format(assembly_type.title()))
    protocol.append("")

    fragments = []
    for part_type in part_order:
        part_name = protocol_data[part_type]
        fragments.append(part_name)
    part_str = "({0}){1}".format(":".join(fragments), protocol_data[VECTOR_TYPE_NAME])

    protocol.append("Entities to assemble: {0}".format(part_str))
    protocol.append("Reaction should be performed as follows:")

    part_order.append(VECTOR_TYPE_NAME)
    for part_type in part_order:
        part_name = protocol_data[part_type]
        protocol.append("\t75 ng of {0}".format(part_name))

    for enzyme in get_enzymes_for_protocol(protocol_data, part_order):
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


def get_enzymes_for_protocol(protocol_data, part_order):
    'it gets the necessary enzymes'

    vector = Feature.objects.using(DB).get(uniquename=protocol_data[VECTOR_TYPE_NAME])
    vec_enzyme_in = vector.enzyme_in
    vec_enzyme_out = vector.enzyme_out
    enzymes = set(vec_enzyme_in)

    for part_type in part_order:
        part_name = protocol_data[part_type]
        part = Feature.objects.using(DB).get(uniquename=part_name)
        enzyme_outs = part.enzyme_out
        if vec_enzyme_in not in enzyme_outs:
            for enzyme_out in enzyme_outs:
                if enzyme_out not in  vec_enzyme_out:
                    enzymes.add(enzyme_out)
                    break
    return list(enzymes)


def multipartite_protocol_view(request):
    "it returns the protocol "
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    part_order = [p[0] for p in PARTS_TO_ASSEMBLE[request.POST['multi_type']]]
    protocol = write_protocol(request.POST, 'multipartite', part_order)
    return HttpResponse(protocol, mimetype='text/plain')


class MultipartiteFormFreeInitial(forms.Form):
    vectors = Feature.objects.using(DB).filter(type__name=VECTOR_TYPE_NAME)
    choices = vectors_to_choice(vectors)
    vector = forms.CharField(max_length=100, widget=Select(choices=choices))

    def clean_vector(self):
        return create_feature_validator('vector')(self)


def _get_multipartite_free_form(feat_uniquenames):
    form_fields = OrderedDict()
    count = 0
    for feat_uniquename in feat_uniquenames:
        if count == 0:
            part_name = 'vector'
        else:
            part_name = 'part_{0}'.format(count)
        field = forms.CharField(required=True, initial=feat_uniquename)
        form_fields[part_name] = field
        count += 1

    form = type('MultiPartiteFreeValForm', (forms.BaseForm,),
                {'base_fields': form_fields})

    for field_name in form_fields.keys():
        setattr(form, 'clean_{0}'.format(field_name),
                create_feature_validator(field_name))
    return form


def _get_fragments_from_request(request):
    post_data = request.POST
    vector = post_data['vector']
    parts = [post_data[k] for k in sorted(post_data.keys()) if 'part' in k]
    protocol_data = {VECTOR_TYPE_NAME: vector}
    part_order = []
    for part in parts:
        feat = Feature.objects.using(DB).get(uniquename=part)
        protocol_data[feat.type.name] = part
        part_order.append(feat.type.name)
    return protocol_data, part_order


def multipartite_view_free_protocol(request):
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    protocol_data, part_order = _get_fragments_from_request(request)
    protocol = write_protocol(protocol_data, 'multipartite', part_order)
    return HttpResponse(protocol, mimetype='text/plain')


def multipartite_view_free_genbank(request):
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    else:
        msg = "To get the sequence you need first to assemble parts"
        return HttpResponseBadRequest(msg)

    feats = [request_data['vector']]
    for k in sorted(request_data.keys()):
        if 'part_' in k:
            feats.append(request_data[k])

    form_class = _get_multipartite_free_form(feats)
    form = form_class(request_data)
    if form.is_valid():
        last_feat = Feature.objects.using(DB).get(uniquename=feats[-1])
        last_suffix = last_feat.suffix
        if last_suffix == 'CGCT':
            protocol_data, part_order = _get_fragments_from_request(request)
            assembled_seq = assemble_parts(protocol_data, part_order)
            return  HttpResponse(assembled_seq.format('genbank'),
                                 mimetype='text/plain')

    return HttpResponseBadRequest('There was an error in the assembly')


def multipartite_view_free(request, form_num):
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
        form = MultipartiteFormFreeInitial()
        context['form_num'] = '1'
    else:
        if request_data:
            feats = [request_data['vector']]
            for k in sorted(request_data.keys()):
                if 'part_' in k:
                    feats.append(request_data[k])

            form_class = _get_multipartite_free_form(feats)
            form = form_class(request_data)
            if form.is_valid():
                last_feat = Feature.objects.using(DB).get(uniquename=feats[-1])
                last_suffix = last_feat.suffix
                if last_suffix == 'CGCT':
                    used_parts = OrderedDict({'vector': feats[0]})
                    for feat in feats[1:]:
                        feat = Feature.objects.using(DB).get(uniquename=feat)
                        used_parts[feat.type.name] = feat.uniquename
                    return render_to_response('multipartite_free_result.html',
                                              {'used_parts': used_parts,
                                               'multi_type': 'free',
                                               'post_data': form.cleaned_data},
                                      context_instance=RequestContext(request))
                    return  HttpResponse('result', mimetype='text/plain')
                else:
                    # add new_field
                    part_num = len(feats)
                    choices = [('', '')]
                    feats = Feature.objects.using(DB).filter(prefix=last_suffix)
                    feats = feats.exclude(type__name__in=[VECTOR_TYPE_NAME,
                                                          TU_TYPE_NAME,
                                                          MODULE_TYPE_NAME,
                                                          PHRASE_TYPE_NAME])
                    for feat in feats:
                        choices.append((feat.uniquename, feat.uniquename))
                    form.fields['part_{0}'.format(part_num)] = forms.CharField(max_length=100,
                                                widget=Select(choices=choices),
                                                required=True)
                    context['form_num'] = str(part_num)

    if form is None:
        form_class = _get_multipartite_free_form()
        form = form_class()
        context['form_num'] = '1'

    context['form'] = form
    template = 'multipartite_free_template.html'
    mimetype = None
    return render_to_response(template, context, mimetype=mimetype)
