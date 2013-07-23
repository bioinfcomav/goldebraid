# -*- coding: utf-8 -*-
'''
Created on 2013 urt 17

@author: peio
'''
import os
from tempfile import NamedTemporaryFile
from django.db.utils import IntegrityError
from django.http.response import HttpResponseServerError
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from django.forms.widgets import Select
from django.core.context_processors import csrf
from django.template.context import RequestContext
from django.shortcuts import render_to_response, redirect
from django.http import Http404, HttpResponseBadRequest, HttpResponse
from django import forms
from django.conf import settings as proj_settings
from django.db.models import Q

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO

from goldenbraid.models import Feature, Count
from goldenbraid.settings import (PARTS_TO_ASSEMBLE, UT_SUFFIX, UT_PREFIX,
                                  ASSEMBLED_SEQ)
from goldenbraid.tags import (VECTOR_TYPE_NAME, REVERSE, TU_TYPE_NAME,
                              MODULE_TYPE_NAME)
from goldenbraid.views.feature_views import get_prefix_and_suffix_index, \
    add_feature
from goldenbraid.forms import (get_multipartite_form,
                               get_multipartite_free_form,
                               MultipartiteFormFreeInitial,
                               features_to_choices, get_vector_choices)


def assemble_parts(parts, part_types):
    'We build the parts using the form data'
    part_types.append(VECTOR_TYPE_NAME)
    joined_seq = SeqRecord(Seq('', alphabet=generic_dna))
    names = {'parts': []}
    for part_type in part_types:
        part_uniquename = parts[part_type]
        part = Feature.objects.get(uniquename=part_uniquename)
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

    try:
        count = Count.objects.get(name=ASSEMBLED_SEQ)
    except Count.DoesNotExist:
        count = Count.objects.create(name=ASSEMBLED_SEQ, value=1)
    next_value = count.next

    joined_seq.id = ASSEMBLED_SEQ + '_' + next_value
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
    form_class = get_multipartite_form(multi_type, request.user)
    if request_data:
        form = form_class(request_data)

        if form.is_valid():
            multi_form_data = form.cleaned_data
            part_types = [p[0] for p in PARTS_TO_ASSEMBLE[multi_type]]
            assembled_seq = assemble_parts(multi_form_data, part_types)
            filename = assembled_seq.name + '.gb'
            response = HttpResponse(assembled_seq.format('genbank'),
                                    mimetype='text/plain')
            response['Content-Disposition'] = 'attachment; '
            response['Content-Disposition'] += 'filename="{0}"'.format(filename)
            return response
    return HttpResponseBadRequest()


def multipartite_view_add(request):
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST
    if not request_data:
        return render_to_response('goldenbraid_info.html',
                               {'info': "Not enougth data to add the feature"},
                                context_instance=RequestContext(request))
    multi_type = request_data['category']

    allowed_categories = PARTS_TO_ASSEMBLE.keys() + ['free']
    if multi_type is None or multi_type not in allowed_categories:
        return render_to_response('goldenbraid_info.html',
                               {'info': "Not enougth data to add the feature"},
                                context_instance=RequestContext(request))
    if multi_type == 'free':
        part_types = request_data['order'].split(':')
    else:
        part_types = [p[0] for p in PARTS_TO_ASSEMBLE[multi_type]]

    multi_data = {'Vector': request_data['Vector']}
    for part_type in part_types:
        multi_data[part_type] = request_data[part_type]

    assembled_seq = assemble_parts(multi_data, part_types)
    name = request_data['name']
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(assembled_seq.id),
                                    suffix='.gb')
    temp_fhand.write(assembled_seq.format('gb'))
    temp_fhand.flush()
    temp_fhand.seek(0)

    props = {'Description': [request_data['description']],
             'Reference': [request_data['reference']]}
    try:
        feature = add_feature(name=name, type_name=TU_TYPE_NAME,
                              vector=request_data['Vector'],
                              genbank=temp_fhand, props=props,
                              owner=request.user, is_public=False)

    except IntegrityError as error:
        print error
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render_to_response('feature_exists.html',
                                      {},
                            context_instance=RequestContext(request))
        else:
            return HttpResponseServerError()
    except Exception as error:
        print error
        return HttpResponseServerError()
    # if everithing os fine we show the just added feature
    return redirect(feature.url)


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
    form_class = get_multipartite_form(multi_type, request.user)
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
    protocol.append(u"\t3u of T4 ligase")
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
    vector = Feature.objects.get(uniquename=protocol_data[VECTOR_TYPE_NAME])
    protocol.append(lline2.format(vector.resistance[0]))

    protocol = "\n".join(protocol)

    return protocol


def get_enzymes_for_protocol(protocol_data, part_order):
    'it gets the necessary enzymes'

    vector = Feature.objects.get(uniquename=protocol_data[VECTOR_TYPE_NAME])
    vec_enzyme_in = vector.enzyme_in
    vec_enzyme_out = vector.enzyme_out
    enzymes = set(vec_enzyme_in)

    for part_type in part_order:
        part_name = protocol_data[part_type]
        part = Feature.objects.get(uniquename=part_name)
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
    response = HttpResponse(protocol, mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
    return response


def _get_fragments_from_request(request):
    post_data = request.POST
    vector = post_data['vector']
    parts = [post_data[k] for k in sorted(post_data.keys()) if 'part' in k]
    protocol_data = {VECTOR_TYPE_NAME: vector}
    part_order = []
    counter = {}
    for part in parts:
        feat = Feature.objects.get(uniquename=part)
        feat_type = feat.type.name
        if feat_type in protocol_data:
            if feat_type not in counter:
                counter[feat_type] = 1
            counter[feat_type] += 1
            feat_type += '.' + str(counter[feat_type])

        protocol_data[feat_type] = part
        part_order.append(feat_type)
    return protocol_data, part_order


def multipartite_view_free_protocol(request):
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    protocol_data, part_order = _get_fragments_from_request(request)
    protocol = write_protocol(protocol_data, 'multipartite', part_order)
    response = HttpResponse(protocol, mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
    return response


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

    form_class = get_multipartite_free_form(feats)
    form = form_class(request_data)
    if form.is_valid():
        last_feat = Feature.objects.get(uniquename=feats[-1])
        last_suffix = last_feat.suffix
        if last_suffix == 'CGCT':
            protocol_data, part_order = _get_fragments_from_request(request)
            assembled_seq = assemble_parts(protocol_data, part_order)
            response = HttpResponse(assembled_seq.format('genbank'),
                                    mimetype='text/plain')
            filename = assembled_seq.name + '.gb'
            response['Content-Disposition'] = 'attachment; '
            response['Content-Disposition'] += 'filename="{0}"'.format(filename)
            return response

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
    user = request.user
    if form_num is None:
        form = MultipartiteFormFreeInitial()
        form.fields['vector'].widget.choices = get_vector_choices(request.user)
        context['form_num'] = '1'
    else:
        if request_data:
            feats = [request_data['vector']]
            for k in sorted(request_data.keys()):
                if 'part_' in k:
                    feats.append(request_data[k])

            form_class = get_multipartite_free_form(feats)
            form = form_class(request_data)
            if form.is_valid():
                last_feat = Feature.objects.get(uniquename=feats[-1])
                last_suffix = last_feat.suffix
                if last_suffix == str(Seq(UT_SUFFIX).reverse_complement()):
                    last_suffix = UT_PREFIX
                if last_suffix == 'CGCT':
                    used_parts = OrderedDict({'Vector': feats[0]})
                    part_order = []
                    counters = {}
                    for feat in feats[1:]:
                        feat = Feature.objects.get(uniquename=feat)
                        feat_type = feat.type.name
                        if feat_type in used_parts:
                            if feat_type not in counters:
                                counters[feat_type] = 1
                            counters[feat_type] += 1
                            feat_type = '{0}.{1}'.format(feat_type, counters[feat_type])

                        used_parts[feat_type] = feat.uniquename
                        part_order.append(feat_type)

                    return render_to_response('multipartite_free_result.html',
                                              {'used_parts': used_parts,
                                               'multi_type': 'free',
                                               'post_data': form.cleaned_data,
                                               'order': ":".join(part_order)},
                                      context_instance=RequestContext(request))
                else:
                    # add new_field
                    part_num = len(feats)
                    feats = Feature.objects.filter(prefix=last_suffix)
                    if user.is_staff:
                        pass
                    elif user.is_authenticated():
                        feats = feats.filter(Q(featureperm__owner__username=user) |
                                             Q(featureperm__is_public=True))

                    else:
                        feats = feats.filter(featureperm__is_public=True)

                    feats = feats.exclude(type__name__in=[VECTOR_TYPE_NAME,
                                                          TU_TYPE_NAME,
                                                          MODULE_TYPE_NAME])
                    choices = features_to_choices(feats)
                    form.fields['part_{0}'.format(part_num)] = forms.CharField(max_length=100,
                                                widget=Select(choices=choices),
                                                required=True)
                    context['form_num'] = str(part_num)

    if form is None:
        form_class = get_multipartite_free_form()
        form = form_class()
        context['form_num'] = '1'

    context['form'] = form
    template = 'multipartite_free_template.html'
    mimetype = None
    return render_to_response(template, context, mimetype=mimetype)
