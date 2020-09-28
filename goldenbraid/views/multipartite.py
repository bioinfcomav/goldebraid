# -*- coding: utf-8 -*-

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

import textwrap
import json
import os
from os.path import join
import re
from tempfile import NamedTemporaryFile
from django.db.utils import IntegrityError
from django.http.response import HttpResponseServerError
from goldenbraid.sbol import convert_to_sbol
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from django.forms.widgets import Select
from django.template.context_processors import csrf
from django.template.context import RequestContext
from django.shortcuts import render, redirect
from django.http import Http404, HttpResponseBadRequest, HttpResponse
from django import forms
from django.conf import settings as proj_settings
from django.db.models import Q


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import (SeqIO, SeqFeature)
from Bio.SeqFeature import FeatureLocation, SeqFeature, ExactPosition

from goldenbraid.models import Feature, Count
from goldenbraid.settings import (PARTS_TO_ASSEMBLE, UT_SUFFIX, UT_PREFIX,
                                  ASSEMBLED_SEQ, REBASE_FILE,
                                  CRYSPER_TARGETS_TO_DOMESTICATE,
                                  CRYSPR_MULTIPLEX_EDITING_LEVEL_MINUS_ONE,
                                  CRYSPR_MULTIPLEX_REGULATION_LEVEL_MINUS_ONE
                                  )
from goldenbraid.tags import (VECTOR_TYPE_NAME, REVERSE, TU_TYPE_NAME,
                              MODULE_TYPE_NAME, TRNA,
                              SCAFFOLD, CRISPR_MULTIPLEXING_TARGET,
                              CRISPR_EDITING, CRISPR_REGULATION, GRNA_CAS12_DICOT,
                              FUNGAL_TU_TYPE_NAME, REVERSE_MARKER, FLANK_5UTR,
                              FORWARD_MARKER, FLANK_3UTR, 
                              FUNGAL_PROM_5UTR, FUNGAL_CDS,
                              FUNGAL_UTR3_TERM,
                              FUNGAL_KNOCK_OUT,
                              TARGET_CAS12A, PROM_CAS12)
from goldenbraid.utils import get_prefix_and_suffix_index, parse_rebase_file, has_rec_sites
from goldenbraid.views.feature import add_feature
from goldenbraid.forms.assemblers import (get_multipartite_form,
                                          get_multipartite_free_form,
                                          MultipartiteFormFreeInitial,
                                          features_to_choices,
                                          get_vector_choices,
                                          Cas12aSingleTUForm,
                                          Cas12aMultiplexTUForm,
                                          populate_cas12a_single_tu_form,
                                          populate_cas12a_multiplex_tu_form,
                                          FungalTUForm,
                                          populate_fungal_tu_form,
                                          GeneDisruptionForm,
                                          populate_gene_disruption_form)

TMP_DIR = join("/home/golden/devel/gbdb", "files")


def fungal_gene_disruption_view(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == "GET":
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        form = GeneDisruptionForm(request_data)
        if form.is_valid():
            cleaned_data = form.cleaned_data
            part_types = [REVERSE_MARKER,
                          FLANK_5UTR,
                          FORWARD_MARKER,
                          FLANK_3UTR]
            multi_form_data = OrderedDict([(REVERSE_MARKER, cleaned_data["negative_marker"]),
                                           (FLANK_5UTR, cleaned_data["five_prime_flank"]),
                                           (FORWARD_MARKER, cleaned_data["positive_marker"]),
                                           (FLANK_3UTR, cleaned_data["three_prime_flank"]),
                                           ("Vector", cleaned_data["vector"])])
            context.update({'used_parts': multi_form_data, 'multi_type': FUNGAL_KNOCK_OUT,
                            'part_types': part_types, "Vector": cleaned_data["vector"]})
            return render(request, 'fungal_gene_disruption_results.html', context=context)


    form = GeneDisruptionForm()
    populate_gene_disruption_form(form, request.user)
    context.update({'form': form})
    content_type = None
    return render(request, "fungal_gene_disruption.html", context=context)


def crispr_view_cas12a_multiplexing_TU(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        form = Cas12aMultiplexTUForm(request_data)
        if form.is_valid():
            part_types = ["promoter", "target"]
            cleaned_data = form.cleaned_data
            multi_form_data = OrderedDict([("part_1", cleaned_data["promoter"]),
                                           ("part_2", cleaned_data["target"]),
                                           ("vector", cleaned_data["vector"])])
            order = "part_1:part_2:Vector"
            context.update({'used_parts': multi_form_data, 'multi_type': "free",
                            'part_types': part_types, "order": order, "Vector": cleaned_data["vector"]})
            return render(request, 'cas12_multiplex_tu_results.html', context=context)


    form = Cas12aMultiplexTUForm()
    populate_cas12a_multiplex_tu_form(form, request.user)

    context.update({'form': form})
    template = 'cas12a_multiplex_tu_template.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def _get_position_feature(position):
    position_features = Feature.objects.filter(type__name=position).filter(featureperm__owner__username="admin")
    position_features = position_features.filter(prefix='CTCG')
    # There should be only one feature for each position, specified by admin
    if len(position_features) != 1:
        msg = "Only one level -1 for this position should exist"
        raise RuntimeError(msg)
    else:
        minus_1_feature = position_features[0]
        owner = str(minus_1_feature.owner)
        if owner != "admin":
            msg = "This level -1 wasn't provided by admin"
            raise (RuntimeError(msg))
        else:
            gb_path = os.path.join(proj_settings.MEDIA_ROOT,
                                   minus_1_feature.genbank_file.name)
            part_record = SeqIO.read(gb_path, 'gb')
    return part_record, part_record.id


def _left_border_is_malformed(seq):
    left_border_features = []
    for feature in seq.features:
        try:
            if feature.qualifiers['label'][0] == "LB":
                left_border_features.append(feature)
        except:
            continue
    return len(left_border_features) > 1


def _fix_joined_seq(seq):
    features = []
    border_positions = {"five_prime": [], "three_prime": []}
    qualifiers = {}
    for feature in seq.features:

        if "label" not in feature.qualifiers:
            features.append(feature)
            continue

        if feature.qualifiers['label'][0] != "LB":
            features.append(feature)
        else:
            five_prime = feature.location.start
            three_prime = feature.location.end
            border_positions["five_prime"].append(five_prime)
            border_positions["three_prime"].append(three_prime)

            if not qualifiers:
                qualifiers = {key: value for key, value in feature.qualifiers.items()}
    left_border_five_prime = ExactPosition(min(border_positions["five_prime"]))
    left_border_three_prime = ExactPosition(max(border_positions["three_prime"]))
    location = FeatureLocation(left_border_five_prime,
                                          left_border_three_prime)
    features.append(SeqFeature(location, strand=1,
                                          type="repeat_region",
                                          qualifiers=qualifiers))
    seq.features = [feature for feature in features]
    return seq


def _get_sub_seq_index(sub_seq, seq):
    start_idx = seq.find(sub_seq)
    if start_idx == -1:
        msg = "Part not found"
        raise RuntimeError(msg)
    else:
        end_idx = start_idx + len(sub_seq)
    return start_idx, end_idx


def reverse_complement(seq):
    code = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(code[nucl] for nucl in seq[::-1])


def _parse_rebase_file(enzyme):
    restriction_site = parse_rebase_file(REBASE_FILE)[enzyme]
    if '^' in restriction_site:
        raise NotImplementedError
    forw_rec_site, cut_site = restriction_site.split('(')
    rev_rec_site = reverse_complement(forw_rec_site)
    forw_cut_delta, rev_cut_delta = cut_site.rstrip(')').split('/')
    forw_cut_delta, rev_cut_delta = int(forw_cut_delta), int(rev_cut_delta)
    return {forw_rec_site: forw_cut_delta, rev_rec_site: rev_cut_delta}


def _get_parts_from_minus_one_part(level_minus_one_part, parts, enzymes=None):
    minus_one_part, uniquename = _get_position_feature(level_minus_one_part)
    if enzymes is None:
        enzymes = ['BsmBI']

    rebase_file = _parse_rebase_file(enzymes[0])

    rec_sites = sorted([site for site in rebase_file.keys()])
    rec_sites_regex = '(' + '|'.join(rec_sites) + ')'
    rec_sites_regex = re.compile(rec_sites_regex, flags=re.IGNORECASE)
    splitted_parts = rec_sites_regex.split(str(minus_one_part.seq))

    # we need to look at the sequence as circular
    # so we are appending the end of the sequence to the beggining
    # this way we can find the first restriction site and the tRNA

    first_part = splitted_parts[-1] + splitted_parts[0]
    first_part = rec_sites_regex.split(str(first_part))
    splitted_parts_circular = first_part + splitted_parts[1:]

    if len(splitted_parts_circular) != 9:
        msg = "This level -1 doesn't have 4 restriction sites"
        raise RuntimeError(msg)
    else:
        trna = splitted_parts_circular[2]
        trna_left_rec_site = splitted_parts_circular[1]
        trna_right_rec_site = splitted_parts_circular[3]
        trna_start_idx, trna_end_idx = _get_sub_seq_index(trna,
                                                          minus_one_part.seq)

        trna_start_idx += rebase_file[trna_left_rec_site]
        trna_end_idx = trna_end_idx - rebase_file[trna_right_rec_site]
        trna_record = minus_one_part[trna_start_idx:trna_end_idx]


        

        scaffold = splitted_parts_circular[6]
        scaffold_left_rec_site = splitted_parts_circular[5]
        scaffold_right_rec_site = splitted_parts_circular[7]
        scaffold_start_idx, scaffold_end_idx = _get_sub_seq_index(scaffold,
                                                                  minus_one_part.seq)

        scaffold_start_idx += rebase_file[scaffold_left_rec_site]
        scaffold_end_idx = scaffold_end_idx - rebase_file[scaffold_right_rec_site]
        scaffold_record = minus_one_part[scaffold_start_idx:scaffold_end_idx]
    parts[TRNA]["seq"] = trna_record
    parts[TRNA]["uniquename"] = trna_record.id
    parts[SCAFFOLD]["seq"] = scaffold_record
    parts[SCAFFOLD]["uniquename"] = scaffold_record.id
    return parts


def _get_target_part(parts):
    target_uniquename = parts[CRISPR_MULTIPLEXING_TARGET]["uniquename"]
    target = Feature.objects.get(uniquename=target_uniquename)
    parts[CRISPR_MULTIPLEXING_TARGET]["feature"] = target
    return parts


def _get_vector_part(parts):
    feature = Feature.objects.get(uniquename=parts[VECTOR_TYPE_NAME]["uniquename"])
    parts[VECTOR_TYPE_NAME]["feature"] = feature
    return parts


def level_0_genbank(request):
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    part_types = PARTS_TO_ASSEMBLE['crispr_multiplexing']
    if VECTOR_TYPE_NAME not in part_types:
        part_types.append(VECTOR_TYPE_NAME)
    used_parts = OrderedDict()

    if request_data:
        for part_type in part_types:
            used_parts[part_type] = {"uniquename": request_data[part_type]}
        assembled_seq = assemble_unit_zero(used_parts, part_types)
        filename = assembled_seq.name + '.gb'
        response = HttpResponse(assembled_seq.format('genbank'),
                                content_type='text/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}"'.format(filename)
        return response
    return HttpResponseBadRequest()


def _get_parts_for_level_0(request):
    trna, trna_id = _get_position_feature(request.POST[u'tRNA'])
    target = request.POST[u'CRISPR Multiplexing Target']
    scaffold, scaffold_id = _get_position_feature(request.POST[u'scaffold'])
    vector = request.POST[u'Vector']
    return [trna_id, target, scaffold_id, vector]


def level_0_protocol(request):
    "it returns the protocol "
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    part_order = _get_parts_for_level_0(request)
    protocol = write_protocol_for_level_0(request.POST, 'CRISPR Multipartite', part_order)
    response = HttpResponse(protocol, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
    return response


def write_protocol_for_level_0(protocol_data, assembly_type, part_order):
    "it writes the protocol in a variable"
    protocol = []
    protocol.append("{0} Assembly Protocol".format(assembly_type))
    protocol.append("")

    fragments = []
    for part in part_order[0:-1]:
        part_name = part
        fragments.append(part_name)
    part_str = "({0}){1}".format(":".join(fragments),
                                 protocol_data[VECTOR_TYPE_NAME])

    protocol.append("Entities to assemble: {0}".format(part_str))
    protocol.append("Reaction should be performed as follows:")
    protocol.append("\t75 ng of {}".format(part_order[0]))
    protocol.append("\t2 ng of {}".format(part_order[1]))
    protocol.append("\t75 ng of {}".format(part_order[-1]))
    protocol.append("\t5-10u of {0}".format('BsmBI'))
    protocol.append("")
    protocol.append(u"\t3u of T4 ligase")
    protocol.append("")
    protocol.append(u"\t1 microlitre Ligase Buffer")
    protocol.append("")
    protocol.append(u"Final volume: 10 microlitre")
    protocol.append("")

    lline0 = "We use Promega T4 DNA ligase(M180B), NEB BsaI (R0535S or R0535L), "
    lline0 += "NEB BtgZI (R0703S) and fermentas BsmBI/Esp3I (ER0451). We haven't "
    lline0 += "tried other enzymes suppliers but they will problably work as well"
    protocol.append(lline0)
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


def level_0_sbol(request):
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    part_types = PARTS_TO_ASSEMBLE['crispr_multiplexing']
    part_types.append(VECTOR_TYPE_NAME)
    used_parts = OrderedDict()

    if request_data:
        for part_type in part_types:
            used_parts[part_type] = {"uniquename": request_data[part_type]}
        assembled_seq = assemble_unit_zero(used_parts, part_types)
        filename = assembled_seq.name + '.xml'
        response = HttpResponse(convert_to_sbol(assembled_seq),
                                content_type='xml/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}"'.format(filename)
        return response
    return HttpResponseBadRequest()


def level_0_add(request):
    context = {}
    context.update(csrf(request))
    request_data = request.POST
    if not request_data:
        context.update({'info': "Not enough data to add the feature"})
        return render(request, 'goldenbraid_info.html',
                      context=context)

    part_types = [item for item in PARTS_TO_ASSEMBLE['crispr_multiplexing']]
    if VECTOR_TYPE_NAME not in part_types:
        part_types.append(VECTOR_TYPE_NAME)
    used_parts = OrderedDict()

    if request_data:
        for part_type in part_types:
            used_parts[part_type] = {"uniquename": request_data[part_type]}
    assembled_seq = assemble_unit_zero(used_parts, part_types)

    name = request_data['name']
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(assembled_seq.id),
                                    suffix='.gb', dir=TMP_DIR)
    temp_fhand.write(assembled_seq.format('gb').encode())
    temp_fhand.flush()
    temp_fhand.seek(0)
    props = {'Description': [textwrap.fill(request_data['description'], width=70)],
             'Reference': [request_data['reference']]}

    type_name = request_data[TRNA]

    try:
        feature = add_feature(name=name, type_name=type_name,
                              vector=request_data['Vector'],
                              genbank=temp_fhand, props=props,
                              owner=request.user, is_public=False)

    except IntegrityError as error:
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render(request, 'feature_exists.html',
                          context={})
        else:
            return HttpResponseServerError()
    except RuntimeError as error:
        if 'standar GBCloning categories. You should uso Other 'in str(error):
            info = str(error) + '.\n'
            info += 'Can not use this tool to add this feature.\n'
            info += 'You can download the genbank and add using the Add Feature tool.\n'
            info += 'Sorry for the inconvenience\n'
            context.update({'info': info})
            return render(request, 'goldenbraid_info.html',
                          context=context)
    except Exception as error:
        return HttpResponseServerError(str(error))
    # if everithing os fine we show the just added feature
    return redirect(feature.url)



def _get_cas12_minus_one_parts(level_minus_one_part, enzyme=None):
    parts = []
    if enzyme is not None:
        rebase_file = _parse_rebase_file(enzyme)
        rec_sites = sorted([site for site in rebase_file.keys()])
        rec_sites_regex = '(' + '|'.join(rec_sites) + ')'
        rec_sites_regex = re.compile(rec_sites_regex, flags=re.IGNORECASE)
        nonrec_minus_one = rec_sites_regex.split(str(level_minus_one_part.seq))[2]
        start_idx, end_idx = _get_sub_seq_index(nonrec_minus_one, level_minus_one_part.seq)
        minus_one_record = level_minus_one_part[start_idx+1:end_idx-5]
    else:
        minus_one_record = level_minus_one_part

    rec_sites_regex = "NNNNNNNNNNNNNNNNNNNN"
    rec_sites_regex = re.compile(rec_sites_regex, flags=re.IGNORECASE)
    splitted_parts = rec_sites_regex.split(str(minus_one_record.seq))
    for part in splitted_parts:
        part_start_idx, part_end_idx = _get_sub_seq_index(part, level_minus_one_part.seq)
        part_record = level_minus_one_part[part_start_idx:part_end_idx]
        parts.append(part_record)
    return parts


def get_crispr_cas12_multiplexing_scaffold(category):
    print(category)
    scaffold = Feature.objects.filter(Q(type__name=category) & Q(featureperm__owner__username="admin") &
                                      Q(prefix__iexact="CTCG") & Q(suffix__iexact="TGAG"))
    print(scaffold)
    # There should be only one feature for each position, specified by admin
    if len(scaffold) != 1:
        msg = "Only one level -1 for this position should exist"
        raise RuntimeError(msg)
    else:
        minus_1_feature = scaffold[0]
        owner = str(minus_1_feature.owner)
        if owner != "admin":
            msg = "This level -1 wasn't provided by admin"
            raise (RuntimeError(msg))
        else:
            gb_path = os.path.join(proj_settings.MEDIA_ROOT,
                                   minus_1_feature.genbank_file.name)
            part_record = SeqIO.read(gb_path, 'gb')
    return part_record, part_record.id


def _merge_scaffold_and_targets(scaffold_parts, targets, category=None):
    joined_seq = SeqRecord(Seq('', alphabet=generic_dna))
    for scaffold_part in scaffold_parts:
        joined_seq += scaffold_part
        target_feature_start = len(joined_seq)
        if targets:
            joined_seq += targets.pop(0)
            target_feature_end = len(joined_seq)
            feature_location = FeatureLocation(target_feature_start, target_feature_end)
            feature = SeqFeature(feature_location, type="target", qualifiers={"note": category})
            joined_seq.features.append(feature)

    return joined_seq


def crispr_view_cas12_single_TU(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        form = Cas12aSingleTUForm(request_data)
        if form.is_valid():
            part_types = ["promoter", "target", "signal"]
            cleaned_data = form.cleaned_data
            used_parts = OrderedDict([("PROM DPolIII+DRCas12 (A1-A2-A3-B1-B2e)", cleaned_data["promoter"]),
                                      ("D CAS12a target (B3e-B4-B5e)", cleaned_data["target"]),
                                      ("3_prime processing (B6c-C1)", cleaned_data["signal"]),
                                      ("vector", cleaned_data["vector"])])
            multi_form_data = OrderedDict([("part_1", cleaned_data["promoter"]),
                                           ("part_2", cleaned_data["target"]),
                                           ("part_3", cleaned_data["signal"]),
                                           ("vector", cleaned_data["vector"])])
            order = "PROM DPolIII+DRCas12 (A1-A2-A3-B1-B2e):D CAS12a target (B3e-B4-B5e):3_prime processing (B6c-C1)"
            context.update({'used_parts': used_parts, 'multi_type': "free",
                            'part_types': part_types, "order": order, "Vector": cleaned_data["vector"],
                            "multiform_data": multi_form_data})
            return render(request, 'cas12_single_tu_results.html', context=context)


    form = Cas12aSingleTUForm()
    populate_cas12a_single_tu_form(form, request.user)

    context.update({'form': form})
    template = 'cas12a_single_tu_template.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def assemble_cas12_unit_zero(targets, category):
    targets_copy = targets.copy()
    scaffold, scaffold_id = get_crispr_cas12_multiplexing_scaffold(category)
    scaffold_parts = _get_cas12_minus_one_parts(scaffold)
    scaffold_seq = _merge_scaffold_and_targets(scaffold_parts, targets, category=category)
    scaffold_parts_non_rec = _get_cas12_minus_one_parts(scaffold, enzyme='BsmBI')
    scaffold_seq_non_rec = _merge_scaffold_and_targets(scaffold_parts_non_rec, targets_copy, category=category)
    vector = Feature.objects.get(Q(prefix__iexact="TGAG") & Q(suffix__iexact="CTCG") & Q(name__iexact="pUPD2"))
    vector_gb_path = os.path.join(proj_settings.MEDIA_ROOT, "genbank",
                                  vector.genbank_file.name.replace("genbank/", ""))
    vector_record = SeqIO.read(vector_gb_path, 'gb')
    if has_rec_sites(str(scaffold_seq_non_rec.seq)):
        raise RuntimeError
    
    vector_enzyme = vector.enzyme_in[0]
    vector_seq = vector.residues
    pref_idx, suf_idx = get_prefix_and_suffix_index(vector_seq, vector_enzyme)[:2]
    vector_sub_seq = vector_record[pref_idx:]
    vector_sub_seq += vector_record[:suf_idx]
    merged_seq = scaffold_seq_non_rec + vector_sub_seq

    try:
        count = Count.objects.get(name=ASSEMBLED_SEQ)
    except Count.DoesNotExist:
        count = Count.objects.create(name=ASSEMBLED_SEQ, value=1)
    next_value = count.next

    merged_seq.id = ASSEMBLED_SEQ + '_' + next_value
    merged_seq.name = merged_seq.id
    merged_seq.description = "({0}){1}".format(scaffold_id,
                                               vector.uniquename)
    scaffold_seq.name = merged_seq.name

    return scaffold_seq, merged_seq


def assemble_unit_zero(parts, part_types):
    joined_seq = SeqRecord(Seq('', alphabet=generic_dna))
    names = {'parts': []}
    position_part = parts[TRNA]["uniquename"]
    parts = _get_parts_from_minus_one_part(position_part, parts)
    parts = _get_target_part(parts)
    parts = _get_vector_part(parts)
    for part_type in part_types:
        print(part_type)
        if "feature" in parts[part_type]:
            feature = parts[part_type]["feature"]
            seq = Seq(feature.residues)
            gb_path = os.path.join(proj_settings.MEDIA_ROOT, "genbank",
                                   feature.genbank_file.name.replace("genbank/", ""))
            part_record = SeqIO.read(gb_path, 'gb')

        else:
            part_record = parts[part_type]["seq"]
            seq = part_record.seq

        if part_type == CRISPR_MULTIPLEXING_TARGET:
            pref_idx = 0
            suf_idx = len(seq) - len(feature.suffix)
            names['parts'].append(feature.uniquename)
        elif part_type in [TRNA, SCAFFOLD]:
            # suffix was removed when extracted from level -1 vector
            # so we must conserve all the remaining sequence
            pref_idx = 0
            suf_idx = len(seq)
            names['parts'].append(part_record.id)
        else:
            if part_type == VECTOR_TYPE_NAME:
                enzyme = feature.enzyme_in[0]
                names['vector'] = feature.uniquename
            pref_idx, suf_idx = get_prefix_and_suffix_index(seq, enzyme)[:2]

        if part_type == VECTOR_TYPE_NAME and feature.direction == REVERSE:
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
            # we have to annotate the trna keeping in mind that
            # first 5 nucleotides of target are actually trna
            # also preffix of trna must not be annotated
            if part_type == CRISPR_MULTIPLEXING_TARGET:
                start_feature = ExactPosition(8)
                end_feature = ExactPosition(len(joined_seq) - len(part_sub_seq) + 4)
                location = FeatureLocation(start_feature,
                                           end_feature)
                feature =(SeqFeature(location, strand=1,
                                     type="misc_region",
                                     qualifiers={"note": "tRNA"}))
                joined_seq.features.append(feature)

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

        # crispr targets are the unique parts that are not vectors
        # and does not have a vector
        if part.vector is None and part.type.name != VECTOR_TYPE_NAME:
            pref_idx = 0
            suf_idx = len(seq) - len(part.suffix)
            names['parts'].append(part.uniquename)
        else:
            if part.type.name == VECTOR_TYPE_NAME:
                enzyme = part.enzyme_in[0]
                names['vector'] = part.uniquename
            else:
                names['parts'].append(part.uniquename)
                enzyme = part.enzyme_out[0]
            print(enzyme)
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
            print(part.uniquename)
            if suf_idx >= pref_idx:
                print(pref_idx, suf_idx)
                part_sub_seq = part_record[pref_idx:suf_idx]
            else:
                part_sub_seq = part_record[pref_idx:]
                part_sub_seq += part_record[:suf_idx]
            joined_seq += part_sub_seq

    if _left_border_is_malformed(joined_seq):
        joined_seq = _fix_joined_seq(joined_seq)

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

def cas12a_multiplexing_view_genbank(request):
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        targets = eval(request_data["targets"])
        category = request_data["category"]
        seq_name = request_data["seq_name"]
        _, assembled_seq = assemble_cas12_unit_zero(targets, category)
        filename = seq_name + '.gb'
        response = HttpResponse(assembled_seq.format('genbank'),
                                content_type='text/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}"'.format(filename)
        return response
    return HttpResponseBadRequest()

def cas12a_multiplexing_view_sbol(request, multi_type=None):
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        targets = eval(request_data["targets"])
        category = request_data["category"]
        seq_name = request_data["seq_name"]
        _, assembled_seq = assemble_cas12_unit_zero(targets, category)
        filename = seq_name + '.xml'
        response = HttpResponse(convert_to_sbol(assembled_seq),
                                content_type='xml/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}"'.format(filename)
        return response
    return HttpResponseBadRequest()


def cas12a_multiplexing_view_add(request):
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        targets = eval(request_data["targets"])
        category = request_data["category"]
        seq_name = request_data["seq_name"]
        _, assembled_seq = assemble_cas12_unit_zero(targets, category)
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(assembled_seq.id),
                                    suffix='.gb', dir=TMP_DIR)
    temp_fhand.write(assembled_seq.format('gb').encode())
    temp_fhand.flush()
    temp_fhand.seek(0)
    props = {'Description': [textwrap.fill(request_data['description'], width=70)],
             'Reference': [request_data['reference']]}
    try:
        feature = add_feature(name=assembled_seq.id, type_name=category,
                              vector='pUPD2',
                              genbank=temp_fhand, props=props,
                              owner=request.user, is_public=False)

    except IntegrityError as error:
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render(request, 'feature_exists.html',
                          context={})
        else:
            return HttpResponseServerError()

    except RuntimeError as error:
        if 'standar GBCloning categories. You should uso Other 'in str(error):
            info = str(error) + '.\n'
            info += 'Can not use this tool to add this feature.\n'
            info += 'You can download the genbank and add using the Add Feature tool.\n'
            info += 'Sorry for the inconvenience\n'
            context.update({'info': info})
            return render(request, 'goldenbraid_info.html',
                          context=context.flatten())
    except Exception as error:
        return HttpResponseServerError(str(error))
    # if everything is fine we show the just added feature

    return redirect(feature.url)


def multipartite_view_genbank(request, multi_type=None):
    'view of the multipartite tool'
    print(multi_type)

    if multi_type is None:
        return render(request, 'multipartite_initial.html',
                      context={})
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
                                    content_type='text/plain')
            response['Content-Disposition'] = 'attachment; '
            response['Content-Disposition'] += 'filename="{0}"'.format(filename)
            return response
    return HttpResponseBadRequest()


def multipartite_view_sbol(request, multi_type=None):
    'view of the multipartite tool'

    if multi_type is None:
        return render(request, 'multipartite_initial.html',
                      context={})
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
            filename = assembled_seq.name + '.xml'
            response = HttpResponse(convert_to_sbol(assembled_seq),
                                    content_type='xml/plain')
            response['Content-Disposition'] = 'attachment; '
            response['Content-Disposition'] += 'filename="{0}"'.format(filename)
            return response
    return HttpResponseBadRequest()


def multipartite_view_add(request):
    context = {}
    context.update(csrf(request))
    request_data = request.POST
    if not request_data:
        context.update({'info': "Not enough data to add the feature"})
        return render(request, 'goldenbraid_info.html',
                      context=context)
    multi_type = request_data['category']
    allowed_categories = list(PARTS_TO_ASSEMBLE.keys()) + ['free']
    if multi_type is None or multi_type not in allowed_categories:
        context.update({'info': "Not enough data to add the feature"})
        return render(request, 'goldenbraid_info.html', context=context)
    if multi_type == 'free':
        part_types = request_data['order'].split(':')
        print(2)
        print(request_data['order'])
    else:
        part_types = [p[0] for p in PARTS_TO_ASSEMBLE[multi_type]]

    multi_data = {'Vector': request_data['Vector']}
    for part_type in part_types:
        multi_data[part_type] = request_data[part_type]
    assembled_seq = assemble_parts(multi_data, part_types)
    name = request_data['name']
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(assembled_seq.id),
                                    suffix='.gb', dir=TMP_DIR)
    temp_fhand.write(assembled_seq.format('gb').encode())
    temp_fhand.flush()
    temp_fhand.seek(0)
    props = {'Description': [textwrap.fill(request_data['description'], width=70)],
             'Reference': [request_data['reference']]}
    if multi_type == FUNGAL_TU_TYPE_NAME:
        type_name = FUNGAL_TU_TYPE_NAME
    elif multi_type == FUNGAL_KNOCK_OUT:
        type_name = FUNGAL_KNOCK_OUT
    else:
        type_name = TU_TYPE_NAME
    try:
        feature = add_feature(name=name, type_name=type_name,
                              vector=request_data['Vector'],
                              genbank=temp_fhand, props=props,
                              owner=request.user, is_public=False)

    except IntegrityError as error:
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render(request, 'feature_exists.html',
                          context={})
        else:
            return HttpResponseServerError()
    except RuntimeError as error:
        print(str(error))
        if 'standar GBCloning categories. You should uso Other 'in str(error):
            info = str(error) + '.\n'
            info += 'Can not use this tool to add this feature.\n'
            info += 'You can download the genbank and add using the Add Feature tool.\n'
            info += 'Sorry for the inconvenience\n'
            context.update({'info': info})
            return render(request, 'goldenbraid_info.html',
                          context=context)
    except Exception as error:
        return HttpResponseServerError(str(error))
    # if everithing os fine we show the just added feature
    return redirect(feature.url)

def multipartite_fungal_tu_view(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        form = FungalTUForm(request_data)
        if form.is_valid():
            used_parts = OrderedDict()
            multi_form_data = form.cleaned_data
            used_parts = {FUNGAL_PROM_5UTR: multi_form_data["promoter"],
                          FUNGAL_CDS: multi_form_data["cds"],
                          FUNGAL_UTR3_TERM: multi_form_data["ter"],
                          VECTOR_TYPE_NAME: multi_form_data["vector"]}
            context.update({'used_parts': used_parts, 'multi_type': FUNGAL_TU_TYPE_NAME})
            return render(request, 'multipartite_result.html', context=context)
    else:
        form = FungalTUForm()
        populate_fungal_tu_form(form, request.user)
    context.update({'form': form})
    template = "fungal_tu.html"
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def multipartite_view(request, multi_type=None):
    context = {}
    if multi_type is None:
        return render(request, 'multipartite_initial.html', context=context)
    elif multi_type not in PARTS_TO_ASSEMBLE.keys():
        return Http404
    context = {}
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
            context.update({'used_parts': used_parts, 'multi_type': multi_type})
            return render(request, 'multipartite_result.html', context=context)
    else:
        form = form_class()

    context.update({'form': form})
    template = 'multipartite_template.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def write_protocol(protocol_data, assembly_type, part_order):
    "it writes the protocol in a variable"
    protocol = []
    protocol.append("{0} Assembly Protocol".format(assembly_type.title()))
    protocol.append("")

    fragments = []
    for part_type in part_order:
        part_name = protocol_data[part_type]
        fragments.append(part_name)
    part_str = "({0}){1}".format(":".join(fragments),
                                 protocol_data[VECTOR_TYPE_NAME])

    protocol.append("Entities to assemble: {0}".format(part_str))
    protocol.append("Reaction should be performed as follows:")

    part_order.append(VECTOR_TYPE_NAME)
    for part_type in part_order:
        part_name = protocol_data[part_type]
        quantity = '2' if part_type in CRYSPER_TARGETS_TO_DOMESTICATE else "75"
        protocol.append("\t{} ng of {}".format(quantity, part_name))
    for enzyme in get_enzymes_for_protocol(protocol_data, part_order):

        protocol.append("\t5-10u of {0}".format(enzyme))
    protocol.append("")
    protocol.append(u"\t3u of T4 ligase")
    protocol.append("")
    protocol.append(u"\t1 microlitre Ligase Buffer")
    protocol.append("")
    protocol.append(u"Final volume: 10 microlitre")
    protocol.append("")

    lline0 = "We use Promega T4 DNA ligase(M180B), NEB BsaI (R0535S or R0535L), "
    lline0 += "NEB BtgZI (R0703S) and fermentas BsmBI/Esp3I (ER0451). We haven't "
    lline0 += "tried other enzymes suppliers but they will problably work as well"
    protocol.append(lline0)
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

def _btgzi_must_be_added(vector, part_order, protocol_data):
    if "omega" in vector.uniquename:
        for part in part_order:
            part_name = protocol_data[part]

            print(part, part_name)
            part = Feature.objects.get(uniquename=part_name)
            if part.type.name != VECTOR_TYPE_NAME and part.vector is None:
                continue
            elif part.type.name != VECTOR_TYPE_NAME:
                print(part.vector.uniquename)
                if part.vector.uniquename == "pUPD2":
                    return True
        return False
    else:
        return False


def get_enzymes_for_protocol(protocol_data, part_order):
    'it gets the necessary enzymes'
    vector = Feature.objects.get(uniquename=protocol_data[VECTOR_TYPE_NAME])
    vec_enzyme_in = vector.enzyme_in
    vec_enzyme_out = vector.enzyme_out
    enzymes = set(vec_enzyme_in)

    if _btgzi_must_be_added(vector, part_order, protocol_data):
        enzymes.add("BtgZI")

    for part_type in part_order:
        part_name = protocol_data[part_type]
        part = Feature.objects.get(uniquename=part_name)
        # this conditions are only met by crispr targets
        if part.type.name != VECTOR_TYPE_NAME and part.vector is None:
            continue
        enzyme_outs = part.enzyme_out
        if vec_enzyme_in not in enzyme_outs:
            for enzyme_out in enzyme_outs:
                if enzyme_out not in vec_enzyme_out:
                    enzymes.add(enzyme_out)
                    break
    print(enzymes)
    return list(enzymes)


def multipartite_protocol_view(request):
    "it returns the protocol "
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    part_order = [p[0] for p in PARTS_TO_ASSEMBLE[request.POST['multi_type']]]
    protocol = write_protocol(request.POST, 'multipartite', part_order)
    response = HttpResponse(protocol, content_type='text/plain')
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
    print(protocol_data, part_order)
    protocol = write_protocol(protocol_data, 'multipartite', part_order)
    response = HttpResponse(protocol, content_type='text/plain')
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
                                    content_type='text/plain')
            filename = assembled_seq.name + '.gb'
            response['Content-Disposition'] = 'attachment; '
            response['Content-Disposition'] += 'filename="{0}"'.format(filename)
            return response

    return HttpResponseBadRequest('There was an error in the assembly')


def multipartite_view_free_sbol(request):
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
            response = HttpResponse(convert_to_sbol(assembled_seq),
                                    content_type='xml/plain')
            filename = assembled_seq.name + '.xml'
            response['Content-Disposition'] = 'attachment; '
            response['Content-Disposition'] += 'filename="{0}"'.format(filename)
            return response

    return HttpResponseBadRequest('There was an error in the assembly')


def multipartite_view_free(request, form_num):
    context = {}
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
                last_suffix = last_feat.suffix.upper()
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
                    context.update({'used_parts': used_parts, 'multi_type': 'free',
                                    'post_data': form.cleaned_data, 'order': ":".join(part_order)})
                    return render(request, 'multipartite_free_result.html', context=context)

                else:
                    # add new_field
                    part_num = len(feats)
                    feats = Feature.objects.filter(prefix__iexact=last_suffix)
                    if user.is_staff:
                        pass
                    elif user.is_authenticated:
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
    content_type = None
    return render(request, template, context=context, content_type=content_type)
