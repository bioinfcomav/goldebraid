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

from collections import OrderedDict
from os.path import join
from tempfile import NamedTemporaryFile
from django.template.context_processors import csrf
from django.template import RequestContext
from django.contrib.auth.decorators import login_required
from django.http.response import HttpResponseServerError, HttpResponse
from django.shortcuts import render, redirect
from django.db.utils import IntegrityError
from goldenbraid.views.feature import add_feature
from goldenbraid.forms.crispr_multiplexing import (LocationForm,
                                                   TaxaChoiceForm,
                                                   SelectSequenceForm,
                                                   RegulationLocationForm,
                                                   RegulationSequenceForm,
                                                   VectorSelectForm,
                                                   PromoterChoiceForm,
                                                   ConfigurationChoiceForm,
                                                   AutoMultiTypeForm,
                                                   AutoMultiTargetsForm,
                                                   populate_location_form,
                                                   populate_sequence_form,
                                                   populate_regulation_location_form,
                                                   populate_regulation_sequence_form,
                                                   populate_vector_select_form,
                                                   populate_promoter_choice_form,
                                                   populate_configuration_choice_form,
                                                   populate_auto_type_form,
                                                   AutoMultiTUForm)
from goldenbraid.models import Feature
from django.db.models import Q
from goldenbraid.settings import (CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE,
                                  PARTS_TO_ASSEMBLE, AUTO_CRISPR_SELECTION_CHOICES,
                                  CRYSPER_CATEGORIES, AUTO_CRISPR_SELECTION_POSITIONS_DOMEST,
                                  CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO)

from goldenbraid.tags import (MONOCOT_TAXA, DICOT_TAXA,
                              CRISPR_MULTIPLEXING_TARGET,
                              TRNA, SCAFFOLD,
                              VECTOR_TYPE_NAME,
                              PROM_DICOT,
                              PROM_MONOCOT, TU_TYPE_NAME)
from goldenbraid.views.domestication import CRYSPER_PROTOCOL
from goldenbraid.domestication import domestication_crispr
from goldenbraid.views.multipartite import (_get_position_feature, assemble_unit_zero, assemble_parts,
                                            write_protocol_for_level_0,  write_protocol)


TAXAS = [MONOCOT_TAXA, DICOT_TAXA]
LOCATIONS = [CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[item][0] for item in CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE]
TMP_DIR = join("/home/golden/devel/gbdb", "files")


def _get_feature_name(position):
    # TEMPORARY BYPASS
    print(position)

    position_features = Feature.objects.filter(type__name=position).filter(featureperm__owner__username="admin")
    position_features = position_features.filter(prefix='CTCG')
    print(position_features)
    for position_feature in position_features:
        print(position_feature.uniquename)
    # There should be only one feature for each position, specified by admin
    if len(position_features) != 1:
        msg = "Only one level -1 for this position should exist"
        raise RuntimeError(msg)
    else:
        minus_1_feature_name = position_features[0].uniquename

    return minus_1_feature_name


def _get_targets(user):
    targets = Feature.objects.filter(type__name=CRISPR_MULTIPLEXING_TARGET)
    if user.is_authenticated():
        targets = targets.filter(Q(featureperm__owner__username=user) |
                                 Q(featureperm__is_public=True))
    return targets

def cripsr_view_cas9_auto_omega_protocol(request):
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    request_data = request.POST
    protocol = write_protocol(request_data, "bipartite", ['part_1', 'part_2'])
    response = HttpResponse(protocol, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
    return response

def cripsr_view_cas9_auto_omega_add(request):
    context = {}
    context.update(csrf(request))
    user = request.user
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        vector = request_data["Vector"]
        seq = assemble_parts(request_data, ['part_1', 'part_2'])
        props = {'Description': [textwrap.fill(request_data['description'], width=70)],
                 'Reference': [request_data['reference']]}
        name = request_data['name']
        temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(seq.id), suffix='.gb', dir=TMP_DIR)
        temp_fhand.write(seq.format('gb').encode())
        temp_fhand.flush()
        temp_fhand.seek(0)
        try:
            feature = add_feature(name=name, type_name=TU_TYPE_NAME,
                                  vector=vector, genbank=temp_fhand, props=props,
                                  owner=request.user, is_public=False)
            protocol = write_protocol(request_data, "bipartite", ['part_1', 'part_2'])
        except IntegrityError as error:
            if 'feature already in db' in str(error):
                # TODO choose a template
                return render_to_response('feature_exists.html',
                                          {}, context_instance=RequestContext(request))
            else:
                return HttpResponseServerError()
        except Exception as error:
            return HttpResponseServerError(error)

        return redirect(feature.url)

@login_required
def crispr_view_cas9_multiplexing_auto_add_tu(request):
    context = {}
    context.update(csrf(request))
    user = request.user
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        form = AutoMultiTUForm(request_data)
        if form.is_valid():
            cleaned_data = form.cleaned_data
            vector = cleaned_data['vector']
            auto_tu = request_data['auto_TU']
            transcription_unit = cleaned_data['transcription_unit']
            context = {'Vector': vector, 'part_1': transcription_unit,
                       'part_2': auto_tu}
            return render(request, "cas9_auto_results.html", context=context)
        

def crispr_view_cas9_multiplexing_auto_protocol(request):
    merged_protocols = ""
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    protocols = eval(request.POST['protocols'])
    for protocol in protocols:
        merged_protocols += protocol
    response = HttpResponse(merged_protocols, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
    return response


@login_required
def get_already_created_seq(request, category, seq,  prefix, suffix):
    query = Feature.objects.filter(Q(type__name__icontains=category) & 
                                   Q(residues__iexact=str(seq.seq)) & 
                                   Q(prefix__iexact=prefix) &
                                   Q(suffix__iexact=suffix) &
                                   Q(featureperm__owner__username=request.user))
    if len(query) > 0:
        return query[0]
    else:
        return None


@login_required
def crispr_view_cas9_multiplexing_auto(request, section):
    protocols = []
    target_protocols = []
    level_0_protocols = []
    context = OrderedDict()
    context.update(csrf(request))
    user = request.user
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    form = None
    if section is None:
        form = AutoMultiTypeForm()
        populate_auto_type_form(form)
        context['section'] = 'select_type'
    elif section == 'select_type':
        form = AutoMultiTypeForm(request_data)
        if form.is_valid():
            cleaned_data = form.cleaned_data
            multiplexing_type = cleaned_data["multiplexing_type"]
            context['section'] = multiplexing_type
            form = AutoMultiTargetsForm(category=multiplexing_type)

    elif section in list(AUTO_CRISPR_SELECTION_CHOICES.keys()):
        form = AutoMultiTargetsForm(request_data, category=section)
        if form.is_valid():
            cleaned_data = form.cleaned_data
            target_prefix = 'GTGC'
            target_suffix = 'GTTT'
            for part in PARTS_TO_ASSEMBLE[section][1:]:
                category = part[0]
                target = cleaned_data[category]
                new_seq = domestication_crispr(target, CRISPR_MULTIPLEXING_TARGET, target_prefix, target_suffix)
                forw_dim = str(new_seq.seq[:-4])
                rev_dim = str(new_seq[4:].reverse_complement().seq)
                feature = get_already_created_seq(request,CRISPR_MULTIPLEXING_TARGET, new_seq, 
                                                       target_prefix, target_suffix)
                if feature is None:
                    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(new_seq.id),
                                                    suffix='.gb', dir=TMP_DIR)
                    
                    temp_fhand.write(new_seq.format('gb').encode())
                    temp_fhand.flush()
                    temp_fhand.seek(0)
                    props = {'Description': ["auto_domestication, {}".format(category)],
                             'Reference': ['']}
                    try:
                        feature_category = CRISPR_MULTIPLEXING_TARGET
                        feature = add_feature(name="auto_{}".format(new_seq.id), type_name=CRISPR_MULTIPLEXING_TARGET,
                                              vector=None, genbank=temp_fhand, props=props,
                                              owner=request.user, is_public=False,
                                              prefix=target_prefix, suffix=target_suffix)
                    except IntegrityError as error:
                        if 'feature already in db' in str(error):
                        # TODO choose a template
                            return render(request, 'feature_exists.html', context=context)
                        else:
                            return HttpResponseServerError()
                    except Exception as error:
                        raise
                        return (error)

                target_feature = feature.uniquename
      
                if "domesticated_targets" in context:
                    context["domesticated_targets"][category] = feature.uniquename
                else:
                    context["domesticated_targets"] = OrderedDict({category: feature.uniquename})

                protocol = "Domesticated target {} {}\n\n".format(category, feature.uniquename) + CRYSPER_PROTOCOL.format(forw_primer=forw_dim, rev_primer=rev_dim)
                target_protocols.append(protocol)
                protocol = ""
                
                part_types = [item for item in PARTS_TO_ASSEMBLE['crispr_multiplexing']]
                part_types.append(VECTOR_TYPE_NAME)
                _, scaffold_id = _get_position_feature(category)
                vector = Feature.objects.get(Q(prefix__iexact="TGAG") & Q(suffix__iexact="CTCG") & Q(name__iexact="pUPD2"))
                used_parts = OrderedDict([('tRNA', {'uniquename': category}),
                                          (CRISPR_MULTIPLEXING_TARGET, {"uniquename": feature.uniquename}),
                                          ('scaffold', {'uniquename': category}),
                                          (VECTOR_TYPE_NAME, {'uniquename': vector.uniquename})])
                assembled_seq = assemble_unit_zero(used_parts, part_types)
                prefix = CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[category][1]
                suffix = CRYSPR_MULTIPLEX_CATEGORIES_LEVEL_ZERO[category][2]
                feature = get_already_created_seq(request, category, assembled_seq, prefix, suffix)
                
                if feature is None:
                    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(assembled_seq.id),
                                                    suffix='.gb', dir=TMP_DIR)
                    temp_fhand.write(assembled_seq.format('gb').encode())
                    temp_fhand.flush()
                    temp_fhand.seek(0)
                    props = {'Description': ["auto_level0, {}".format(category)],
                             'Reference': [""]}

                    try:
                        feature = add_feature(name="auto_{}".format(assembled_seq.id), type_name=category,
                                              vector=vector.uniquename,
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
                    
                if "level0_sequences" in context:
                    context["level0_sequences"][category] = feature.uniquename

                else:
                    context["level0_sequences"] = OrderedDict({category: feature.uniquename})
                protocol_data = {"Vector": vector.uniquename}
                _, position_id = _get_position_feature(category)
                part_order = [position_id, target_feature, vector.uniquename]
                protocol = write_protocol_for_level_0(protocol_data,
                                                      category,
                                                      part_order)
                level_0_protocols.append("\n\n"+ "Assembly of {}\n".format(feature.uniquename) +protocol)
                protocol = ""

            part_types = [p[0] for p in PARTS_TO_ASSEMBLE[section]]
            multi_data = {'Vector': 'pDGB3_alpha2'}
            if "dicot" in section:
                multi_data[PROM_DICOT] = "GB1001"
            elif "monocot" in section:
                multi_data[PROM_MONOCOT] = "GB1184"
            for part_type in part_types[1:]:
                multi_data[part_type] = context['level0_sequences'][part_type]
            assembled_seq = assemble_parts(multi_data, part_types)
            #We only do the assembly with pDGB3_alpha2
            prefix = "GTCA"
            suffix = "CGCT" 
            category = TU_TYPE_NAME
            feature =  get_already_created_seq(request, category, assembled_seq, prefix, suffix)
            if feature is None:
                temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(assembled_seq.id),
                                                suffix='.gb', dir=TMP_DIR)
                temp_fhand.write(assembled_seq.format('gb').encode())
                temp_fhand.flush()
                temp_fhand.seek(0)
                props = {'Description': ["auto_TU, {}".format(assembled_seq.id)],
                         'Reference': [""]}
                try:
                    feature = add_feature(name="auto_{}".format(assembled_seq.id), type_name=category,
                                          vector=multi_data["Vector"],
                                          genbank=temp_fhand, props=props,
                                          owner=request.user, is_public=False)
                except IntegrityError as error:
                    print("already")
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
                        info += 'Sorry for the inconvenience\n'
                        context.update({'info': info})
                        return render(request, 'goldenbraid_info.html',
                                      context=context)
                    elif "rec site" in str(error):
                        info = "The final assembly have restriction sites" + '.\n'
                        info += 'Can not use this tool to add this feature.\n'
                        info += 'Sorry for the inconvenience\n'
                    context.update({'info': info})
                    return render(request, 'goldenbraid_info.html',
                                  context=context)

                except Exception as error:
                    return HttpResponseServerError(str(error))
            
            context["auto_TU"] = feature.uniquename
            protocol = write_protocol(multi_data, 'multipartite', part_types[:-1])
            tu_protocol = "\n\nAssembly of {}\n".format(feature.uniquename) + protocol
            protocols = target_protocols + level_0_protocols 
            protocols.append(tu_protocol)
            protocol = ""

            form = AutoMultiTUForm(request_data)
            context['form'] = form
            context['protocols'] = protocols
            context['section'] = None
            return render(request, 'cas9_auto.html',
                          context=context)

    context['form'] = form
    template = 'crispr_multiplex_auto.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def level_0_editing(request, section):
    context = {}
    context.update(csrf(request))
    user = request.user
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    form = None
    if section is None:
        form = TaxaChoiceForm()
        context['section'] = 'Position_Choice'
    elif section == 'Position_Choice':
        form = TaxaChoiceForm(request_data)
        if form.is_valid():
            cleaned_data = form.cleaned_data
            taxa = cleaned_data["taxa_choice"]
            context['section'] = taxa
            form = LocationForm()
            populate_location_form(form, taxa)

    elif section in TAXAS:
        form = LocationForm(request_data)
        if form.is_valid():
            cleaned_data = form.cleaned_data
            taxa = cleaned_data['taxa']
            position = cleaned_data['position']
            context['section'] = cleaned_data['position']
            form = SelectSequenceForm()
            populate_sequence_form(form, taxa, position, user)

    elif section in LOCATIONS:
            form = SelectSequenceForm(request_data)
            if form.is_valid():
                cleaned_data = form.cleaned_data
                position_part = cleaned_data["position"]

                target_part = cleaned_data["target"]
                part_types = PARTS_TO_ASSEMBLE['crispr_multiplexing']
                used_parts = OrderedDict()
                parts_name = OrderedDict()
                for part_type in part_types:
                    if part_type in [TRNA, SCAFFOLD]:
                        used_parts[part_type] = position_part
                        parts_name[part_type] = _get_feature_name(position_part)
                    elif part_type == CRISPR_MULTIPLEXING_TARGET:
                        used_parts[part_type] = target_part
                        parts_name[part_type] = target_part
                used_parts.update({VECTOR_TYPE_NAME: u'pUPD2'})
                context.update({"used_parts": used_parts, "multi_type": "crispr_multiplexing",
                                "parts_name": parts_name})
                return render(request, 'level0_result.html', context=context)

    context['form'] = form
    template = 'crispr_level0_editing.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)



def level_0_regulation(request, section, mode):
    context = {}
    context.update(csrf(request))
    user = request.user
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    form = None

    if section is None:
        form = RegulationLocationForm()
        populate_regulation_location_form(form, mode)
        context['section'] = 'Position_Choice'
    elif section == 'Position_Choice':
        form = RegulationLocationForm(request_data)
        if form.is_valid():
            cleaned_data = form.cleaned_data
            position = cleaned_data["position"]
            context['section'] = position
            form = RegulationSequenceForm()
            populate_regulation_sequence_form(form, position, user)
    elif section in LOCATIONS:
            form = RegulationSequenceForm(request_data)
            if form.is_valid():
                cleaned_data = form.cleaned_data
                position_part = cleaned_data["position"]
                target_part = cleaned_data["target"]
                part_types = PARTS_TO_ASSEMBLE['crispr_multiplexing']
                used_parts = OrderedDict()
                parts_name = OrderedDict()
                for part_type in part_types:
                    if part_type in [TRNA, SCAFFOLD]:
                        used_parts[part_type] = position_part
                        parts_name[part_type] = _get_feature_name(position_part)
                    elif part_type == CRISPR_MULTIPLEXING_TARGET:
                        used_parts[part_type] = target_part
                        parts_name[part_type] = target_part
                used_parts.update({VECTOR_TYPE_NAME: u'pUPD2'})
                context.update({"used_parts": used_parts, "parts_name": parts_name, 
                                "multi_type": "crispr_multiplexing"})
                return render(request, 'level0_result.html', context=context)

    context['form'] = form
    template = 'crispr_level0_regulation.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)




@login_required
def crispr_multiplexing(request, section):
    context = {}
    context.update(csrf(request))
    user = request.user
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    form = None

    if section is None:
        form = VectorSelectForm()
        populate_vector_select_form(form, user)
        context['section'] = 'Promoter_Choice'

    elif section == 'Promoter_Choice':
        form = VectorSelectForm(request_data)
        if form.is_valid():
            cleaned_data = form.cleaned_data
            vector = cleaned_data["vector"]
            context['section'] = 'Config_choice'
            form = PromoterChoiceForm()
            populate_promoter_choice_form(form, vector)

    elif section == 'Config_choice':
        form = PromoterChoiceForm(request_data)
        if form.is_valid():
            cleaned_data = form.cleaned_data
            vector = cleaned_data['vector']
            promoter = cleaned_data['promoter']
            form = ConfigurationChoiceForm()
            populate_configuration_choice_form(form, vector, promoter)


    context['form'] = form
    template = 'crispr_multiplexing.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)
