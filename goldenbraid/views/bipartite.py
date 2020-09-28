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
from django.contrib.auth.decorators import login_required
from django.core.files.temp import NamedTemporaryFile
from django.db.utils import IntegrityError
from django.http.response import HttpResponseServerError
from django.template.context import RequestContext
from django.template.context_processors import csrf
from django.shortcuts import render, redirect
from django.http import HttpResponse, HttpResponseBadRequest

from goldenbraid.views.feature import add_feature
from goldenbraid.views.multipartite import assemble_parts, write_protocol
from goldenbraid.tags import VECTOR_TYPE_NAME, MODULE_TYPE_NAME, FUNGAL_MODULE_TYPE_NAME
from goldenbraid.forms.assemblers import (BipartiteForm1, BipartiteForm2,
                                          get_part2_choices, BipartiteForm3,
                                          get_bipart_vector_choices,
                                          get_part1_choice,
                                          FungalBipartiteForm1,
                                          FungalBipartiteForm2,
                                          FungalBipartiteForm3,
                                          fungal_get_part1_choice,
                                          fungal_get_part2_choices,
                                          fungal_get_bipart_vector_choices)
from goldenbraid.forms.autocomplete_forms import Bipart1Form, Bipart2Form, Bipart3Form
from goldenbraid.sbol import convert_to_sbol

TMP_DIR = join("/home/golden/devel/gbdb", "files")


def fungal_bipartite_view(request, form_num):
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
        form = FungalBipartiteForm1()
        context['form_num'] = '1'
        form.fields['part_1'].widget.choices = fungal_get_part1_choice(request.user)
    elif form_num == '1':
        if request_data:
            form = FungalBipartiteForm1(request_data)
            if form.is_valid():
                form1_data = form.cleaned_data
                form = FungalBipartiteForm2()
                form.fields['part_1'].initial = form1_data['part_1']
                choices_part2 = fungal_get_part2_choices(form1_data['part_1'],
                                                         request.user)
                form.fields['part_2'].widget.choices = choices_part2
                context['form_num'] = '2'
            else:
                context['form_num'] = '1'
    elif form_num == '2':
        if request_data:
            form = FungalBipartiteForm2(request_data)
            if form.is_valid():
                form2_data = form.cleaned_data
                form = FungalBipartiteForm3()
                form.fields['part_1'].initial = form2_data['part_1']
                form.fields['part_2'].initial = form2_data['part_2']
                choices_vector = fungal_get_bipart_vector_choices(form2_data['part_1'],
                                                                  request.user)
                form.fields['Vector'].widget.choices = choices_vector
                context['form_num'] = '3'
    elif form_num == '3':
        if request_data:
            form = FungalBipartiteForm3(request_data)
            if form.is_valid():
                used_parts = OrderedDict()
                cleaned_data = form.cleaned_data
                used_parts['part_1'] = cleaned_data['part_1']
                used_parts['part_2'] = cleaned_data['part_2']
                used_parts[VECTOR_TYPE_NAME] = cleaned_data[VECTOR_TYPE_NAME]
                return render(request, 'fungal_bipartite_result.html',
                              context={'used_parts': used_parts})
            else:
                context['form_num'] = '3'
    if form is None:
        form = FungalBipartiteForm1()
        context['form_num'] = '1'
    context['form'] = form
    template = 'fungal_bipartite_template.html'
    content_type = None
    return render(request, template, context=context.flatten(),
                  content_type=content_type)


# def bipartite_view(request, form_num):
#     context = RequestContext(request)
#     context.update(csrf(request))
#     if request.method == 'POST':
#         request_data = request.POST
#     elif request.method == 'GET':
#         request_data = request.GET
#     else:
#         request_data = None
#     form = None
#     if form_num is None:
#         form = BipartiteForm1()
#         context['form_num'] = '1'
#         form.fields['part_1'].widget.choices = get_part1_choice(request.user)
#     elif form_num == '1':
#         if request_data:
#             form = BipartiteForm1(request_data)

#             if form.is_valid():
#                 form1_data = form.cleaned_data
#                 form = BipartiteForm2()
#                 form.fields['part_1'].initial = form1_data['part_1']
#                 choices_part2 = get_part2_choices(form1_data['part_1'],
#                                                   request.user)
#                 form.fields['part_2'].widget.choices = choices_part2
#                 context['form_num'] = '2'
#             else:
#                 context['form_num'] = '1'
#     elif form_num == '2':
#         if request_data:
#             form = BipartiteForm2(request_data)
#             if form.is_valid():
#                 form2_data = form.cleaned_data
#                 form = BipartiteForm3()
#                 form.fields['part_1'].initial = form2_data['part_1']
#                 form.fields['part_2'].initial = form2_data['part_2']
#                 choices_vector = get_bipart_vector_choices(form2_data['part_1'],
#                                                            request.user)
#                 form.fields['Vector'].widget.choices = choices_vector
#                 context['form_num'] = '3'
#     elif form_num == '3':
#         if request_data:
#             form = BipartiteForm3(request_data)
#             if form.is_valid():
#                 used_parts = OrderedDict()
#                 cleaned_data = form.cleaned_data
#                 used_parts['part_1'] = cleaned_data['part_1']
#                 used_parts['part_2'] = cleaned_data['part_2']
#                 used_parts[VECTOR_TYPE_NAME] = cleaned_data[VECTOR_TYPE_NAME]
#                 return render(request, 'bipartite_result.html',
#                               context={'used_parts': used_parts})
#             else:
#                 context['form_num'] = '3'
#     if form is None:
#         form = BipartiteForm1()
#         context['form_num'] = '1'
#     context['form'] = form
#     template = 'bipartite_template.html'
#     content_type = None
#     return render(request, template, context=context.flatten(),
#                   content_type=content_type)


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
        form = Bipart1Form()
        context['form_num'] = '1'
    elif form_num == '1':
        if request_data:
            form = Bipart1Form(request_data)
            if form.is_valid(): 
                form1_data = form.cleaned_data
                
                form = Bipart2Form(initial={'part1': form1_data['uniquename']})
                context['form_num'] = '2'
            else:
                context['form_num'] = '1'
    elif form_num == '2':
        if request_data:
            form = Bipart2Form(request_data)
            if form.is_valid():
                form2_data = form.cleaned_data
                form = Bipart3Form(initial={'part2': form2_data['uniquename'], 'part1': form2_data['part1']})
                context['form_num'] = '3'
    elif form_num == '3':
        if request_data:
            form = Bipart3Form(request_data)
            if form.is_valid():
                used_parts = OrderedDict()
                cleaned_data = form.cleaned_data
                used_parts['part_1'] = cleaned_data['part1']
                used_parts['part_2'] = cleaned_data['part2']
                used_parts[VECTOR_TYPE_NAME] = cleaned_data['uniquename']
                return render(request, 'bipartite_result.html',
                              context={'used_parts': used_parts})
            else:
                context['form_num'] = '3'
    if form is None:
        form = Bipart1Form()
        context['form_num'] = '1'
    context['form'] = form
    template = 'bipartite_autocomplete.html'
    content_type = None
    return render(request, template, context=context.flatten(),
                  content_type=content_type)


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
            response = HttpResponse(seq.format('genbank'),
                                    content_type='text/plain')
            filename = seq.name + '.gb'
            response['Content-Disposition'] = 'attachment; '
            response['Content-Disposition'] += 'filename="{0}"'.format(filename)
            return response
    return HttpResponseBadRequest()


def bipartite_view_sbol(request):
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
            response = HttpResponse(convert_to_sbol(seq),
                                    content_type='xml/plain')
            filename = seq.name + '.xml'
            response['Content-Disposition'] = 'attachment; '
            response['Content-Disposition'] += 'filename="{0}"'.format(filename)
            return response
    return HttpResponseBadRequest()


def bipartite_view_protocol(request):
    "it returns the protocol "
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    protocol = write_protocol(request.POST, "bipartite", ['part_1', 'part_2'])
    response = HttpResponse(protocol, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
    return response


@login_required
def bipartite_view_add(request):
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST
    if not request_data:
        return render(request, 'goldenbraid_info.html',
                               {'info': "Not enougth data to add the feature"})
    name = request_data['name']
    vector_name = request_data['Vector']
    used_parts = {'Vector': vector_name,
                  'part_1': request_data['part_1'],
                  'part_2': request_data['part_2']}
    seq = assemble_parts(used_parts, ['part_1', 'part_2'])
    props = {'Description': [textwrap.fill(request_data['description'])],
             'Reference': [request_data['reference']]}
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(seq.id), suffix='.gb', dir=TMP_DIR)
    temp_fhand.write(seq.format('gb').encode())
    temp_fhand.flush()
    temp_fhand.seek(0)
    if "fungal" in request_data:
        type_name = FUNGAL_MODULE_TYPE_NAME
    else:
        type_name = MODULE_TYPE_NAME
    try:
        feature = add_feature(name=name, type_name=type_name,
                              vector=vector_name, genbank=temp_fhand, props=props,
                              owner=request.user, is_public=False)

    except IntegrityError as error:
        print(error)
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render_to_response('feature_exists.html',
                                      {}, context_instance=RequestContext(request))
        else:
            return HttpResponseServerError()
    except Exception as error:
        print(error)
        return HttpResponseServerError(error)
    # if everithing os fine we show the just added feature
    return redirect(feature.url)
