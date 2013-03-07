'''
Created on 2013 ots 5

@author: peio
'''

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from django.template.context import RequestContext
from django.core.context_processors import csrf
from django.shortcuts import render_to_response
from django.http import HttpResponse, HttpResponseBadRequest

from goldenbraid.views.multipartite_views import assemble_parts, write_protocol
from goldenbraid.tags import VECTOR_TYPE_NAME
from goldenbraid.forms import (BipartiteForm1, BipartiteForm2,
                               get_part2_choices, BipartiteForm3,
                               get_bipart_vector_choices)


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
                choices_part2 = get_part2_choices(form1_data['part_1'])
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
                choices_vector = get_bipart_vector_choices(form2_data['part_1'])
                form.fields['Vector'].widget.choices = choices_vector
                context['form_num'] = '3'
    elif form_num == '3':
        if request_data:
            form = BipartiteForm3(request_data)
            if form.is_valid():
                used_parts = OrderedDict()
                cleaned_data = form.cleaned_data
                used_parts['part_1'] = cleaned_data['part_1']
                used_parts['part_2'] = cleaned_data['part_2']
                used_parts[VECTOR_TYPE_NAME] = cleaned_data[VECTOR_TYPE_NAME]
                return render_to_response('bipartite_result.html',
                                          {'used_parts': used_parts},
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
            response = HttpResponse(seq.format('genbank'),
                                    mimetype='text/plain')
            response['Content-Disposition'] = 'attachment; '
            response['Content-Disposition'] += 'filename="assembled_seq.gb"'
            return response
    return HttpResponseBadRequest()


def bipartite_view_protocol(request):
    "it returns the protocol "
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    protocol = write_protocol(request.POST, "bipartite", ['part_1', 'part_2'])
    response = HttpResponse(protocol, mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
    return response
