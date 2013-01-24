'''
Created on 2013 urt 17

@author: peio
'''
from django.template.context import RequestContext
from django.shortcuts import render_to_response
from django.http import Http404
from django import forms

from goldenbraid.models import Feature
from goldenbraid.settings import DB
from django.forms.widgets import Select
from django.core.context_processors import csrf

PARTS_TO_ASSEMBLE = {'basic': [('PROM+UTR+ATG', 'GGAG', 'AATG'),
                               ('CDS', 'AATG', 'GCTT'),
                               ('TER', 'GCTT', 'CGCT')],
                     'secreted': [('PROM+UTR+ATG', 'GGAG', 'AATG'),
                                 ('SP', 'AATG', 'AGCC'),
                                 ('CDS', 'AGCC', 'GCTT'),
                                 ('TER', 'GCTT', 'CGCT')]
                     }


def _get_multipartite_form(multi_type):
    'It returns a form for the given multipartite'
    form_fields = {}

    for parts in PARTS_TO_ASSEMBLE[multi_type]:
        choices = []
        for feat in Feature.objects.using(DB).filter(type__name=parts[0],
                                                     prefix=parts[1],
                                                     suffix=parts[2]):
            choices.append((feat.uniquename, feat.uniquename))

        name = parts[0]
        form_fields[name] = forms.CharField(max_length=100,
                                            widget=Select(choices=choices))
    form = type('MultiPartiteForm', (forms.BaseForm,),
                {'base_fields': form_fields})
    return form


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

    if request_data:
        form = _get_multipartite_form(multi_type, request_data)
    else:
        form = _get_multipartite_form(multi_type)

    context['form'] = form

    template = 'multipartite_template.html'
    mimetype = None
    return render_to_response(template, context, mimetype=mimetype)

