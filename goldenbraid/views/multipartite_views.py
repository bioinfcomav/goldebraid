'''
Created on 2013 urt 17

@author: peio
'''
from django.template.context import RequestContext
from django.shortcuts import render_to_response
from django.http import Http404


PARTS_TO_ASSEMBLE = {'basic':['PROM+UTR+ATG', 'CDS', 'TER'],
                     'secreted':['PROM+UTR+ATG', 'CDS', 'TER']}


def multipartite_view(request, multi_type=None):
    'view of the multipartite tool'
    if multi_type is None:
        return render_to_response('multipartite_initial.html', {},
                                  context_instance=RequestContext(request))
    elif multi_type not in PARTS_TO_ASSEMBLE.keys():
        return Http404
