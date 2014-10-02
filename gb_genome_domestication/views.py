'''
Created on 2014 uzt 16

@author: peio
'''
from django.shortcuts import render_to_response, redirect
from django.template.context import RequestContext
from django.core.context_processors import csrf
from django.core.urlresolvers import reverse
# from django.views.decorators.csrf import csrf_exempt

from gb_genome_domestication.models import Feature
from gb_genome_domestication.settings import DB

from restcmd_client.views.tool import run_tool


def search_view(request):
    url = reverse(run_tool, kwargs={'cmd': 'blastplus'})
    return redirect(url)


def feature_view(request, uniquename):
    'The feature view'
    context = RequestContext(request)
    context.update(csrf(request))
    try:
        feature = Feature.objects.using(DB).get(uniquename=uniquename)
    except Feature.DoesNotExist:
        feature = None
    context['feature'] = feature
    return render_to_response('genome_domest_feature_template.html',
                              context)
