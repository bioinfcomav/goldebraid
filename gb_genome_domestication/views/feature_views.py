from django.shortcuts import render_to_response
from gb_genome_domestication.models import Feature
from gb_genome_domestication.settings import DB
from django.template.context import RequestContext
from django.core.context_processors import csrf
# Create your views here.


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
