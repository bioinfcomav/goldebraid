from django.shortcuts import  render_to_response
from gb_genome_domestication.models import Feature
from gb_genome_domestication.settings import DB
# Create your views here.


def feature_view(request, uniquename):
    'The feature view'
    try:
        feature = Feature.objects.using(DB).get(uniquename=uniquename)
    except Feature.DoesNotExist:
        feature = None
    return render_to_response('genome_domest_feature_template.html',
                              {'feature': feature})

