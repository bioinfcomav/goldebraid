'''
Created on 2015 api. 1

@author: peio
'''
import json
from django.http.response import HttpResponse

from goldenbraid.models import Feature


def api_feature_uniquenames_view(request):
    query = Feature.objects.all()

    if request.method == 'GET':
        if u'term' in request.GET:
            term = request.GET['term']
            query = query.filter(uniquename__contains=term)
        if u'limit' in request.GET:
            try:
                limit = int(request.GET[u'limit'])
                query = query[:limit]
            except ValueError:
                pass

    uniquenames = query.values('uniquename')
    uniquenames = [uniqna['uniquename'] for uniqna in uniquenames]
    return HttpResponse(json.dumps(uniquenames),
                        content_type='application/json')


def get_all_children(feature):
    children = set()
    for child in feature.children:
        if child:
            children.add(child)
            children.update(get_all_children(child))
    return children


def api_features_children(request):

    request_data = request.GET
    request_dict = dict(request_data)
    features = []
    if 'features[]' in request_dict:
        for feat_uniquename in request_dict['features[]']:
            try:
                feat = Feature.objects.get(uniquename=feat_uniquename)
            except Feature.DoesNotExist:
                feat = None
            if feat is not None:
                features.append(feat)
    children = set()
    for feature in features:
        children.update(get_all_children(feature))

    return HttpResponse(json.dumps([c.uniquename for c in children]),
                        content_type='application/json')
