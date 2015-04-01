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
