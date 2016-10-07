'''
Created on 2015 api. 1

@author: peio
'''
import json
from tempfile import NamedTemporaryFile

from django.http.response import HttpResponse, Http404, HttpResponseForbidden
from django.db.models import Q
from Bio.SeqIO import read
from goldenbraid.models import Feature, ExperimentPropExcel, ExperimentKeyword
from goldenbraid.excel import draw_combined_graph
from goldenbraid.sbol import convert_to_sbol


def feature_sbol(request, uniquename):
    user = request.user
    print(uniquename)
    query = Feature.objects.filter(uniquename=uniquename)
    print(query)
    if not user.is_staff:
        query = query.filter(Q(featureperm__is_public=True) |
                             Q(featureperm__owner=user))
    print(query)
    if query:
        query[0].genbank
        seq = read(query[0].genbank_fileo, 'gb')
        response = HttpResponse(convert_to_sbol(seq), content_type='xml/plain')
        filename = seq.name + '.xml'
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}"'.format(filename)
        return response

    else:
        return Http404


def feature_uniquenames(request):
    user = request.user
    query = Feature.objects.all()
    if not user.is_staff:
        query = query.filter(Q(featureperm__is_public=True) |
                             Q(featureperm__owner=user))
    if request.method == 'GET':
        if u'term' in request.GET:
            term = request.GET['term']
            query = query.filter(uniquename__icontains=term)
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


def _get_all_children(feature):
    children = set()
    for child in feature.children:
        if child:
            children.add(child)
            children.update(_get_all_children(child))
    return children


def _feature_children(request):
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
        children.update(_get_all_children(feature))

    return children, features


def features_children(request):
    children = _feature_children(request)[0]
    return HttpResponse(json.dumps([c.uniquename for c in children]),
                        content_type='application/json')


def features_key_elements(request):
    children, features = _feature_children(request)
    all_feats = list(children) + features
    return HttpResponse(json.dumps([c.uniquename for c in all_feats]),
                        content_type='application/json')


def excel_image(request, excel_id):
    try:
        exp_excel = ExperimentPropExcel.objects.get(experiment_prop_excel_id=excel_id)
    except ExperimentPropExcel.DoesNotExist:
        return Http404
    image_content, content_type = exp_excel.drawed_image
    return HttpResponse(image_content, content_type=content_type)


def experiment_keywords(request):
    query = ExperimentKeyword.objects.all()
    if request.method == 'GET':
        if u'term' in request.GET:
            term = request.GET['term']
            query = query.filter(keyword__icontains=term)
        if u'limit' in request.GET:
            try:
                limit = int(request.GET[u'limit'])
                query = query[:limit]
            except ValueError:
                pass
    keywords = list({exp_keywords.keyword for exp_keywords in query})
    return HttpResponse(json.dumps(keywords), content_type='application/json')


def combined_excel_image(request, uniquename, exp_type):
    feat = Feature.objects.get(uniquename=uniquename)
    excel_datas = feat.combined_experiment_excel_data(exp_type)
    out_fhand = NamedTemporaryFile(suffix='.svg')
    draw_combined_graph(excel_datas, out_fhand)
    image_content = open(out_fhand.name).read()
    content_type = 'image/svg+xml'
    return HttpResponse(image_content, content_type=content_type)
