from django.shortcuts import render
from rest_framework import viewsets, status,  generics
from django_filters import rest_framework as filters
from goldenbraid.models import Feature, FeaturePerm
from django.contrib.auth.models import User
from goldenbraid.serializers import FeatureSerializer
from django.db.models import Q
from django_filters.rest_framework import DjangoFilterBackend
from goldenbraid.forms.feature import  SPECIAL_SEARCH_CATEGORIES


def search_feature_index(request):
    request_data = request.GET
    if "category" not in request_data:
        category = "all"
        prefix = None
        suffix = None
    else:
        category , prefix, suffix = request_data["category"].split(",")
    print(category)
    return render(request, 'search_feature_index.html', context={"category": category,
                                                                 "prefix": prefix,
                                                                 "suffix": suffix})

class FeatureFilter(filters.FilterSet):
    type__name = filters.CharFilter(field_name="type__name", lookup_expr='iexact')
    prefix = filters.CharFilter(field_name="prefix", lookup_expr='iexact')
    suffix = filters.CharFilter(field_name="suffix", lookup_expr='iexact')


    class Meta:
        model = Feature
        fields = ['type__name', 'prefix', 'suffix']
        together = ['prefix', 'suffix']


def check_validity(queryset):
    kmers = []
    for feature in queryset:
        if feature.prefix not in kmers:
            kmers.append(feature.prefix)
        if feature.suffix not in kmers:
            kmers.append(feature.suffix)
    print(kmers)
    



class FeatureViewSet(viewsets.ModelViewSet):
    serializer_class = FeatureSerializer
    http_method_names = ['get', 'head']
    #filter_backends = [filters.SearchFilter]
    # filter_backends = (DjangoFilterBackend,)
    # filter_class = FeatureFilter
    # filter_set = ['type__name', 'prefix', 'suffix']
    # search_filters = ['type__name', 'prefix', 'suffix']

    def get_queryset(self):
        prefix = self.request.query_params.get('prefix', None)
        suffix = self.request.query_params.get('suffix', None)
        category = self.request.query_params.get('type__name', None)
        only_owned = self.request.query_params.get('onlymyparts', None)
        user = self.request.user
        if only_owned == "True":
            queryset = Feature.objects.filter(featureperm__owner__username=user)
        elif user.is_staff:
            queryset = Feature.objects.all()
        else: 
            queryset = Feature.objects.filter(Q(featureperm__is_public=True) |
                                              Q(featureperm__owner__username=user))
        if category == "all":
            return queryset

        if prefix != "None" and suffix != "None":
            queryset = queryset.filter(Q(prefix=prefix) &
                                       Q(suffix=suffix))
        elif category in SPECIAL_SEARCH_CATEGORIES:
            queryset = queryset.filter(type__name=category)


        check_validity(queryset)
        
        return queryset