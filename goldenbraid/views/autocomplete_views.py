from Bio.Seq import Seq
from dal import autocomplete
from django.db.models import Q
from goldenbraid.models import Feature
from goldenbraid.settings import (BIPARTITE_ALLOWED_PARTS, SITE_A, SITE_B, SITE_C, UT_PREFIX, UT_SUFFIX)
from goldenbraid.tags import VECTOR_TYPE_NAME
from goldenbraid.utils import filter_feature_by_user_indexed
from django.utils.html import format_html
from haystack.query import SearchQuerySet



class BipartitePart1Autocomplete(autocomplete.Select2QuerySetView):

    def get_result_label(self, result):
        if result.vector_prefix == SITE_B and result.vector_suffix == SITE_A:
            direction = "Forward part"
        if result.vector_prefix == Seq(SITE_A).reverse_complement() and result.vector_suffix == Seq(SITE_B).reverse_complement():
            direction = "Reverse part"
        return format_html("{} - {} - {} - {}", direction, result.uniquename, result.name, result.vector)

    def get_result_value(self, result):
        return str(result.uniquename)

    def get_queryset(self):
        queryset = SearchQuerySet().models(Feature).filter(Q(type__name__in=BIPARTITE_ALLOWED_PARTS)
                                                           & Q(prefix=SITE_A) & Q(suffix=SITE_C))
        queryset = filter_feature_by_user_indexed(queryset, self.request.user)
        
        if self.q:
            queryset = queryset.filter(Q(name__contains=self.q) | Q(uniquename__contains=self.q) | Q(vector__name__contains=self.q))
        return queryset


class BipartitePart2Autocomplete(autocomplete.Select2QuerySetView):


    def get_result_label(self, result):
        if result.vector_prefix == SITE_B and result.vector_suffix == SITE_A:
            direction = "Forward part"
        if result.vector_prefix == Seq(SITE_A).reverse_complement() and result.vector_suffix == Seq(SITE_B).reverse_complement():
            direction = "Reverse part"
        return format_html("{} - {} - {} - {}", direction, result.uniquename, result.name, result.vector)
    
    def get_result_value(self, result):
        return str(result.uniquename)

    def get_queryset(self):
        part1 = Feature.objects.get(uniquename=self.forwarded.get('part1', None))
        part1_enzyme_out = part1.enzyme_out[0]
        queryset = SearchQuerySet().models(Feature).filter(Q(type__name__in=BIPARTITE_ALLOWED_PARTS)
                                                           & Q(prefix=SITE_C) & Q(suffix=SITE_B) & Q(enzyme_out=part1_enzyme_out))
        queryset = filter_feature_by_user_indexed(queryset, self.request.user)
        
        if self.q:
            queryset = queryset.filter(Q(name__contains=self.q) | Q(uniquename__contains=self.q))
        return queryset


class BipartitePart3Autocomplete(autocomplete.Select2QuerySetView):


    def get_result_label(self, result, direction = None):
        
        if result.prefix == UT_SUFFIX and result.suffix == UT_PREFIX:
            direction = "Forward part"
        if result.prefix == Seq(UT_PREFIX).reverse_complement() and result.suffix == Seq(UT_SUFFIX).reverse_complement():
            direction = "Reverse part"
        if direction is not None:
            return format_html("{} - {} - {}", direction, result.uniquename, result.name)

    def get_result_value(self, result):
        return str(result.uniquename)

    def get_queryset(self):
        part1 = Feature.objects.get(uniquename=self.forwarded.get('part1', None))
        part1_enzyme_out = part1.enzyme_out[0]
        queryset = SearchQuerySet().models(Feature).filter(Q(type__name=VECTOR_TYPE_NAME) & Q(enzyme_in=part1_enzyme_out))
        queryset = filter_feature_by_user_indexed(queryset, self.request.user)
        

        if self.q:
            queryset = queryset.filter(Q(name__contains=self.q) | Q(uniquename__contains=self.q))
        return queryset  