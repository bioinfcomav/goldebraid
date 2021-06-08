from Bio.Seq import Seq
from dal import autocomplete
from django.db.models import Q
from goldenbraid.models import Feature
from goldenbraid.settings import (BIPARTITE_ALLOWED_PARTS, SITE_A, SITE_B, SITE_C, UT_PREFIX, UT_SUFFIX)
from goldenbraid.tags import VECTOR_TYPE_NAME, ENZYME_IN_TYPE_NAME
from goldenbraid.utils import filter_feature_by_user_perms
from django.utils.html import format_html



class BipartitePart1Autocomplete(autocomplete.Select2QuerySetView):

    def get_result_label(self, result):
        if result.vector.prefix == SITE_B and result.vector.suffix == SITE_A:
            direction = "Forward part"
        if result.vector.prefix == Seq(SITE_A).reverse_complement() and result.vector.suffix == Seq(SITE_B).reverse_complement():
            direction = "Reverse part"
        return format_html("{} - {} - {} - {}", direction, result.uniquename, result.name, result.vector)

    def get_result_value(self, result):
        return str(result.uniquename)

    def get_queryset(self):
        user = self.request.user
        if user.is_staff:
            queryset = Feature.objects.filter(Q(type__name__in=BIPARTITE_ALLOWED_PARTS)
                                              & Q(prefix=SITE_A) & Q(suffix=SITE_C))
        if user.is_authenticated:
             queryset = Feature.objects.filter(Q(featureperm__owner__username=user) |
                                               Q(featureperm__is_public=True))
             queryset = queryset.filter(Q(type__name__in=BIPARTITE_ALLOWED_PARTS)
                                        & Q(prefix=SITE_A) & Q(suffix=SITE_C))
        else:
            queryset = Feature.objects.filter(featureperm__is_public=True)
            queryset.filter(Q(type__name__in=BIPARTITE_ALLOWED_PARTS)
                            & Q(prefix=SITE_A) & Q(suffix=SITE_C))
        

        #queryset = filter_feature_by_user_perms(queryset, self.request.user)

        if self.q:
            queryset = queryset.filter(Q(name__icontains=self.q) | Q(uniquename__icontains=self.q) | Q(vector__name__icontains=self.q))
        if self.request.user.is_staff:
            queryset = [feature for feature in queryset if feature.direction is not None]
        queryset = sorted(queryset, key=lambda t: t.direction)
        return queryset


class BipartitePart2Autocomplete(autocomplete.Select2QuerySetView):


    def get_result_label(self, result):
        if result.vector.prefix == SITE_B and result.vector.suffix == SITE_A:
            direction = "Forward part"
        if result.vector.prefix == Seq(SITE_A).reverse_complement() and result.vector.suffix == Seq(SITE_B).reverse_complement():
            direction = "Reverse part"
        return format_html("{} - {} - {} - {}", direction, result.uniquename, result.name, result.vector)
    
    def get_result_value(self, result):
        return str(result.uniquename)

    def get_queryset(self):
        part1 = Feature.objects.get(uniquename=self.forwarded.get('part1', None))
        part1_enzyme_out = part1.enzyme_out
        queryset = Feature.objects.filter(Q(type__name__in=BIPARTITE_ALLOWED_PARTS)
                                          & Q(prefix=SITE_C) & Q(suffix=SITE_B))
        

        queryset = filter_feature_by_user_perms(queryset, self.request.user)
        
        selectable_uniquenames = [feature.uniquename for feature in queryset if feature.enzyme_out == part1_enzyme_out]
        queryset = queryset.filter(uniquename__in=selectable_uniquenames)
        

        

        if self.q:
            queryset = queryset.filter(Q(name__icontains=self.q) | Q(uniquename__icontains=self.q))
        if self.request.user.is_staff:
            queryset = [feature for feature in queryset if feature.direction is not None]
       
        queryset = sorted(queryset, key=lambda t: t.direction)
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
        queryset = Feature.objects.filter(type__name=VECTOR_TYPE_NAME)
        
        print("potential parts")
        print(len(queryset))

        queryset = filter_feature_by_user_perms(queryset, self.request.user)
        print("user perms")
        print(len(queryset))
        
        
        queryset = queryset.filter(featureprop__type__name=ENZYME_IN_TYPE_NAME,
                                 featureprop__value=part1_enzyme_out)
        print("same enzyme")
        print(len(queryset))

        if self.q:
            queryset = queryset.filter(Q(name__icontains=self.q) | Q(uniquename__icontains=self.q))
        return queryset  


# def get_bipart_vector_choices(part_uniquename, user):
#     part = Feature.objects.get(uniquename=part_uniquename)
#     part_enzyme_out = part.enzyme_out[0]

#     vectors = Feature.objects.filter(type__name=VECTOR_TYPE_NAME)
#     vectors = vectors.filter(featureprop__type__name=ENZYME_IN_TYPE_NAME,
#                              featureprop__value=part_enzyme_out)
#     vectors = filter_feature_by_user_perms(vectors, user)

#     return _vectors_to_choice(vectors)


# def _vectors_to_choice(vectors):
#     "it returns the given vectors but prepared to use as choices in a select"
#     for_vectors = vectors.filter(prefix=UT_SUFFIX, suffix=UT_PREFIX)
#     rev_vectors = vectors.filter(prefix=Seq(UT_PREFIX).reverse_complement(),
#                                  suffix=Seq(UT_SUFFIX).reverse_complement())
#     for_vector_choices = features_to_choices(for_vectors, blank_line=False)
#     rev_vector_choices = features_to_choices(rev_vectors, blank_line=False)
#     vector_choices = (('', ''),
#                       ('Forward vectors', for_vector_choices),
#                       ('Reverse vectors', rev_vector_choices))

#     return vector_choices