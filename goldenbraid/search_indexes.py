from django.db import models
from django.db.models.fields import CharField
from haystack import indexes
from goldenbraid.models import Feature


class FeatureIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    uniquename = indexes.CharField(model_attr='uniquename')
    name = indexes.CharField(model_attr='name')
    prefix = indexes.CharField(model_attr='prefix')
    suffix = indexes.CharField(model_attr='suffix')
    type = indexes.CharField(model_attr="type")
    vector = indexes.CharField(model_attr="vector", null=True)
    vector_prefix = indexes.CharField(model_attr="vector", null=True)
    vector_suffix = indexes.CharField(model_attr="vector", null=True)
    enzyme_out = indexes.MultiValueField(null=True)
    enzyme_in = indexes.CharField(null=True)
    owner = indexes.CharField(null=True)
    is_public = indexes.BooleanField(null=True)
    content_auto = indexes.EdgeNgramField(use_template=True)



    def get_model(self):
        return Feature

    def prepare_type(self, obj):
        return obj.type.name
    
    def prepare_vector(self, obj):
        try:
            return obj.vector.name
        except:
            return None
    
    def prepare_vector_prefix(self, obj):
        try:
            return obj.vector.prefix
        except:
            return None

    def prepare_vector_suffix(self, obj):
        try:
            return obj.vector.suffix
        except:
            return None

    def prepare_enzyme_out(self, obj):
        if obj.enzyme_out:
            return obj.enzyme_out
        else:
            return None

    def prepare_enzyme_in(self, obj):
        if obj.enzyme_in:
            return obj.enzyme_in
        else:
            return None

    def prepare_owner(self, obj):
        try:
            return obj.owner
        except:
            return None

    def prepare_is_public(self, obj):
        try:
            return obj.is_public
        except:
            return None

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        return Feature.objects.all()