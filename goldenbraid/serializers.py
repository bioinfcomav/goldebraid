from rest_framework import serializers
from goldenbraid.models import Feature
from datetime import datetime
from django.contrib.auth import get_user_model


class FeatureSerializer(serializers.ModelSerializer):
    type = serializers.ReadOnlyField(source='type.name')

    def to_representation(self, instance):
        representation = super().to_representation(instance)
        representation['genbank_file'] = representation['genbank_file'].split("/")[-1]
        representation['timecreation'] = representation['timecreation'].split("T")[0]
        representation['gb_category'] = str(representation['gb_category']).split(" ")[0]

        return representation


    class Meta:
        model = Feature
        fields = ('uniquename', 'name', 'type', 'prefix', 'suffix', 'field_description', 'gb_category', 'field_owner', 'timecreation', 'genbank_file')
