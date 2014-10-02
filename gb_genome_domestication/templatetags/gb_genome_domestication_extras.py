
from django import template

from gb_genome_domestication.models import Feature

register = template.Library()


def rec_sites(name):
    try:
        feat = Feature.objects.get(uniquename=name)
        num_rec_sites = str(feat.num_rec_sites)
    except Feature.DoesNotExist:
        num_rec_sites = 'Unknown'
    return str(num_rec_sites)

register.filter('rec_sites', rec_sites)


def species(name):
    try:
        feat = Feature.objects.get(uniquename=name)
        species_ = feat.species
    except Feature.DoesNotExist:
        species_ = 'Unknown'
    return species_.replace('_', ' ').capitalize()

register.filter('species', species)
