# Copyright 2013 Diego Orzaez, Univ.Politecnica Valencia, Consejo Superior de
# Investigaciones Cientificas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from collections import OrderedDict
from operator import itemgetter

from django import forms
from django.forms.widgets import Select
from django.db.models import Q

from Bio.Seq import Seq

from goldenbraid.models import Feature
from goldenbraid.tags import (VECTOR_TYPE_NAME, ENZYME_IN_TYPE_NAME,
                              PROM_CAS12, TARGET_CAS12A, CLEAVAGE_SIGNAL,
                              TARGET_DICOT, TARGET_MONOCOT, TER_CRYSPER,
                              FUNGAL_PROM_5UTR, FUNGAL_CDS, FUNGAL_UTR3_TERM,
                              REVERSE_MARKER, FORWARD_MARKER,
                              FLANK_5UTR, FLANK_3UTR)
from goldenbraid.settings import (PARTS_TO_ASSEMBLE, UT_SUFFIX,
                                  UT_PREFIX, SITE_B, SITE_A, SITE_C,
                                  BIPARTITE_ALLOWED_PARTS,
                                  FUNGAL_BIPARTITE_ALLOWED_PARTS)
from goldenbraid.forms.feature import (features_to_choices,
                                       create_feature_validator)
from goldenbraid.utils import filter_feature_by_user_perms


def get_vector_choices(user):
    vectors = Feature.objects.filter(type__name=VECTOR_TYPE_NAME)
    vectors = filter_feature_by_user_perms(vectors, user)
    vector_choices = _vectors_to_choice(vectors)
    return vector_choices


def _vectors_to_choice(vectors):
    "it returns the given vectors but prepared to use as choices in a select"
    for_vectors = vectors.filter(prefix=UT_SUFFIX, suffix=UT_PREFIX)
    rev_vectors = vectors.filter(prefix=Seq(UT_PREFIX).reverse_complement(),
                                 suffix=Seq(UT_SUFFIX).reverse_complement())
    for_vector_choices = features_to_choices(for_vectors, blank_line=False)
    rev_vector_choices = features_to_choices(rev_vectors, blank_line=False)
    vector_choices = (('', ''),
                      ('Forward vectors', for_vector_choices),
                      ('Reverse vectors', rev_vector_choices))

    return vector_choices


def get_multipartite_form(multi_type, user):
    'It returns a form for the given multipartite'
    form_fields = OrderedDict()

    part_defs = PARTS_TO_ASSEMBLE[multi_type]
    print(multi_type)
    for parts in part_defs:
        features = Feature.objects.filter(type__name=parts[0],
                                          prefix=parts[1],
                                          suffix=parts[2])
        if parts[0] == TER_CRYSPER:
            if multi_type == 'gRNA_dicot' or multi_type == 'gRNA_monocot':
                features = Feature.objects.filter(uniquename="GB0645")
            elif multi_type == 'gRNA_monocot_reg' or multi_type == 'gRNA_dicot_reg':
                features = Feature.objects.filter(Q(type__name=TER_CRYSPER) & ~Q(uniquename="GB0645"))
        features = filter_feature_by_user_perms(features, user)
        choices = features_to_choices(features)
        name = parts[0]
        form_fields[name] = forms.CharField(max_length=100,
                                            widget=Select(choices=choices))

    # last we need to add the vector to the form
    vector_choices = get_vector_choices(user)
    form_fields[VECTOR_TYPE_NAME] = forms.CharField(max_length=100,
                                                    widget=Select(choices=vector_choices))

    form = type('MultiPartiteForm', (forms.BaseForm,),
                {'base_fields': form_fields})
    for field_name in form_fields.keys():
        setattr(form, 'clean_{0}'.format(field_name),
                create_feature_validator(field_name))
    return form

class Cas12aSingleTUForm(forms.Form):
    promoter = forms.CharField(max_length=100, widget=Select())
    target = forms.CharField(max_length=100, widget=Select())
    signal = forms.CharField(max_length=100, widget=Select())
    vector = forms.CharField(max_length=100, widget=Select())



def populate_cas12a_single_tu_form(form, user):
    promoter_choices = [("", "")]
    promoters = Feature.objects.filter(Q(type__name__iexact=PROM_CAS12)).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))
    for promoter in promoters:
        promoter_choices.append((promoter.uniquename, promoter.name))

    target_choices = [("", "")]
    cas12_targets = Feature.objects.filter(type__name__iexact=TARGET_CAS12A).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))
    for target in cas12_targets:
        target_choices.append((target.uniquename, target.name))

    signal_choices = [("", "")]
    signals = Feature.objects.filter(type__name__iexact=CLEAVAGE_SIGNAL).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))
    for signal in signals:
        signal_choices.append((signal.uniquename, signal.name))

    form.fields['promoter'].widget.choices = promoter_choices
    form.fields['target'].widget.choices = target_choices
    form.fields['signal'].widget.choices = signal_choices
    form.fields['vector'].widget.choices = get_vector_choices(user)

class FungalTUForm(forms.Form):
    promoter = forms.CharField(max_length=100, widget=Select())
    cds = forms.CharField(max_length=100, widget=Select())
    ter = forms.CharField(max_length=100, widget=Select())
    vector = forms.CharField(max_length=100, widget=Select())

def populate_fungal_tu_form(form, user):
    promoters = Feature.objects.filter(Q(type__name__iexact=FUNGAL_PROM_5UTR)).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))
    cdss = Feature.objects.filter(Q(type__name__iexact=FUNGAL_CDS)).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))
    ters = Feature.objects.filter(Q(type__name__iexact=FUNGAL_UTR3_TERM)).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))

    promoter_choices = [("", "")]
    for promoter in promoters:
        promoter_choices.append((promoter.uniquename, promoter.name))

    cds_choices = [("", ""), ("GB0211", "pHygromycin"), ("GB0053", "pYFP")]
    for cds in cdss:
        cds_choices.append((cds.uniquename, cds.name))

    ter_choices = [("", "")]
    for ter in ters:
        ter_choices.append((ter.uniquename, ter.name))

    form.fields['promoter'].widget.choices = promoter_choices
    form.fields['cds'].widget.choices = cds_choices
    form.fields['ter'].widget.choices = ter_choices
    form.fields['vector'].widget.choices = get_vector_choices(user)

class GeneDisruptionForm(forms.Form):
    negative_marker = forms.CharField(max_length=100, widget=Select())
    five_prime_flank = forms.CharField(max_length=100, widget=Select())
    positive_marker = forms.CharField(max_length=100, widget=Select())
    three_prime_flank = forms.CharField(max_length=100, widget=Select())
    vector = forms.CharField(max_length=100, widget=Select())

def populate_gene_disruption_form(form, user):
    reverse_markers = Feature.objects.filter(Q(type__name__iexact=REVERSE_MARKER)).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))
    five_flanks = Feature.objects.filter(Q(type__name__iexact=FLANK_5UTR)).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))
    forward_markers = Feature.objects.filter(Q(type__name__iexact=FORWARD_MARKER)).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))
    three_flanks = Feature.objects.filter(Q(type__name__iexact=FLANK_3UTR)).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))


    reverse_markers_choices = [("", "")]
    for reverse_marker in reverse_markers:
        reverse_markers_choices.append((reverse_marker.uniquename, reverse_marker.name))

    five_flank_choices = [("", "")]
    for five_flank in five_flanks:
        five_flank_choices.append((five_flank.uniquename, five_flank.name))

    forward_markers_choices = [("", "")]
    for forward_marker in forward_markers:
        forward_markers_choices.append((forward_marker.uniquename, forward_marker.name))

    three_flank_choices = [("", "")]
    for three_flank in three_flanks:
        three_flank_choices.append((three_flank.uniquename, three_flank.name))

    form.fields["negative_marker"].widget.choices = reverse_markers_choices
    form.fields["five_prime_flank"].widget.choices = five_flank_choices
    form.fields["positive_marker"].widget.choices = forward_markers_choices
    form.fields["three_prime_flank"].widget.choices = three_flank_choices
    form.fields["vector"].widget.choices = get_vector_choices(user)


class Cas12aMultiplexTUForm(forms.Form):
    promoter = forms.CharField(max_length=100, widget=Select())
    target = forms.CharField(max_length=100, widget=Select())
    vector = forms.CharField(max_length=100, widget=Select())


def populate_cas12a_multiplex_tu_form(form, user):
    promoter_choices = [("", "")]
    promoter = Feature.objects.get(uniquename="GB1443")
    promoter_choices.append((promoter.uniquename, promoter.name))

    target_choices = [("", "")]
    cas12_targets = Feature.objects.filter((Q(prefix__iexact="AGAT") & Q(suffix__iexact="CGCT"))).filter(Q(featureperm__owner__username=user) | Q(featureperm__is_public=True))
    for target in cas12_targets:
        target_choices.append((target.uniquename, target.name))


    form.fields['promoter'].widget.choices = promoter_choices
    form.fields['target'].widget.choices = target_choices
    form.fields['vector'].widget.choices = get_vector_choices(user)       


class MultipartiteFormFreeInitial(forms.Form):
    vector = forms.CharField(max_length=100, widget=Select(choices=[]))

    def clean_vector(self):
        return create_feature_validator('vector')(self)


def get_multipartite_free_form(feat_uniquenames):
    form_fields = OrderedDict()
    count = 0
    for feat_uniquename in feat_uniquenames:
        if count == 0:
            part_name = 'vector'
        else:
            part_name = 'part_{0}'.format(count)
        field = forms.CharField(required=True, initial=feat_uniquename)
        field.widget.attrs['readonly'] = True
        form_fields[part_name] = field
        count += 1

    form = type('MultiPartiteFreeValForm', (forms.BaseForm,),
                {'base_fields': form_fields})

    for field_name in form_fields.keys():
        setattr(form, 'clean_{0}'.format(field_name),
                create_feature_validator(field_name))
    return form


## Bipartite fungal ##
def fungal_get_part1_choice(user):
    _bi_parts = Feature.objects.filter(type__name__in=FUNGAL_BIPARTITE_ALLOWED_PARTS)
    print(_bi_parts)
    for part in _bi_parts:
        print(part.prefix, part.suffix, part.type.name)
    _parts = _bi_parts.filter(prefix=SITE_A, suffix=SITE_C)
    _parts = filter_feature_by_user_perms(_parts, user)

    parts_forw = _parts.filter(vector__prefix=SITE_B, vector__suffix=SITE_A)
    parts_rev = _parts.filter(vector__prefix=Seq(SITE_A).reverse_complement(),
                              vector__suffix=Seq(SITE_B).reverse_complement())
    part_forw_choices = features_to_choices(parts_forw, blank_line=False)
    part_rev_choices = features_to_choices(parts_rev, blank_line=False)
    part_choices = (('', ''),
                    ('Forward parts', part_forw_choices),
                    ('Reverse parts', part_rev_choices))
    return part_choices


class FungalBipartiteForm1(forms.Form):
    part_1 = forms.CharField(max_length=100, widget=Select(choices=[]))

    def clean_part_1(self):
        return create_feature_validator('part_1')(self)


class FungalBipartiteForm2(forms.Form):
    part_1 = forms.CharField(max_length=100)
    part_1.widget.attrs['readonly'] = True

    part_2 = forms.CharField(max_length=100, widget=Select(choices=[]))

    def clean_part_2(self):
        return create_feature_validator('part_2')(self)


class FungalBipartiteForm3(forms.Form):
    part_1 = forms.CharField(max_length=100)
    part_1.widget.attrs['readonly'] = True

    part_2 = forms.CharField(max_length=100)
    part_2.widget.attrs['readonly'] = True

    Vector = forms.CharField(max_length=100, widget=Select(choices=[]))

    def clean_part_1(self):
        return create_feature_validator('part_1')(self)

    def clean_part_2(self):
        return create_feature_validator('part_2')(self)

    def clean_Vector(self):
        return create_feature_validator('Vector')(self)


def fungal_get_part2_choices(part1_uniquename, user):
    part1 = Feature.objects.get(uniquename=part1_uniquename)
    part1_enzyme_out = part1.enzyme_out
    bi_parts = Feature.objects.filter(type__name__in=FUNGAL_BIPARTITE_ALLOWED_PARTS)
    parts = bi_parts.filter(prefix=SITE_C, suffix=SITE_B)
    parts = filter_feature_by_user_perms(parts, user)

    parts_forw = parts.filter(vector__prefix=SITE_B, vector__suffix=SITE_A)
    parts_rev = parts.filter(vector__prefix=Seq(SITE_A).reverse_complement(),
                             vector__suffix=Seq(SITE_B).reverse_complement())
    part_forw_choices = []
    for part in parts_forw:
        if part.enzyme_out == part1_enzyme_out:
            uniquename = part.uniquename #.encode('utf-8')
            if part.name:
                show = u'{0} - {1}'.format(uniquename, part.name)
            else:
                show = uniquename
            part_forw_choices.append((uniquename, show))

    part_rev_choices = []
    for part in parts_rev:
        if part.enzyme_out == part1_enzyme_out:
            uniquename = part.uniquename #.encode('utf-8')
            if part.name:
                show = u'{0} - {1}'.format(uniquename, part.name)
            else:
                show = uniquename
            part_rev_choices.append((uniquename, show))

    part_forw_choices = sorted(part_forw_choices, key=itemgetter(0))
    part_rev_choices = sorted(part_rev_choices, key=itemgetter(0))

    part_choices = (('', ''),
                    ('Forward parts', part_forw_choices),
                    ('Reverse parts', part_rev_choices))
    return part_choices


def fungal_get_bipart_vector_choices(part_uniquename, user):
    part = Feature.objects.get(uniquename=part_uniquename)
    part_enzyme_out = part.enzyme_out[0]

    vectors = Feature.objects.filter(type__name=VECTOR_TYPE_NAME)
    vectors = vectors.filter(featureprop__type__name=ENZYME_IN_TYPE_NAME,
                             featureprop__value=part_enzyme_out)
    vectors = filter_feature_by_user_perms(vectors, user)

    return _vectors_to_choice(vectors)


# # Bipartite ##
def get_part1_choice(user):
    _bi_parts = Feature.objects.filter(type__name__in=BIPARTITE_ALLOWED_PARTS)
    _parts = _bi_parts.filter(prefix=SITE_A, suffix=SITE_C)
    _parts = filter_feature_by_user_perms(_parts, user)

    parts_forw = _parts.filter(vector__prefix=SITE_B, vector__suffix=SITE_A)
    parts_rev = _parts.filter(vector__prefix=Seq(SITE_A).reverse_complement(),
                              vector__suffix=Seq(SITE_B).reverse_complement())
    part_forw_choices = features_to_choices(parts_forw, blank_line=False)
    part_rev_choices = features_to_choices(parts_rev, blank_line=False)
    part_choices = (('', ''),
                    ('Forward parts', part_forw_choices),
                    ('Reverse parts', part_rev_choices))
    return part_choices


class BipartiteForm1(forms.Form):
    part_1 = forms.CharField(max_length=100, widget=Select(choices=[]))

    def clean_part_1(self):
        return create_feature_validator('part_1')(self)


class BipartiteForm2(forms.Form):
    part_1 = forms.CharField(max_length=100)
    part_1.widget.attrs['readonly'] = True

    part_2 = forms.CharField(max_length=100, widget=Select(choices=[]))

    def clean_part_2(self):
        return create_feature_validator('part_2')(self)


class BipartiteForm3(forms.Form):
    part_1 = forms.CharField(max_length=100)
    part_1.widget.attrs['readonly'] = True

    part_2 = forms.CharField(max_length=100)
    part_2.widget.attrs['readonly'] = True

    Vector = forms.CharField(max_length=100, widget=Select(choices=[]))

    def clean_part_1(self):
        return create_feature_validator('part_1')(self)

    def clean_part_2(self):
        return create_feature_validator('part_2')(self)

    def clean_Vector(self):
        return create_feature_validator('Vector')(self)


def get_part2_choices(part1_uniquename, user):
    part1 = Feature.objects.get(uniquename=part1_uniquename)
    part1_enzyme_out = part1.enzyme_out
    bi_parts = Feature.objects.filter(type__name__in=BIPARTITE_ALLOWED_PARTS)
    parts = bi_parts.filter(prefix=SITE_C, suffix=SITE_B)
    parts = filter_feature_by_user_perms(parts, user)

    parts_forw = parts.filter(vector__prefix=SITE_B, vector__suffix=SITE_A)
    parts_rev = parts.filter(vector__prefix=Seq(SITE_A).reverse_complement(),
                             vector__suffix=Seq(SITE_B).reverse_complement())
    part_forw_choices = []
    for part in parts_forw:
        if part.enzyme_out == part1_enzyme_out:
            uniquename = part.uniquename #.encode('utf-8')
            if part.name:
                show = u'{0} - {1}'.format(uniquename, part.name)
            else:
                show = uniquename
            part_forw_choices.append((uniquename, show))

    part_rev_choices = []
    for part in parts_rev:
        if part.enzyme_out == part1_enzyme_out:
            uniquename = part.uniquename #.encode('utf-8')
            if part.name:
                show = u'{0} - {1}'.format(uniquename, part.name)
            else:
                show = uniquename
            part_rev_choices.append((uniquename, show))

    part_forw_choices = sorted(part_forw_choices, key=itemgetter(0))
    part_rev_choices = sorted(part_rev_choices, key=itemgetter(0))

    part_choices = (('', ''),
                    ('Forward parts', part_forw_choices),
                    ('Reverse parts', part_rev_choices))
    return part_choices


def get_bipart_vector_choices(part_uniquename, user):
    part = Feature.objects.get(uniquename=part_uniquename)
    part_enzyme_out = part.enzyme_out[0]

    vectors = Feature.objects.filter(type__name=VECTOR_TYPE_NAME)
    vectors = vectors.filter(featureprop__type__name=ENZYME_IN_TYPE_NAME,
                             featureprop__value=part_enzyme_out)
    vectors = filter_feature_by_user_perms(vectors, user)

    return _vectors_to_choice(vectors)
