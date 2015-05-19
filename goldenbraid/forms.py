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
from django.core.exceptions import ValidationError
from django.forms.widgets import Select
from django.forms.utils import ErrorDict
from django.db.models import Q

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from goldenbraid.utils import has_rec_sites
from goldenbraid.models import Cvterm, Feature
from goldenbraid.tags import (VECTOR_TYPE_NAME, ENZYME_IN_TYPE_NAME,
                              MODULE_TYPE_NAME, OTHER_TYPE_NAME, TU_TYPE_NAME,
                              TER_CRYSPER, PROM_MONOCOT, PROM_DICOT,
                              TARGET_DICOT, TARGET_MONOCOT)
from goldenbraid.settings import (PARTS_TO_ASSEMBLE, UT_SUFFIX,
                                  UT_PREFIX, SITE_B, SITE_A, SITE_C,
                                  BIPARTITE_ALLOWED_PARTS, CATEGORIES,
                                  SEARCH_MENU_TYPE_CHOICES,
                                  MINIMUN_PCR_LENGTH, CRYSPER_CATEGORIES,
                                  CRYSPER_TARGETS_TO_DOMESTICATE)


def get_vector_choices(user):
    vectors = Feature.objects.filter(type__name=VECTOR_TYPE_NAME)
    if user.is_authenticated():
        vectors = vectors.filter(Q(featureperm__owner__username=user) |
                                 Q(featureperm__is_public=True))
    else:
        vectors = vectors.filter(featureperm__is_public=True)

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


def features_to_choices(features, blank_line=True):
    choices = [('', '')] if blank_line else []

    for feat in features:
        uniquename = feat.uniquename.encode('utf-8')
        if feat.name:
            show = u'{0} - {1}'.format(feat.uniquename, feat.name)
        else:
            show = uniquename
        choices.append((uniquename, show))
    choices = sorted(choices, key=itemgetter(0))
    return choices


def _prepare_feature_kind():
    'It prepares the feature kind select choices to put in the type widget'
    kinds = [cat[0] for cat in CATEGORIES.values()]
    kinds.extend([TU_TYPE_NAME, OTHER_TYPE_NAME, MODULE_TYPE_NAME,
                  PROM_DICOT, PROM_MONOCOT, TER_CRYSPER])

    feature_kinds = [(kind, kind.replace('_', ' ')) for kind in kinds]

    feature_kinds.insert(0, ('', ''))  # no kind
    return feature_kinds


def _prepare_feature_kind_old():
    'It prepares the feature kind select choices to put in the type widget'
    if SEARCH_MENU_TYPE_CHOICES:
        kinds = SEARCH_MENU_TYPE_CHOICES
    else:
        kinds = Feature.objects.distinct('type').values('type__name')
        kinds = [kind['type__name'] for kind in kinds]
    if VECTOR_TYPE_NAME in kinds:
        kinds.pop(kinds.index(VECTOR_TYPE_NAME))

    feature_kinds = [(kind, kind.replace('_', ' ')) for kind in kinds]

    feature_kinds.insert(0, ('', ''))  # no kind
    return feature_kinds


def get_all_vectors_as_choices(user):
    vectors = Feature.objects.filter(type__name=VECTOR_TYPE_NAME)
    if user.is_authenticated():
        vectors = vectors.filter(Q(featureperm__owner__username=user) |
                             Q(featureperm__is_public=True))
    else:
        vectors = vectors.filter(featureperm__is_public=True)
    return features_to_choices(vectors, blank_line=True)


class VectorForm(forms.Form):
    'Form to add vectors to db'
    name = forms.CharField(max_length=255, required=False)
    description = forms.CharField(max_length=255)
    reference = forms.CharField(max_length=255, required=False)
    enzyme_in = forms.CharField(max_length=255)
    enzyme_out = forms.CharField(max_length=255)
    resistance = forms.CharField(max_length=255)
    gbfile_label = 'Select a GenBank-formatted local file on your computer'
    gbfile = forms.FileField(label=gbfile_label, required=True)


class FeatureForm(forms.Form):
    'Form to add features to db'
    name = forms.CharField(max_length=255, required=False)
    description = forms.CharField(max_length=255, required=False)
    reference = forms.CharField(max_length=255, required=False)

    type_choices = _prepare_feature_kind()
    type = forms.CharField(max_length=100, widget=Select(choices=type_choices))

    vector = forms.CharField(max_length=100, widget=Select(choices=[]))

    gbfile_label = 'Select a GenBank-formatted local file on your computer'
    gbfile = forms.FileField(label=gbfile_label, required=True)

    def clean_type(self):
        'It validates the type field'
        type_str = self.cleaned_data['type']
        try:
            Cvterm.objects.get(name=type_str)
        except Cvterm.DoesNotExist:
            raise ValidationError('This type does not exist in the database')
        return type_str

    def clean_vector(self):
        '''It validates the vector.

        If the feature is a vector it does not validate anything
        if feature is not a vector if validates that the vector
        is in the database'''
        vector = self.cleaned_data['vector']
        error_in_type = self.errors.get('type', False)
        if error_in_type:
            return vector
        type_str = self.cleaned_data['type']

        if type_str != VECTOR_TYPE_NAME:
            try:
                vector_type = Cvterm.objects.get(name=VECTOR_TYPE_NAME)
                Feature.objects.get(uniquename=vector,
                                              type=vector_type)
            except Feature.DoesNotExist:
                raise ValidationError('The given vector does not exist')

        else:
            if vector:
                raise ValidationError('A vector does not have a vector')

        return vector


def create_feature_validator(field_name):

    def validator(self):
        uniquename_str = self.cleaned_data[field_name]
        try:
            Feature.objects.get(uniquename=uniquename_str)
        except Feature.DoesNotExist:
            raise ValidationError('This feature does not exist in the database')
        return uniquename_str

    return validator


def get_multipartite_form(multi_type, user):
    'It returns a form for the given multipartite'
    form_fields = OrderedDict()

    part_defs = PARTS_TO_ASSEMBLE[multi_type]
    for parts in part_defs:
        features = Feature.objects.filter(type__name=parts[0],
                                          prefix=parts[1], suffix=parts[2])
        if user.is_staff:
            pass
        else:
            if user.is_authenticated():
                features = features.filter(Q(featureperm__owner__username=user) |
                                           Q(featureperm__is_public=True))
            else:
                features = features.filter(featureperm__is_public=True)

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


### Bipartite ######
def get_part1_choice(user):
    _bi_parts = Feature.objects.filter(type__name__in=BIPARTITE_ALLOWED_PARTS)
    _parts = _bi_parts.filter(prefix=SITE_A, suffix=SITE_C)
    if user.is_authenticated():
        _parts = _parts.filter(Q(featureperm__owner__username=user) |
                               Q(featureperm__is_public=True))

    else:
        _parts = _parts.filter(featureperm__is_public=True)
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
    if user.is_authenticated():
        parts = parts.filter(Q(featureperm__owner__username=user) |
                             Q(featureperm__is_public=True))
    else:
        parts = parts.filter(featureperm__is_public=True)

    parts_forw = parts.filter(vector__prefix=SITE_B, vector__suffix=SITE_A)
    parts_rev = parts.filter(vector__prefix=Seq(SITE_A).reverse_complement(),
                             vector__suffix=Seq(SITE_B).reverse_complement())
    part_forw_choices = []
    for part in parts_forw:
        if part.enzyme_out == part1_enzyme_out:
            uniquename = part.uniquename.encode('utf-8')
            if part.name:
                show = u'{0} - {1}'.format(uniquename, part.name)
            else:
                show = uniquename
            part_forw_choices.append((uniquename, show))

    part_rev_choices = []
    for part in parts_rev:
        if part.enzyme_out == part1_enzyme_out:
            uniquename = part.uniquename.encode('utf-8')
            if part.name:
                show = u'{0} - {1}'.format(uniquename, part.name)
            else:
                show = uniquename
            part_rev_choices.append((uniquename, show))

    part_forw_choices = sorted(part_forw_choices, key=itemgetter(0))
    part_rev_choices = sorted(part_rev_choices, key=itemgetter(0))

    part_choices = (('', ''), ('Forward parts', part_forw_choices),
                    ('Reverse parts', part_rev_choices))
    return part_choices


def get_bipart_vector_choices(part_uniquename, user):
    part = Feature.objects.get(uniquename=part_uniquename)
    part_enzyme_out = part.enzyme_out[0]

    vectors = Feature.objects.filter(type__name=VECTOR_TYPE_NAME)
    vectors = vectors.filter(featureprop__type__name=ENZYME_IN_TYPE_NAME,
                             featureprop__value=part_enzyme_out)
    if user.is_authenticated():
        vectors = vectors.filter(Q(featureperm__owner__username=user) |
                                 Q(featureperm__is_public=True))
    else:
        vectors = vectors.filter(featureperm__is_public=True)

    return _vectors_to_choice(vectors)


# Domesticator #
class DomesticationForm(forms.Form):
    super_choice = [('', '')]
    choices = []
    for category_name in CATEGORIES.keys():
        choices.append((category_name, category_name))
    super_choice.append(('GB parts', choices))
    # Add crysper Proms and TERM
    choices = []
    for category_name in CRYSPER_CATEGORIES.keys():
        if category_name in CRYSPER_TARGETS_TO_DOMESTICATE:
            continue
        choices.append((category_name, category_name))
    super_choice.append(('Crysper parts', choices))

    intron_label = 'The secuence has introns in lowercase'
    with_intron = forms.BooleanField(label=intron_label, required=False)
    category = forms.CharField(max_length=100,
                               label='Choose a category to domesticate to',
                               widget=Select(choices=super_choice), required=False)
    seq = forms.FileField(max_length=100,
                          label='Add a genbank or a fast file', required=False)
    residues = forms.CharField(widget=forms.HiddenInput(), required=False)
    prefix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)
    suffix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)

    def clean_category(self):
        category_name = self.cleaned_data['category']
        if not category_name:
            del self.cleaned_data['category']
            return
        if category_name not in CATEGORIES.keys():
            raise ValidationError('You must choose a valid category')
        return category_name

    def _seq_validation(self, seq):
        with_intron = self.cleaned_data['with_intron']
        if not _seq_is_dna(seq.seq):
            msg = 'The given file contains seqs with not allowed nucleotides'
            msg += ' ATGC'
            raise ValidationError(msg)
        if len(seq) < MINIMUN_PCR_LENGTH + 20:
            msg = 'Given seq must be at least 70 base pairs'
            raise ValidationError(msg)
        if self._data_in(self.cleaned_data, 'category'):
            category = self.cleaned_data['category']
            if category in ('CDS (B3-B4-B5)', 'SP (B3)', 'NTAG (B2)',
                            'CDS (B3-B4)'):
                if not _seq_has_codon_start(seq.seq, with_intron):
                    msg = 'The provided seq must start with start codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in ('CDS (B3-B4-B5)', 'CDS (B4-B5)', 'CTAG (B5)'):
                if not _seq_has_codon_end(seq.seq, with_intron):
                    msg = 'The provided seq must end with a end codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in ('CDS (B3-B4-B5)', 'SP (B3)', 'NTAG (B2)',
                            'CDS (B3-B4)', 'CDS (B4-B5)', 'CTAG (B5)'):
                if not _is_seq_3_multiple(seq.seq, with_intron):
                    msg = 'The provided seq must be multiple of three in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in ('goi (B2-B3)'):
                if len(seq) > 500:
                    msg = 'The provided seq must have less than 500 nucleotides in'
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in ('CDS (B3-B4)', 'NTAG (B2)', 'SP (B3)'):
                if _seq_has_codon_end(seq.seq, with_intron):
                    msg = 'The provided seq must not end with a stop codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
        return seq

    def clean_seq(self):
        inmemoryfile = self.cleaned_data['seq']
        if not inmemoryfile:
            del self.cleaned_data['seq']
            return inmemoryfile
        content = inmemoryfile.chunks().next()
        self.cleaned_data['seq'].seek(0)
        if content.startswith('LOCUS') or content.startswith('>'):
            pass
        else:
            msg = 'The given file must be a fasta or a genbank file'
            raise ValidationError(msg)
        format_ = 'fasta' if content.startswith('>') else 'genbank'
        seq = SeqIO.read(self.cleaned_data['seq'], format_)
        seq = self._seq_validation(seq)
        return seq

    def clean_residues(self):
        residues = self.cleaned_data['residues']
        if not residues:
            del self.cleaned_data['residues']
            return residues
        seq = SeqRecord(seq=Seq(residues))
        seq = self._seq_validation(seq)
        return seq

    def _clean_customtags(self, kind):
        tag = self.cleaned_data[kind]
        if not tag:
            del self.cleaned_data[kind]
            return tag
        if len(tag) != 4:
            raise ValidationError('{0} tag must be of length 4'.format(kind))
        if not _seq_is_dna(tag):
            msg = 'The given tag seqs with not allowed nucleotides: ATGC'
            raise ValidationError(msg)
        return tag

    def clean_suffix(self):
        return self._clean_customtags('suffix')

    def clean_prefix(self):
        return self._clean_customtags('prefix')

    def full_clean(self):
        """
        Cleans all of self.data and populates self._errors and
        self.cleaned_data.
        """
        self._errors = ErrorDict()
        if not self.is_bound:  # Stop further processing.
            return
        self.cleaned_data = {}
        # If the form is permitted to be empty, and none of the form data has
        # changed from the initial data, short circuit any validation.
        if self.empty_permitted and not self.has_changed():
            return
        self._clean_fields()
        self._clean_form()
        self._post_clean()
        # custom validations
        self._seqform_validation()
        self._prefix_sufix_category_validation()
        if self._errors:
            del self.cleaned_data

    @staticmethod
    def _data_in(dictionary, key):
        return True if key in dictionary and dictionary[key] else False

    def _seqform_validation(self):
        cleaned_data = self.cleaned_data
        if 'seq' not in self._errors and 'residues' not in self._errors:
            try:
                if (not self._data_in(cleaned_data, 'seq') and
                        not self._data_in(cleaned_data, 'residues')):
                    raise ValidationError('Fasta or genbank File Required')

            except ValidationError, e:
                self._errors['seq'] = self.error_class(e.messages)
                for name in ('seq', 'residues'):
                    if name in self.cleaned_data:
                        del self.cleaned_data[name]
                return

            try:
                if (self._data_in(cleaned_data, 'seq') and
                        self._data_in(cleaned_data, 'residues')):
                    raise ValidationError('Form can not accept File and residues')

            except ValidationError, e:
                self._errors['seq'] = self.error_class(e.messages)
                for name in ('seq', 'residues'):
                    if name in self.cleaned_data:
                        del self.cleaned_data[name]
                return

    def _prefix_sufix_category_validation(self):
        cleaned_data = self.cleaned_data
        try:
            if (not self._data_in(cleaned_data, 'category') and
                not self._data_in(cleaned_data, 'suffix') and
                    not self._data_in(cleaned_data, 'prefix')):
                msg = 'At least we need category or prefix/suffix pair'
                raise ValidationError(msg)
        except ValidationError, e:
            name = 'category'
            self._errors[name] = self.error_class(e.messages)
            if name in self.cleaned_data:
                del self.cleaned_data[name]
            return

        try:
            if (self._data_in(cleaned_data, 'category') and
                (self._data_in(cleaned_data, 'suffix') or
                 self._data_in(cleaned_data, 'prefix'))):
                msg = 'Can not use category and prefix/suffix simoultaneously'
                raise ValidationError(msg)
        except ValidationError, e:
            name = 'category'
            self._errors[name] = self.error_class(e.messages)
            if name in self.cleaned_data:
                del self.cleaned_data[name]
            return

        try:
            if (not self._data_in(cleaned_data, 'category') and
                (not self._data_in(cleaned_data, 'suffix') or
                 not self._data_in(cleaned_data, 'prefix'))):
                msg = 'You must provide prefix and suffix together'
                raise ValidationError(msg)
        except ValidationError, e:
            if not self._data_in(cleaned_data, 'suffix'):
                name = 'suffix'
            else:
                name = 'prefix'
            self._errors[name] = self.error_class(e.messages)
            if name in self.cleaned_data:
                del self.cleaned_data[name]


def _seq_is_dna(string):
    len_sum = sum([string.count(l.upper()) + string.count(l.lower()) for l
                                                                  in ('ATCG')])
    return False if len_sum != len(string) else True


# this function is repeated in domestication module. But it can not be imported
# because cross imports. I copied here as it is very simple
def get_upper_nucls(seq):
    return ''.join([nucl for nucl in seq if nucl.isupper()])


def _seq_has_codon_start(seq, with_intron):
    if with_intron:
        seq = get_upper_nucls(seq)
    start = str(seq[:3].upper())
    return True if start == 'ATG' else False


def _seq_has_codon_end(seq, with_intron):
    if with_intron:
        seq = get_upper_nucls(seq)
    end = str(seq[-3:].upper())
    return True if end in ('TAG', 'TAA', 'TGA') else False


def _is_seq_3_multiple(seq, with_intron):
    if with_intron:
        seq = get_upper_nucls(seq)
    return True if divmod(len(seq), 3)[1] == 0 else False


class DomesticationAddForm(DomesticationForm):
    name = forms.CharField(max_length=255, required=False)
    description = forms.CharField(max_length=255, required=False)
    reference = forms.CharField(max_length=255, required=False)


# feature_management
class FeatureManagementForm(forms.Form):
    feature = forms.CharField(max_length=30, widget=forms.HiddenInput())
    action = forms.CharField(max_length=30, widget=forms.HiddenInput())

    def clean_action(self):
        action = self.cleaned_data['action']
        if action in ('delete', 'make_public', 'make_private'):
            return action
        raise ValidationError('action must be delete or make_public')

    def clean_feature(self):
        uniquename = self.cleaned_data['feature']
        return create_feature_validator('feature')(self)


class DomesticationCrysperForm(DomesticationForm):
    choices = [('', '')]
    for category_name in CRYSPER_TARGETS_TO_DOMESTICATE:
        choices.append((category_name, category_name))
    category = forms.CharField(max_length=100,
                               label='Choose a category to domesticate to',
                               widget=Select(choices=choices), required=False)
    seq = forms.CharField(max_length=20, label='Add your sequence',
                          required=True)
    prefix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)
    suffix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)

    def clean_seq(self):
        seq = self.cleaned_data['seq']
        if len(seq) != 20:
            raise ValidationError('Crysper seq len must be 20 nucleotides')

        if has_rec_sites(seq):
            msg = 'This secuence can not be domesticated.'
            msg += 'It has internal restriction sites'
            raise ValidationError(msg)
        return seq

    def clean_residues(self):
        pass

    def clean_category(self):
        category_name = self.cleaned_data['category']
        if not category_name:
            del self.cleaned_data['category']
            return
        if category_name not in CRYSPER_CATEGORIES.keys():
            raise ValidationError('You must choose a valid category')
        return category_name

    def _seqform_validation(self):
        cleaned_data = self.cleaned_data
        if 'seq' not in cleaned_data:
            return
        seq = cleaned_data['seq']
        category = cleaned_data['category']
        if category == TARGET_DICOT and str(seq[0]).upper() != 'G':
            msg = 'First nucleotide must be G for target dicot category'
            self._errors['seq'] = self.error_class([msg])
        if category == TARGET_MONOCOT and str(seq[0]).upper() != 'A':
            msg = 'First nucleotide must be A for target monocot category'
            self._errors['seq'] = self.error_class([msg])


