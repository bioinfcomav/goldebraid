from django import forms
from django.forms.widgets import Select
from django.core.exceptions import ValidationError
from django.forms.utils import ErrorDict

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from goldenbraid.settings import (CATEGORIES, MINIMUN_PCR_LENGTH,
                                  CRYSPER_TARGETS_TO_DOMESTICATE,
                                  CRYSPER_CATEGORIES)
from goldenbraid.utils import has_rec_sites
from goldenbraid.tags import TARGET_DICOT, TARGET_MONOCOT


class DomesticationForm(forms.Form):
    choices = [('', '')]
    for category_name in CATEGORIES.keys():
        choices.append((category_name, category_name))
    intron_label = 'The secuence has introns in lowercase'
    with_intron = forms.BooleanField(label=intron_label, required=False)
    category = forms.CharField(max_length=100,
                               label='Choose a category to domesticate to',
                               widget=Select(choices=choices), required=False)
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
            if category in ('13-14-15-16 (CDS)', '13 (SP)', '12 (NT)',
                             '13-14-15 (CDS)'):
                if not _seq_has_codon_start(seq.seq, with_intron):
                    msg = 'The provided seq must start with start codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in ('13-14-15-16 (CDS)', '14-15-16 (CDS)', '16 (CT)'):
                if not _seq_has_codon_end(seq.seq, with_intron):
                    msg = 'The provided seq must end with a end codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in ('13-14-15-16 (CDS)', '13 (SP)', '12 (NT)',
                            '13-14-15 (CDS)', '14-15-16 (CDS)', '16 (CT)'):
                if not _is_seq_3_multiple(seq.seq, with_intron):
                    msg = 'The provided seq must be multiple of three in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in ('12-13 (GOI)'):
                if len(seq) > 500:
                    msg = 'The provided seq must have less than 500 nucleotides in'
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in ('13-14-15 (CDS)', '12 (NT)', '13 (SP)'):
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
        self._multi_field_validation()
        if self._errors:
            del self.cleaned_data

    @staticmethod
    def _data_in(dictionary, key):
        return True if key in dictionary and dictionary[key] else False

    def _multi_field_validation(self):
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


class DomesticationCrisprForm(DomesticationForm):
    choices = [('', '')]
    for category_name in CRYSPER_TARGETS_TO_DOMESTICATE:
        choices.append((category_name, category_name))
    category = forms.CharField(max_length=100,
                               label='Choose a category to domesticate to',
                               widget=Select(choices=choices), required=False)
    seq = forms.CharField(max_length=30, label='Add your sequence',
                          required=True)
    prefix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)
    suffix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)

    def clean_seq(self):
        seq = self.cleaned_data['seq']
        category = self.cleaned_data.get('category', None)
        if len(seq) < 20:
            raise ValidationError('Seq length must be at least 20')

        if len(seq) != 20 and category is not None:
            msg = 'CRISPR target seq length must be 20 nucleotides'
            raise ValidationError(msg)

        if not _seq_is_dna(seq):
            raise ValidationError('Seq must only contain ACTG')
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
