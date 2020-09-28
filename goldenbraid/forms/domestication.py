from tempfile import NamedTemporaryFile
from django import forms
from django.forms.widgets import Select
from django.core.exceptions import ValidationError
from django.forms.utils import ErrorDict

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from goldenbraid.settings import (CATEGORIES, MINIMUN_PCR_LENGTH,
                                  CRYSPER_TARGETS_TO_DOMESTICATE,
                                  CRISPR_CAS9_SINGLE_TO_DOMESTICATE,
                                  CRYSPER_CATEGORIES,
                                  OPTIONAL_DOMEST_ENZYMES,
                                  MANDATORY_DOMEST_ENZYMES,
                                  CRISPR_MULTIPLEXING_TARGET,
                                  FUNGAL_CATEGORIES)
from goldenbraid.utils import has_rec_sites
from goldenbraid.tags import (TARGET_DICOT, TARGET_MONOCOT, CDS, CDS1, NTAG,
                              CDS1_CDS2, CDS2_CTAG, CTAG, TARGET_CAS12A, FUNGAL_CDS,
                              FORWARD_MARKER, REVERSE_MARKER)

def cas12a_seq_is_valid(seq):
        if len(seq) < 20 or len(seq) > 23:
            raise ValidationError('Seq length must be between 20 and 23 nucleotides')

        if not _seq_is_dna(seq):
            raise ValidationError('Seq must only contain ACTG')
        if has_rec_sites(seq):
            msg = 'This sequence can not be domesticated.'
            msg += ' It has internal restriction sites'
            raise ValidationError(msg)
        return True


class DomesticationCas12Mult2XForm(forms.Form):
    target_1 = forms.CharField(max_length=30, label='Add target 1',
                          required=True)
    target_2 = forms.CharField(max_length=30, label='Add target 2',
                          required=True)
    prefix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)
    suffix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)


    def clean_target_1(self):
        target_1 = self.cleaned_data['target_1']
        if cas12a_seq_is_valid(target_1):
            return target_1

    def clean_target_2(self):
        target_2 = self.cleaned_data['target_2']
        if cas12a_seq_is_valid(target_2):
            return target_2


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
        if self._errors:
            del self.cleaned_data

class DomesticationCas12Mult3XForm(DomesticationCas12Mult2XForm):
    target_3 = forms.CharField(max_length=30, label='Add target 3',
                               required=True)
    def clean_target_3(self):
        target_3 = self.cleaned_data['target_3']
        if cas12a_seq_is_valid(target_3):
            return target_3


class DomesticationCas12Mult4XForm(DomesticationCas12Mult3XForm):
    target_4 = forms.CharField(max_length=30, label='Add target 4',
                               required=True)
    def clean_target_4(self):
        target_4 = self.cleaned_data['target_4']
        if cas12a_seq_is_valid(target_4):
            return target_4

class DomesticationCas12Mult5XForm(DomesticationCas12Mult4XForm):
    target_5 = forms.CharField(max_length=30, label='Add target 5',
                               required=True)
    def clean_target_5(self):
        target_5 = self.cleaned_data['target_5']
        if cas12a_seq_is_valid(target_5):
            return target_5

class DomesticationCas12Mult6XForm(DomesticationCas12Mult5XForm):
    target_6 = forms.CharField(max_length=30, label='Add target 6',
                               required=True)
    def clean_target_6(self):
        target_6 = self.cleaned_data['target_6']
        if cas12a_seq_is_valid(target_6):
            return target_6


class DomesticationForm(forms.Form):
    super_choice = [('', '')]
    choices = []
    for category_name in CATEGORIES.keys():
        choices.append((category_name, category_name))
    super_choice.append(('GB parts', choices))
    # Add crispr Proms and TERM
    choices = []
    for category_name in CRYSPER_CATEGORIES.keys():
        if category_name in CRYSPER_TARGETS_TO_DOMESTICATE:
            continue
        choices.append((category_name, category_name))
    super_choice.append(('Crispr parts', choices))
    intron_label = 'The sequence has introns in lowercase'
    with_intron = forms.BooleanField(label=intron_label, required=False)
    category = forms.CharField(max_length=100,
                               label='Choose a category to domesticate to',
                               widget=Select(choices=super_choice),
                               required=False)
    seq = forms.FileField(max_length=100,
                          label='Add a genbank or a fast file', required=False)
    residues = forms.CharField(widget=forms.HiddenInput(), required=False)
    prefix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)
    suffix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)

    enzyme_choices = [(enzy, enzy)for enzy in MANDATORY_DOMEST_ENZYMES]
    enzyme_choices += [(enzy, enzy)for enzy in OPTIONAL_DOMEST_ENZYMES]

    enzymes = forms.MultipleChoiceField(widget=forms.CheckboxSelectMultiple,
                                        choices=enzyme_choices,
                                        label='Enzymes', required=False,
                                        initial=MANDATORY_DOMEST_ENZYMES)

    def clean_category(self):
        category_name = self.cleaned_data['category']
        if not category_name:
            del self.cleaned_data['category']
            return
        categories_to_check = list(CATEGORIES.keys()) + list(CRYSPER_CATEGORIES.keys())
        if category_name not in categories_to_check:
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
            if category in (CDS, CDS1, NTAG, CDS1_CDS2):
                if not _seq_has_codon_start(seq.seq, with_intron):
                    msg = 'The provided seq must start with start codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in (CDS, CTAG):
                if not _seq_has_codon_end(seq.seq, with_intron):
                    msg = 'The provided seq must end with a end codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in (CDS, CDS1, NTAG, CDS1_CDS2, CDS2_CTAG, CTAG):
                if not _is_seq_3_multiple(seq.seq, with_intron):
                    msg = 'The provided seq must be multiple of three in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)

            if category in (CDS1_CDS2, NTAG, CDS1):
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
        content = inmemoryfile.chunks().__next__()
        self.cleaned_data['seq'].seek(0)
        if content.startswith(b'LOCUS') or content.startswith(b'>'):
            pass
        else:
            msg = 'The given file must be a fasta or a genbank file'
            raise ValidationError(msg)
        format_ = "fasta" if content.startswith(b'>') else "gb"
        try:
            temp_file = NamedTemporaryFile()
            temp_file.write(self.cleaned_data['seq'].read())
            temp_file.flush()
            seq = SeqIO.read(temp_file.name, format=format_)

        except ValueError:
            raise ValidationError('File malformed')

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
        print(self._errors)
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

            except ValidationError as e:
                self._errors['seq'] = self.error_class(e.messages)
                for name in ('seq', 'residues'):
                    if name in self.cleaned_data:
                        del self.cleaned_data[name]
                return

            try:
                if (self._data_in(cleaned_data, 'seq') and
                        self._data_in(cleaned_data, 'residues')):
                    raise ValidationError('Form can not accept File and residues')

            except ValidationError as e:
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
        except ValidationError as e:
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
        except ValidationError as e:
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
        except ValidationError as e:
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



class FungalDomesticationForm(forms.Form):
    super_choice = [('', '')]
    choices = []
    for category_name in FUNGAL_CATEGORIES.keys():
        choices.append((category_name, category_name))
    super_choice.append(('GB parts', choices))
    noncoding_label = 'The noncoding sequence is in lowercase'
    with_noncoding_seq = forms.BooleanField(label=noncoding_label, required=False)
    category = forms.CharField(max_length=100,
                               label='Choose a category to domesticate to',
                               widget=Select(choices=super_choice),
                               required=False)
    seq = forms.FileField(max_length=100,
                          label='Add a genbank or a fast file', required=False)
    residues = forms.CharField(widget=forms.HiddenInput(), required=False)
    prefix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)
    suffix = forms.CharField(max_length=4,
                             label='custom prefix', required=False)

    enzyme_choices = [(enzy, enzy)for enzy in MANDATORY_DOMEST_ENZYMES]
    enzyme_choices += [(enzy, enzy)for enzy in OPTIONAL_DOMEST_ENZYMES]

    enzymes = forms.MultipleChoiceField(widget=forms.CheckboxSelectMultiple,
                                        choices=enzyme_choices,
                                        label='Enzymes', required=False,
                                        initial=MANDATORY_DOMEST_ENZYMES)

    def clean_category(self):
        category_name = self.cleaned_data['category']
        if not category_name:
            del self.cleaned_data['category']
            return
        categories_to_check = list(FUNGAL_CATEGORIES.keys())
        if category_name not in categories_to_check:
            raise ValidationError('You must choose a valid category')
        return category_name

    def _seq_validation(self, seq):
        with_noncoding_seq = self.cleaned_data['with_noncoding_seq']
        if not _seq_is_dna(seq.seq):
            msg = 'The given file contains seqs with not allowed nucleotides'
            msg += ' ATGC'
            raise ValidationError(msg)
        if len(seq) < MINIMUN_PCR_LENGTH + 20:
            msg = 'Given seq must be at least 70 base pairs'
            raise ValidationError(msg)
        if self._data_in(self.cleaned_data, 'category'):
            category = self.cleaned_data['category']
            if category == FORWARD_MARKER:
                upper_seq = get_upper_nucls(seq.seq)
                if not upper_seq:
                    msg = 'CDS part of the sequence should be in uppercase '
                    msg += 'in order to use this category'
                    raise ValidationError(msg)
            if category in (FUNGAL_CDS, FORWARD_MARKER):
                if not _seq_has_codon_start(seq.seq,with_noncoding_seq):
                    msg = 'The provided seq must start with start codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in (FUNGAL_CDS, FORWARD_MARKER):
                if not _seq_has_codon_end(seq.seq, with_noncoding_seq):
                    msg = 'The provided seq must end with a end codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category in (FUNGAL_CDS, FORWARD_MARKER):
                if not _is_seq_3_multiple(seq.seq, with_noncoding_seq):
                    msg = 'The provided seq must be multiple of three in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
            if category == REVERSE_MARKER:
                reverse_seq = seq.seq.reverse_complement()
                upper_seq = get_upper_nucls(seq.seq)
                if not upper_seq:
                    msg = 'CDS part of the sequence should be in uppercase '
                    msg += 'in order to use this category'
                    raise ValidationError(msg)
                if not _seq_has_codon_start(reverse_seq, with_noncoding_seq):
                    msg = 'The provided seq must start with start codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
                if not _seq_has_codon_end(reverse_seq, with_noncoding_seq):
                    msg = 'The provided seq must end with a end codon in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)
                if not _is_seq_3_multiple(reverse_seq, with_noncoding_seq):
                    msg = 'The provided seq must be multiple of three in '
                    msg += 'order to use as choosen category'
                    raise ValidationError(msg)

        return seq

    def clean_seq(self):
        inmemoryfile = self.cleaned_data['seq']
        if not inmemoryfile:
            del self.cleaned_data['seq']
            return inmemoryfile
        content = inmemoryfile.chunks().__next__()
        self.cleaned_data['seq'].seek(0)
        if content.startswith(b'LOCUS') or content.startswith(b'>'):
            pass
        else:
            msg = 'The given file must be a fasta or a genbank file'
            raise ValidationError(msg)
        format_ = "fasta" if content.startswith(b'>') else "gb"
        try:
            temp_file = NamedTemporaryFile()
            temp_file.write(self.cleaned_data['seq'].read())
            temp_file.flush()
            seq = SeqIO.read(temp_file.name, format=format_)

        except ValueError:
            raise ValidationError('File malformed')

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
        print(self._errors)
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

            except ValidationError as e:
                self._errors['seq'] = self.error_class(e.messages)
                for name in ('seq', 'residues'):
                    if name in self.cleaned_data:
                        del self.cleaned_data[name]
                return

            try:
                if (self._data_in(cleaned_data, 'seq') and
                        self._data_in(cleaned_data, 'residues')):
                    raise ValidationError('Form can not accept File and residues')

            except ValidationError as e:
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
        except ValidationError as e:
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
        except ValidationError as e:
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
        except ValidationError as e:
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

    def clean_description(self):
        description = self.cleaned_data['description']
        description = textwrap.fill(description, width=70)
        return description


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
        print(self.cleaned_data)
        category = self.cleaned_data.get('category', None)
        print(category)
        if category is None:
            print(len(seq))
            if len(seq) > 23 or len(seq) < 20:
                msg = 'To domesticate with the given target type, the CRISPR target '
                msg += 'size must between 20 and 23 nt long'
                raise ValidationError(msg)

        elif len(seq) != 20:
            msg = 'CRISPR target seq length must be 20 nucleotides'
            raise ValidationError(msg)
        
        if len(seq) < 20:
            raise ValidationError('Seq length must be at least 20')

        if not _seq_is_dna(seq):
            raise ValidationError('Seq must only contain ACTG')
        if has_rec_sites(seq):
            msg = 'This sequence can not be domesticated.'
            msg += ' It has internal restriction sites'
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

    def _multi_field_validation(self):
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


class DomesticationCas9SingleCrisprForm(DomesticationCrisprForm):
    choices = [('', '')]
    for category_name in CRISPR_CAS9_SINGLE_TO_DOMESTICATE:
        choices.append((category_name, category_name))
    category = forms.CharField(max_length=100,
                               label='Choose a category to domesticate to',
                               widget=Select(choices=choices), required=False)

class DomesticationCas9MultiplexCrisprForm(DomesticationCrisprForm):
    choices = [(CRISPR_MULTIPLEXING_TARGET, CRISPR_MULTIPLEXING_TARGET)]
    category = forms.CharField(max_length=100,
                               label='Choose a category to domesticate to',
                               widget=Select(choices=choices), required=False)
    
    def clean_seq(self):
        seq = self.cleaned_data['seq']
        print(self.cleaned_data)
        category = self.cleaned_data.get('category', None)
        print(category)
        if len(seq) != 20:
            msg = 'CRISPR target seq length must be 20 nucleotides'
            raise ValidationError(msg)

        if not _seq_is_dna(seq):
            raise ValidationError('Seq must only contain ACTG')
        if has_rec_sites(seq):
            msg = 'This sequence can not be domesticated.'
            msg += ' It has internal restriction sites'
            raise ValidationError(msg)
        return seq

    class Meta:
        exclude = ['category', ]

class DomesticationCas12SingleCrisprForm(DomesticationCrisprForm):
    choices = [(CRISPR_MULTIPLEXING_TARGET, CRISPR_MULTIPLEXING_TARGET)]
    category = forms.CharField(max_length=100,
                               label='Choose a category to domesticate to',
                               widget=Select(choices=choices), required=False)
    
