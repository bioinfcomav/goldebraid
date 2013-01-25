import re

from django.template.context import RequestContext
from django.core.context_processors import csrf
from django import forms
from django.shortcuts import render_to_response
from django.core.exceptions import ValidationError
from django.db.utils import IntegrityError
from django.http import HttpResponseServerError
from Bio import SeqIO
from Bio.Seq import Seq

from goldenbraid.models import Cvterm, Feature, Db, Dbxref, Featureprop
from goldenbraid.settings import DB, REBASE_FILE
from goldenbraid.tags import (GOLDEN_DB, VECTOR_TYPE_NAME,
                              DESCRIPTION_TYPE_NAME, ENZYME_IN_TYPE_NAME,
                              ENZYME_OUT_TYPE_NAME, RESISTANCE_TYPE_NAME,
                              REFERENCE_TYPE_NAME)


def _parse_rebase_file(fpath):
    'It parses the rebase enzyme file and return a list with all the enzymes'
    enzymes = {}
    enz_name = None
    for line in  open(fpath):
        line = line.strip()
        if not line:
            continue
        if line.startswith('<1>'):
            enz_name = line[3:]
        if line.startswith('<3>'):
            if enz_name is None:
                raise RuntimeError()
            enzymes[enz_name] = line[3:]
            enz_name = None
    return enzymes


class FeatureForm(forms.Form):
    'Form to add features to db'
    name = forms.CharField(max_length=255, required=False)
    description = forms.CharField(max_length=255, required=False)
    reference = forms.CharField(max_length=255, required=False)
    type = forms.CharField()
    vector = forms.CharField(required=False)
    enzyme_in = forms.CharField(required=False)
    # TODO we have to change the widgte to allow multiple values
    # Now we do spliting the text with a comma
    enzyme_out = forms.CharField(required=False)
    resistance = forms.CharField(required=False)

    gbfile_label = 'Select a GenBank-formatted local file on your computer'
    gbfile = forms.FileField(label=gbfile_label, required=True)

#    def clean_uniquename(self):
#        'It checks that the unique name is unique in database'
#        uniquename = self.cleaned_data['uniquename']
#        try:
#            Feature.objects.using(DB).get(uniquename=uniquename)
#        except Feature.DoesNotExist:
#            return uniquename
#        raise ValidationError('There is already a feature with this uniquename')

    def clean_type(self):
        'It validates the type field'
        type_str = self.cleaned_data['type']
        try:
            Cvterm.objects.using(DB).get(name=type_str)
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
                vector_type = Cvterm.objects.using(DB).get(name=VECTOR_TYPE_NAME)
                Feature.objects.using(DB).get(uniquename=vector,
                                              type=vector_type)
            except Feature.DoesNotExist:
                raise ValidationError('The given vector does not exist')

        else:
            if vector:
                raise ValidationError('A vector does not have a vector')

        return vector

    def _validate_enzyme(self, kind):
        '''It validates the vector.

        If the feature is a vector it does not validate anything
        if feature is not a vector it validates that the vector
        is in the database'''
        # TODO. change how we deal with two enzymes for the same field
        enzymes = self.cleaned_data['enzyme_{}'.format(kind)].split(',')
        error_in_type = self.errors.get('type', False)
        error_in_vector = self.errors.get('vector', False)

        if error_in_type or error_in_vector:
            return enzymes

        type_ = self.cleaned_data['type']
        if type_ == VECTOR_TYPE_NAME:
            if enzymes[0] == u'':
                err = 'A vector must have a enzyme {}'.format(kind)
                raise ValidationError(err)
            existing_enzymes = _parse_rebase_file(REBASE_FILE)
            errors = []
            for enzyme in enzymes:
                if enzyme not in existing_enzymes.keys():
                    err = 'This enzyme: {} is not a known enzyme'
                    err = err.format(enzyme)
                    errors.append(err)
            if errors:
                raise ValidationError('\n'.join(errors))

        else:
            if enzymes:
                err = 'Only vectors have enzyme {}'.format(kind)
                raise ValidationError(err)

        return enzymes

    def clean_enzyme_in(self):
        '''It validates the in enzyme'''
        return self._validate_enzyme('in')

    def clean_enzyme_out(self):
        '''It validates the out enzyme'''
        return self._validate_enzyme('out')

    # Validate that a vector must have a resistance,
    # only vectors have resistance and
    # if the resistance exists or not

    def clean_resistance(self):
        '''It validates the in resistance'''
        resistance = self.cleaned_data['resistance']
        error_in_type = self.errors.get('type', False)
        error_in_vector = self.errors.get('vector', False)

        if error_in_type or error_in_vector:
            return resistance

        type_ = self.cleaned_data['type']
        if type_ == VECTOR_TYPE_NAME:
            if not resistance:
                raise ValidationError('A vector must have a resistance')

        else:
            if resistance:
                raise ValidationError('Only vectors have resistance')
        return resistance


def _search_rec_sites(seq, rec_site):
    """It looks for the rec sites in the string"""
    # look for the rec_site in the seq
    elong_size = 10
    seq_elonged = seq + seq[:elong_size]
    residues = str(seq_elonged)

    finded_site_indexes = [m.start() for m in re.finditer(rec_site.upper(),
                                                          residues.upper())]
    corrected_site_indexes = set()
    for site in finded_site_indexes:
        if site > len(seq):
            site -= len(seq)
        corrected_site_indexes.add(site)

    if len(corrected_site_indexes) > 2:
        raise RuntimeError('rec site found more than twice')
    corrected_site_indexes = list(corrected_site_indexes)
    corrected_site_indexes.sort()
    return corrected_site_indexes


def _choose_rec_sites(forward_sites, rev_sites):
    'It chooses the forward and reverse site'
    len_for = len(forward_sites)
    len_rev = len(rev_sites)
    if len_for == len_rev and len_for == 1:
        return forward_sites[0], rev_sites[0]
    elif len_for == len_rev and len_for > 2:
        msg = "We can't have this number of sites: {0}"
        msg = msg.format(len_for + len_rev)
        raise RuntimeError(msg)
    elif len_for < len_rev:
        forw_site = forward_sites[0]
        rev_site = None
    else:
        rev_site = rev_sites[0]
        forw_site = None

    all_sites = rev_sites + forward_sites
    all_sites.sort()
    if rev_site is None:
        index_in_all = all_sites.index(forw_site)
        try:
            rev_site = all_sites[index_in_all + 1]
        except IndexError:
            rev_site = all_sites[0]

    elif forw_site is None:
        index_in_all = all_sites.index(rev_site)
        try:
            forw_site = all_sites[index_in_all - 1]
        except IndexError:
            forw_site = all_sites[len(all_sites) - 1]

    return forw_site, rev_site


def get_prefix_and_suffix_index(seq, enzyme):
    'it gets the prefix and the suffix indexes of the feature seq'
    restriction_site = _parse_rebase_file(REBASE_FILE)[enzyme]
    if '^' in restriction_site:
        raise NotImplementedError
    rec_site, cut_site = restriction_site.split('(')
    forw_cut_delta, rev_cut_delta = cut_site.rstrip(')').split('/')
    forw_cut_delta, rev_cut_delta = int(forw_cut_delta), int(rev_cut_delta)

    forw_sites = _search_rec_sites(seq, rec_site)
    rec_seq = Seq(rec_site)
    rec_seq.reverse_complement()
    rev_sites = _search_rec_sites(seq, str(rec_seq.reverse_complement()))
    forw_site, rev_site = _choose_rec_sites(forw_sites, rev_sites)
    prefix_index, suffix_index = _pref_suf_index_from_rec_sites(seq,
                                                                forw_site,
                                                                rev_site,
                                                                rec_site,
                                                                forw_cut_delta,
                                                                rev_cut_delta)

    return prefix_index, suffix_index, rev_cut_delta - forw_cut_delta


def get_prefix_and_suffix(seq, enzyme):
    'it gets the prefix and the suffix of the feature seq'
    prefix_index, suffix_index, prefix_size = get_prefix_and_suffix_index(seq,
                                                                        enzyme)
    return _get_pref_suff_from_index(seq, prefix_index, suffix_index,
                                     prefix_size)


def _get_pref_suff_from_index(seq, prefix_index, suffix_index, prefix_size):

    prefix = seq[prefix_index:prefix_index + prefix_size]
    suffix = seq[suffix_index:suffix_index + prefix_size]

    if len(suffix) < prefix_size:
        remaining = prefix_size - len(suffix)
        suffix += seq[0:remaining]

    return str(prefix), str(suffix)


def _pref_suf_index_from_rec_sites(seq, forw_site, rev_site, rec_site,
                                   forw_cut_delta, rev_cut_delta):

    prefix_index = forw_site + len(rec_site) + forw_cut_delta
    if prefix_index >= len(seq):
        prefix_index = prefix_index - len(seq)

    suffix_index = rev_site - rev_cut_delta
    if suffix_index < 0:
        suffix_index = len(seq) - abs(suffix_index)
    return prefix_index, suffix_index


def add_feature(database, name, type_name, vector, genbank, props):
    'it adds a feature to the database'
    seq = SeqIO.read(genbank, 'gb')
    residues = str(seq.seq)
    name = name
    uniquename = seq.id
    type_ = Cvterm.objects.using(database).get(name=type_name)
    db = Db.objects.using(database).get(name=GOLDEN_DB)
    try:
        dbxref = Dbxref.objects.using(database).create(db=db,
                                                       accession=uniquename)
    except IntegrityError as error:
        raise IntegrityError('feature already in db' + str(error))

    vector_type = Cvterm.objects.using(database).get(name=VECTOR_TYPE_NAME)
    if vector and type_ == vector_type:
        # already checked in form validation
        raise RuntimeError("a vector feature can't  have a vector")
    if vector:
        vector = Feature.objects.using(database).get(uniquename=vector,
                                               type=vector_type)
    else:
        vector = None
    try:
        feature = Feature.objects.using(database).create(uniquename=uniquename,
                                                   name=name, type=type_,
                                                   residues=residues,
                                                  dbxref=dbxref, vector=vector)
    except IntegrityError as error:
        raise IntegrityError('feature already in db' + str(error))
    for type_name, values in props.items():
        try:
            prop_type = Cvterm.objects.using(DB).get(name=type_name)
        except Cvterm.DoesNotExist:
            msg = 'Trying to add a property which cvterm does not exist: {0}'
            msg = msg.format(type_name)
            raise RuntimeError(msg)
        rank = 0
        for value in values:
            Featureprop.objects.using(DB).create(feature=feature,
                                                 type=prop_type,
                                                 value=value, rank=rank)
            rank += 1
    if type_ == vector_type:
        enzyme = props[ENZYME_IN_TYPE_NAME][0]
    else:
        enzyme = vector.enzyme_out[0]
    prefix, suffix = get_prefix_and_suffix(residues, enzyme)
    feature.prefix = prefix
    feature.suffix = suffix
    feature.save(using=database)

    return feature


def add_feature_from_form(form_data):
    'With this function we add a feature to the database'
    props = {}
    vector_type = Cvterm.objects.using(DB).get(name=VECTOR_TYPE_NAME)
    feature_type_name = form_data['type']
    if form_data['description']:
        props[DESCRIPTION_TYPE_NAME] = [form_data['description']]
    if form_data['reference']:
        props[REFERENCE_TYPE_NAME] = [form_data['reference']]
    if feature_type_name == vector_type.name:
        props[ENZYME_IN_TYPE_NAME] = form_data['enzyme_in']
        props[ENZYME_OUT_TYPE_NAME] = form_data['enzyme_out']
        props[RESISTANCE_TYPE_NAME] = [form_data['resistance']]

    feature = add_feature(database=DB,
                          name=form_data['name'], type_name=feature_type_name,
                          vector=form_data['vector'],
                          genbank=form_data['gbfile'],
                          props=props)

    return feature


def add_feature_view(request):
    'The add feature view'
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    else:
        request_data = None

    if request_data:
        form = FeatureForm(request_data, request.FILES)
        if form.is_valid():
            feat_form_data = form.cleaned_data
            try:
                feature = add_feature_from_form(feat_form_data)
            except IntegrityError as error:
                if 'feature already in db' in error:
                    # TODO choose a template
                    return render_to_response('feature_template.html',
                                              {},
                                    context_instance=RequestContext(request))
                else:
                    return HttpResponseServerError()
            except Exception as error:
                return HttpResponseServerError()
            # if everithing os fine we show the just added feature
            return render_to_response('feature_template.html',
                                          {'feature': feature},
                                          context_instance=RequestContext(request))

    else:
        form = FeatureForm()
    context['form'] = form
    template = 'feature_add_template.html'
    return render_to_response(template, context)


def feature_view(request, uniquename):
    'The feature view'
    try:
        feature = Feature.objects.using(DB).get(uniquename=uniquename)
    except Feature.DoesNotExist:
        feature = None
    return render_to_response('feature_template.html', {'feature': feature},
                              context_instance=RequestContext(request))

