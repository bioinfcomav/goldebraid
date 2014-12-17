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

import re
import os

from django.template.context import RequestContext
from django.core.context_processors import csrf
from django.shortcuts import render_to_response, redirect
from django.db.utils import IntegrityError
from django.http import HttpResponseServerError
from django.core.files import File
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required
from django.db import transaction
from django.http.response import HttpResponseBadRequest
from django.contrib.admin.views.decorators import staff_member_required
from django.core.exceptions import MultipleObjectsReturned

from Bio import SeqIO
from Bio.Seq import Seq

from goldenbraid.models import (Cvterm, Feature, Db, Dbxref, Featureprop,
                                FeaturePerm, FeatureRelationship)
from goldenbraid.settings import REBASE_FILE
from goldenbraid.tags import (GOLDEN_DB, VECTOR_TYPE_NAME,
                              DESCRIPTION_TYPE_NAME, ENZYME_IN_TYPE_NAME,
                              REFERENCE_TYPE_NAME, ENZYME_OUT_TYPE_NAME,
                              RESISTANCE_TYPE_NAME, DERIVES_FROM)
from goldenbraid.forms import (FeatureForm, FeatureManagementForm,
                               get_all_vectors_as_choices, VectorForm)


def parse_rebase_file(fpath):
    'It parses the rebase enzyme file and return a list with all the enzymes'
    enzymes = {}
    enz_name = None
    for line in open(fpath):
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
    if not corrected_site_indexes:
        raise RuntimeError("No rec_site")
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
    restriction_site = parse_rebase_file(REBASE_FILE)[enzyme]
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
    try:
        prefix_index, suffix_index, prefix_size = get_prefix_and_suffix_index(seq, enzyme)
    except RuntimeError as error:
        if str(error) == "No rec_site":
            return None, None
        else:
            raise
    return _get_pref_suff_from_index(seq, prefix_index, suffix_index,
                                     prefix_size)


def _get_pref_suff_from_index(seq, prefix_index, suffix_index, prefix_size):

    prefix = seq[prefix_index:prefix_index + prefix_size]
    suffix = seq[suffix_index:suffix_index + prefix_size]

    if len(suffix) < prefix_size:
        remaining = prefix_size - len(suffix)
        suffix += seq[0:remaining]
    if len(prefix) < prefix_size:
        remaining = prefix_size - len(prefix)
        prefix += seq[0:remaining]

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


def get_or_create_feature_relationship(object_, subject):
    derives_from = Cvterm.objects.get(name=DERIVES_FROM)
    try:
        FeatureRelationship.objects.get(subject=subject,
                                        type=derives_from, object=object_)
    except FeatureRelationship.DoesNotExist:
        FeatureRelationship.objects.create(subject=subject,
                                           type=derives_from, object=object_)


def add_relations(feature, seq):
    children = _parse_children_relations_from_gb(seq)
    if not children:
        return

    for child in children:
        child = Feature.objects.get(uniquename=child)
        get_or_create_feature_relationship(object_=feature, subject=child)


def _parse_children_relations_from_gb(seq):
    definition = seq.description
    if '(' in definition and ')' in definition:
        match = re.match('\((.+)\)', definition)
        return match.group(1).split(',')
    else:
        return None


def add_feature(name, type_name, vector, genbank, props, owner,
                is_public=False):
    'it adds a feature to the database'
    feature = None
    #transaction.set_autocommit(False)
    try:
        with transaction.atomic():
            seq = SeqIO.read(genbank, 'gb')
            residues = str(seq.seq)
            name = name
            uniquename = seq.id
            type_ = Cvterm.objects.get(name=type_name)
            db = Db.objects.get(name=GOLDEN_DB)
            genbank_file = File(genbank)
            try:
                dbxref = Dbxref.objects.create(db=db, accession=uniquename)
            except IntegrityError as error:
                raise IntegrityError('feature already in db' + str(error))
            vector_type = Cvterm.objects.get(name=VECTOR_TYPE_NAME)
            if vector and type_ == vector_type:
                # already checked in form validation
                raise RuntimeError("a vector feature can't  have a vector")
            if vector:
                vector = Feature.objects.get(uniquename=vector,
                                             type=vector_type)
            else:
                vector = None
            try:
                user = User.objects.get(username=owner)
            except User.DoesNotExist:
                raise RuntimeError('the given user does not exist')
            try:
                feature = Feature.objects.create(uniquename=uniquename,
                                                 name=name, type=type_,
                                                 residues=residues,
                                                 dbxref=dbxref, vector=vector,
                                                 genbank_file=genbank_file)

            except IntegrityError as error:
                raise IntegrityError('feature already in db' + str(error))

            FeaturePerm.objects.create(feature=feature, owner=user,
                                       is_public=is_public)

            feature.add_relations(seq)

            for type_name, values in props.items():
                try:
                    prop_type = Cvterm.objects.get(name=type_name)
                except Cvterm.DoesNotExist:
                    msg = 'Trying to add a property which cvterm does not exist: {0}'
                    msg = msg.format(type_name)
                    raise RuntimeError(msg)
                except MultipleObjectsReturned:
                    for p in Cvterm.objects.filter(name=type_name):
                        print p.name
                        print p.cvterm_id
                        print p.definition
                    print "type_name", type_name
                    print "feature", feature.uniquename
                    raise
                rank = 0
                for value in values:
                    Featureprop.objects.create(feature=feature, type=prop_type,
                                               value=value, rank=rank)
                    rank += 1
            if type_ == vector_type:
                enzyme = props[ENZYME_IN_TYPE_NAME][0]
            else:
                enzyme = vector.enzyme_out[0]
            prefix, suffix = get_prefix_and_suffix(residues, enzyme)
            if prefix is None or suffix is None:
                raise RuntimeError('The given vector is not compatible with this part')

            feature.prefix = prefix
            feature.suffix = suffix
            feature.save()

    except (IntegrityError, RuntimeError):
        if feature:
            os.remove(feature.genbank_file.path)
        transaction.rollback()
        raise
    return feature


def add_vector_from_form(form_data, user):
    'With this function we add a feature to the database'
    props = {}
    props[DESCRIPTION_TYPE_NAME] = [form_data['description']]
    props[ENZYME_IN_TYPE_NAME] = [form_data['enzyme_in']]
    props[ENZYME_OUT_TYPE_NAME] = [form_data['enzyme_out']]
    props[RESISTANCE_TYPE_NAME] = [form_data['resistance']]

    if form_data['reference']:
        props[REFERENCE_TYPE_NAME] = [form_data['reference']]

    vector = add_feature(name=form_data['name'], type_name=VECTOR_TYPE_NAME,
                         vector=None, genbank=form_data['gbfile'],
                         props=props, owner=user)

    return vector


def add_feature_from_form(form_data, user):
    'With this function we add a feature to the database'
    props = {}
    feature_type_name = form_data['type']
    if form_data['description']:
        props[DESCRIPTION_TYPE_NAME] = [form_data['description']]
    if form_data['reference']:
        props[REFERENCE_TYPE_NAME] = [form_data['reference']]

    feature = add_feature(name=form_data['name'], type_name=feature_type_name,
                          vector=form_data['vector'],
                          genbank=form_data['gbfile'],
                          props=props, owner=user)

    return feature


@staff_member_required
def add_vector_view(request):
    'The add feature view'
    context = RequestContext(request)
    context.update(csrf(request))
    request_data_post = request.POST if request.method == 'POST' else None
#    request_data_get = request.GET if request.method == 'GET' else None

    if request_data_post:
        form = VectorForm(request_data_post, request.FILES)
        if form.is_valid():
            vector_form_data = form.cleaned_data
            try:
                vector = add_vector_from_form(vector_form_data, request.user)
            except IntegrityError as error:
                if 'feature already in db' in str(error):
                    # TODO choose a template
                    return render_to_response('feature_exists.html', {},
                                    context_instance=RequestContext(request))
                else:
                    return HttpResponseServerError(str(error))

            except Exception as error:
                return render_to_response('goldenbraid_info.html',
                                          {'title': 'Error',
                                           'info': str(error)},
                                      context_instance=RequestContext(request))
            # if everithing os fine we show the just added feature
            return redirect(vector.url)

    else:
        form = VectorForm()

    context['form'] = form
    template = 'vector_add_template.html'
    return render_to_response(template, context)


@login_required
def add_feature_view(request):
    'The add feature view'
    context = RequestContext(request)
    context.update(csrf(request))
    request_data_post = request.POST if request.method == 'POST' else None
#    request_data_get = request.GET if request.method == 'GET' else None

    if request_data_post:
        form = FeatureForm(request_data_post, request.FILES)
        if form.is_valid():
            feat_form_data = form.cleaned_data
            try:
                feature = add_feature_from_form(feat_form_data, request.user)
            except IntegrityError as error:
                if 'feature already in db' in str(error):
                    # TODO choose a template
                    return render_to_response('feature_exists.html', {},
                                    context_instance=RequestContext(request))
                else:
                    return HttpResponseServerError(str(error))

            except Exception as error:
                return render_to_response('goldenbraid_info.html',
                                          {'title': 'Error',
                                           'info': str(error)},
                                      context_instance=RequestContext(request))
            # if everithing os fine we show the just added feature
            return redirect(feature.url)

    else:
        form = FeatureForm()
        form.fields['vector'].widget.choices = get_all_vectors_as_choices(request.user)

    context['form'] = form
    template = 'feature_add_template.html'
    return render_to_response(template, context)


def feature_view(request, uniquename):
    'The feature view'
    try:
        feature = Feature.objects.get(uniquename=uniquename)
    except Feature.DoesNotExist:
        feature = None

    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if feature is None:
        return render_to_response('goldenbraid_info.html',
                                  {'title': 'Feature not exist',
                                   'info': 'This feature ({0}) does not exist in the database'.format(uniquename)},
                                  context_instance=RequestContext(request))

    if not request_data:
        if (feature.is_public or
           (request.user.is_staff or request.user == feature.owner)):
            return render_to_response('feature_template.html',
                                      {'feature': feature},
                                      context_instance=RequestContext(request))
        else:
            return render_to_response('goldenbraid_info.html',
                                      {'title': 'Not Allowed',
                                       'info': 'You are not allowed to view this feature'},
                                      context_instance=RequestContext(request))

    else:
        form = FeatureManagementForm(request_data)
        if form.is_valid():
            form_data = form.cleaned_data
            action = form_data['action']
            feature = Feature.objects.get(uniquename=form_data['feature'])
            if action == 'delete':
                if request.user.is_staff or request.user == feature.owner:
                    file_path = feature.genbank_file.path
                    feature.delete()
                    os.remove(file_path)
                    return render_to_response('goldenbraid_info.html',
                                              {'title': 'Feature deleted',
                                               'info': 'Feature Deleted'},
                                    context_instance=RequestContext(request))
                else:
                    return render_to_response('goldenbraid_info.html',
                                              {'title': 'Not Allowed',
                        'info': 'You are not allowed to delete this feature'},
                                    context_instance=RequestContext(request))

            elif action in 'make_public' or 'make_private':
                if request.user.is_staff:
                    if action == 'make_public' and not feature.is_public:
                        featperm = FeaturePerm.objects.get(feature=feature)
                        featperm.is_public = True
                        featperm.save()
                    elif action == 'make_private' and feature.is_public:
                        featperm = FeaturePerm.objects.get(feature=feature)
                        featperm.is_public = False
                        featperm.save()
                    else:
                        raise RuntimeError('bad conbinations of input request')
                    return render_to_response('feature_template.html',
                                              {'feature': feature,
                                               'info': 'Feature modified'},
                              context_instance=RequestContext(request))
                else:
                    return render_to_response('Goldenbraid_info.html',
                                              {'title': 'Not Allowed',
                        'info': 'You are not allowed to modify this feature'},
                                    context_instance=RequestContext(request))

            else:
                return HttpResponseBadRequest()

