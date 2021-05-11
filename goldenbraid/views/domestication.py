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

import textwrap
from io import StringIO
from tempfile import NamedTemporaryFile
from textwrap import fill, wrap
from os.path import join
from collections import OrderedDict
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

from django.template.context_processors import csrf
from django.template.context import RequestContext
from django.shortcuts import render, redirect
from django.http import HttpResponse, HttpResponseBadRequest, Http404
from django.contrib.auth.decorators import login_required
from django.http.response import HttpResponseServerError
from django.db.utils import IntegrityError

from goldenbraid.domestication import (domesticate, domesticate_for_synthesis,
                                       domestication_crispr, annotate_domesticated_crispr)
from goldenbraid.views.multipartite import assemble_cas12_unit_zero
from goldenbraid.settings import (CATEGORIES, CRYSPER_CATEGORIES,
                                  DOMESTICATED_VECTOR, CRISPR_MULTIPLEXING_TARGET,
                                  TARGET_CAS12A, CAS12_LEVEL_MINUS_ONE_2X,
                                  CAS12_LEVEL_MINUS_ONE_3X,
                                  CAS12_LEVEL_MINUS_ONE_4X,
                                  CAS12_LEVEL_MINUS_ONE_5X,
                                  CAS12_LEVEL_MINUS_ONE_6X,
                                  CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE,
                                  FUNGAL_CATEGORIES)
from goldenbraid.tags import CRISPR_MULTIPLEXING_TARGET, FUNGAL_CDS, REVERSE_MARKER


from goldenbraid.forms.domestication import (DomesticationForm,
                                             DomesticationCas9SingleCrisprForm,
                                             DomesticationCas9MultiplexCrisprForm,
                                             DomesticationCas12SingleCrisprForm,
                                             DomesticationCas12Mult2XForm,
                                             DomesticationCas12Mult3XForm,
                                             DomesticationCas12Mult4XForm,
                                             DomesticationCas12Mult5XForm,
                                             DomesticationCas12Mult6XForm,
                                             FungalDomesticationForm)
from goldenbraid.views.feature import add_feature
from goldenbraid.sbol import convert_to_sbol
from goldenbraid.models import Feature

TMP_DIR = join("/home/golden/devel/gbdb", "files")


def synthesis_view(request):
    return _domestication_view(request, kind='synthesis')


def domestication_view(request):
    return _domestication_view(request, kind='domestication')

def fungal_domestication_view(request):
    return _fungal_domestication_view(request, kind='domestication')


def crispr_view_cas9_single(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        form = DomesticationCas9SingleCrisprForm(request_data, request.FILES)
        if form.is_valid():
            seq = form.cleaned_data['seq']
            category = form.cleaned_data.get('category', None)
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = CRYSPER_CATEGORIES[category][1]
                suffix = CRYSPER_CATEGORIES[category][2]
            new_seq = domestication_crispr(seq, category, prefix, suffix)
            forw_dim = new_seq[:-4]
            rev_dim = new_seq[4:].reverse_complement()
            context = {'category': category, 'prefix': prefix, 'suffix': suffix,
                       'seq': str(new_seq.seq), 'seq_name': new_seq.name,
                       'forw_dim': str(forw_dim.seq), 'rev_dim': str(rev_dim.seq)}
            return render(request, 'crispr_result.html', context=context)
    else:
        form = DomesticationCas9SingleCrisprForm()

    context['form'] = form

    template = 'crispr_cas9_single_template.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def crispr_view_cas9_multiplexing(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        form = DomesticationCas9MultiplexCrisprForm(request_data, request.FILES)
        if form.is_valid():
            seq = form.cleaned_data['seq']
            category = CRISPR_MULTIPLEXING_TARGET
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = CRYSPER_CATEGORIES[category][1]
                suffix = CRYSPER_CATEGORIES[category][2]
            new_seq = domestication_crispr(seq, category, prefix, suffix)
            forw_dim = new_seq[:-4]
            rev_dim = new_seq[4:].reverse_complement()
            context = {'category': category, 'prefix': prefix, 'suffix': suffix,
                       'seq': str(new_seq.seq), 'seq_name': new_seq.name,
                       'forw_dim': str(forw_dim.seq), 'rev_dim': str(rev_dim.seq)}
            return render(request, 'crispr_result.html', context=context)
    else:
        form = DomesticationCas9MultiplexCrisprForm()

    context['form'] = form

    template = 'crispr_cas9_multiplexing_template.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def crispr_view_cas12_multiplexing_2X(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        form = DomesticationCas12Mult2XForm(request_data, request.FILES)
        if form.is_valid():
            targets = [form.cleaned_data['target_1'], form.cleaned_data['target_2']]
            category = CAS12_LEVEL_MINUS_ONE_2X
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[category][1]
                suffix = CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[category][2]
            new_seq, assembled_seq = assemble_cas12_unit_zero(targets.copy(), category)

            context = {'category': category, 'prefix': prefix, 'suffix': suffix,
                       'seq': str(new_seq.seq), 'seq_name': new_seq.name,
                       'reverse_seq': str(new_seq.reverse_complement().seq),
                       'targets': targets}
            return render(request, 'cas12_multiplex_result.html', context=context)
    else:
        form = DomesticationCas12Mult2XForm()

    context['form'] = form

    template = 'cas12_multiplex_2x.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def crispr_view_cas12_multiplexing_3X(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        form = DomesticationCas12Mult3XForm(request_data, request.FILES)
        if form.is_valid():
            targets = [form.cleaned_data['target_1'], form.cleaned_data['target_2'],
                       form.cleaned_data['target_3']]
            category = CAS12_LEVEL_MINUS_ONE_3X
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[category][1]
                suffix = CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[category][2]
            new_seq, assembled_seq = assemble_cas12_unit_zero(targets.copy(), category)

            context = {'category': category, 'prefix': prefix, 'suffix': suffix,
                       'seq': str(new_seq.seq), 'seq_name': new_seq.name,
                       'reverse_seq': str(new_seq.reverse_complement().seq),
                       'targets': targets}
            return render(request, 'cas12_multiplex_result.html', context=context)
    else:
        form = DomesticationCas12Mult3XForm()

    context['form'] = form

    template = 'cas12_multiplex_3x.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def crispr_view_cas12_multiplexing_4X(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        form = DomesticationCas12Mult4XForm(request_data, request.FILES)
        if form.is_valid():
            targets = [form.cleaned_data['target_1'], form.cleaned_data['target_2'],
                       form.cleaned_data['target_3'], form.cleaned_data['target_4']]
            category = CAS12_LEVEL_MINUS_ONE_4X
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[category][1]
                suffix = CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[category][2]
            new_seq, assembled_seq = assemble_cas12_unit_zero(targets.copy(), category)

            context = {'category': category, 'prefix': prefix, 'suffix': suffix,
                       'seq': str(new_seq.seq), 'seq_name': new_seq.name,
                       'reverse_seq': str(new_seq.reverse_complement().seq),
                       'targets': targets}
            return render(request, 'cas12_multiplex_result.html', context=context)
    else:
        form = DomesticationCas12Mult4XForm()

    context['form'] = form

    template = 'cas12_multiplex_4x.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def crispr_view_cas12_multiplexing_5X(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        form = DomesticationCas12Mult5XForm(request_data, request.FILES)
        if form.is_valid():
            targets = [form.cleaned_data['target_1'], form.cleaned_data['target_2'],
                       form.cleaned_data['target_3'], form.cleaned_data['target_4'],
                       form.cleaned_data['target_5']]
            category = CAS12_LEVEL_MINUS_ONE_5X
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[category][1]
                suffix = CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[category][2]
            new_seq, assembled_seq = assemble_cas12_unit_zero(targets.copy(), category)

            context = {'category': category, 'prefix': prefix, 'suffix': suffix,
                       'seq': str(new_seq.seq), 'seq_name': new_seq.name,
                       'reverse_seq': str(new_seq.reverse_complement().seq),
                       'targets': targets}
            return render(request, 'cas12_multiplex_result.html', context=context)
    else:
        form = DomesticationCas12Mult5XForm()

    context['form'] = form

    template = 'cas12_multiplex_5x.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def crispr_view_cas12_multiplexing_6X(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None

    if request_data:
        form = DomesticationCas12Mult6XForm(request_data, request.FILES)
        if form.is_valid():
            targets = [form.cleaned_data['target_1'], form.cleaned_data['target_2'],
                       form.cleaned_data['target_3'], form.cleaned_data['target_4'],
                       form.cleaned_data['target_5'], form.cleaned_data['target_6']]
            category = CAS12_LEVEL_MINUS_ONE_6X
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[category][1]
                suffix = CRISPR_CAS12A_MULTIPLEX_CATEGORIES_LEVEL_MINUS_ONE[category][2]
            new_seq, assembled_seq = assemble_cas12_unit_zero(targets.copy(), category)

            context = {'category': category, 'prefix': prefix, 'suffix': suffix,
                       'seq': str(new_seq.seq), 'seq_name': new_seq.name,
                       'reverse_seq': str(new_seq.reverse_complement().seq),
                       'targets': targets}
            return render(request, 'cas12_multiplex_result.html', context=context)
    else:
        form = DomesticationCas12Mult6XForm()

    context['form'] = form

    template = 'cas12_multiplex_6x.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)

def crispr_view_cas12_single(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        form = DomesticationCas12SingleCrisprForm(request_data, request.FILES)
        if form.is_valid():
            seq = form.cleaned_data['seq']
            category = TARGET_CAS12A
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = CRYSPER_CATEGORIES[category][1]
                suffix = CRYSPER_CATEGORIES[category][2]
            new_seq = domestication_crispr(seq, category, prefix, suffix)
            forw_dim = new_seq[:-4]
            rev_dim = new_seq[4:].reverse_complement()
            context = {'category': category, 'prefix': prefix, 'suffix': suffix,
                       'seq': str(new_seq.seq), 'seq_name': new_seq.name,
                       'forw_dim': str(forw_dim.seq), 'rev_dim': str(rev_dim.seq)}
            return render(request, 'crispr_result.html', context=context)
    else:
        form = DomesticationCas12SingleCrisprForm()

    context['form'] = form

    template = 'crispr_cas12a_single_template.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def crispr_view(request):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        form = DomesticationCrisprForm(request_data, request.FILES)
        if form.is_valid():
            seq = form.cleaned_data['seq']
            category = form.cleaned_data.get('category', None)
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = CRYSPER_CATEGORIES[category][1]
                suffix = CRYSPER_CATEGORIES[category][2]
            new_seq = domestication_crispr(seq, category, prefix, suffix)
            forw_dim = new_seq[:-4]
            rev_dim = new_seq[4:].reverse_complement()
            context = {'category': category, 'prefix': prefix, 'suffix': suffix,
                       'seq': str(new_seq.seq), 'seq_name': new_seq.name,
                       'forw_dim': str(forw_dim.seq), 'rev_dim': str(rev_dim.seq)}
            return render(request, 'crispr_result.html', context=context)
    else:
        form = DomesticationCrisprForm()

    context['form'] = form

    template = 'crispr_template.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def _parse_feature(feature):
    print(feature)
    feature_start = feature.SeqFeature.location.start
    print(feature_start)

def _domestication_view(request, kind):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        form = DomesticationForm(request_data, request.FILES, kind=kind)
        # do domestication
        if form.is_valid():
            seq = form.cleaned_data['seq']
            if seq is None:
                seq = form.cleaned_data['residues']
            category = form.cleaned_data.get('category', None)
            enzymes = form.cleaned_data.get('enzymes', None)
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                try:
                    prefix = CATEGORIES[category][1]
                    suffix = CATEGORIES[category][2]
                except KeyError:
                    prefix = CRYSPER_CATEGORIES[category][1]
                    suffix = CRYSPER_CATEGORIES[category][2]

            with_intron = form.cleaned_data['with_intron']
            with_intron_str = '1' if with_intron else '0'
            if kind == 'domestication':
                try:
                    pcr, new_seq_record = domesticate(seq, category, prefix, suffix, enzymes,
                                                      with_intron=with_intron)
                except RuntimeError as error:
                    context = {'title': 'Can not domesticate sequence', 'info': error}
                    return render(request, 'goldenbraid_info.html',
                                  context=context)
                except ValueError as error:
                    context = {'title': 'Can not domesticate sequence', 
                               'info': error}
                    return render(request, 'goldenbraid_info.html',
                                  context=context)
                new_seq_record_handle = StringIO()
                SeqIO.write(new_seq_record, new_seq_record_handle, "genbank")
                context = {'category': category,
                           'prefix': prefix,
                           'suffix': suffix,
                           'pcrs': pcr,
                           'seq': str(seq.seq),
                           'record': new_seq_record_handle.getvalue(),
                           'seq_name': seq.name,
                           'enzymes': enzymes,
                           'with_intron': with_intron_str}
                return render(request, 'domestication_regular_result.html', context=context)
            elif kind == 'synthesis':
                try:
                    seq_for_syn, prepared_seq = domesticate_for_synthesis(seq,
                                                                          category,
                                                                          prefix,
                                                                          suffix,
                                                                          enzymes,
                                                                          with_intron)
                except RuntimeError as error:
                    context = {'title': 'Can not domesticate sequence',
                               'info': error}
                    return render(request, 'goldenbraid_info.html',
                                  context=context)
                features_qualifiers = []
                for feature in seq.features:
                    features_qualifiers.append(feature.qualifiers)
                prepared_seq_handle = StringIO()
                SeqIO.write(prepared_seq, prepared_seq_handle, "genbank")
                print(prepared_seq_handle)
                context = {'category': category,
                           'prefix': prefix,
                           'suffix': suffix,
                           'seq_syn': seq_for_syn,
                           'seq': str(seq.seq),
                           'record': seq,
                           'record': prepared_seq_handle.getvalue(),
                           'seq_name': prepared_seq.name,
                           'enzymes': enzymes,
                           'with_intron': with_intron_str,
                           'seq_features': seq.features, 
                           'features_qualifiers': features_qualifiers}
                return render(request, 'synthesis_result.html', context=context)
    else:
        form = DomesticationForm()
    context['form'] = form
    context['kind'] = kind

    template = 'domestication_template.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def _fungal_domestication_view(request, kind):
    context = {}
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        form = FungalDomesticationForm(request_data, request.FILES)
        # do domestication
        if form.is_valid():
            seq = form.cleaned_data['seq']
            if seq is None:
                seq = form.cleaned_data['residues']
            category = form.cleaned_data.get('category', None)
            enzymes = form.cleaned_data.get('enzymes', None)
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = FUNGAL_CATEGORIES[category][1]
                suffix = FUNGAL_CATEGORIES[category][2]
            
            if category ==  REVERSE_MARKER:
                reverse_orientation = True
            else:
                reverse_orientation =  False

            with_noncoding_seq = form.cleaned_data['with_noncoding_seq']
            with_noncoding_seq_str = '1' if with_noncoding_seq else '0'
            try:
                pcr = domesticate(seq, category, prefix, suffix, enzymes,
                                  with_intron=with_noncoding_seq,
                                  reverse_orientation=reverse_orientation)[0]
            except RuntimeError as error:
                context = {'title': 'Can not domesticate sequence', 'info': error}
                return render(request, 'goldenbraid_info.html', context=context)
            context = {'category': category,
                       'prefix': prefix,
                       'suffix': suffix,
                       'pcrs': pcr,
                       'seq': str(seq.seq),
                       'seq_name': seq.name,
                       'enzymes': enzymes,
                       'with_noncoding_seq': with_noncoding_seq_str}
            return render(request, 'fungal_domestication_result.html', context=context)
    else:
        form = FungalDomesticationForm()
    context['form'] = form
    context['kind'] = kind

    template = 'fungal_domestication_template.html'
    content_type = None
    return render(request, template, context=context, content_type=content_type)


def fungal_domestication_view_genbank(request):
    def function(seq, category, prefix, suffix, enzymes, with_noncoding_seq, reverse_orientation):
        seq = domesticate(seq, category, prefix, suffix, enzymes,
                          with_intron=with_noncoding_seq, reverse_orientation=reverse_orientation)[1]
        response = HttpResponse(seq.format('genbank'),
                                content_type='text/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}.gb"'.format(seq.id)
        return response
    return _fungal_domestication_view_no_template(request, function)

@login_required
def fungal_domestication_view_add(request):
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST
    # request_data
    seq = request_data['seq']
    category = request_data['category']
    prefix = request_data['prefix']
    suffix = request_data['suffix']
    seq_name = request_data['seq_name']
    name = request_data['name']
    # Enzymes not working ???
    try:
        enzymes = request_data['enzymes']
        # it returns a unicode string, converting to list
        if isinstance(enzymes, unicode):
            enzymes = _parse_unicode_into_list(enzymes)

    except:
        enzymes = None

    with_noncoding_seq = bool(int(request_data['with_noncoding_seq']))
    if category == 'None':
        category_name = 'Other'
    else:
        category_name = FUNGAL_CATEGORIES[category][0]
    if category_name == REVERSE_MARKER:
        reverse_orientation = True
    else:
        reverse_orientation = False
    seq = SeqRecord(Seq(seq), id=seq_name, name=seq_name)
    seq = domesticate(seq, category, prefix, suffix, enzymes=enzymes,
                      with_intron=with_noncoding_seq, reverse_orientation=reverse_orientation)[1]
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(seq.id),
                                    suffix='.gb', dir=TMP_DIR)
    temp_fhand.write(seq.format('gb').encode())
    temp_fhand.flush()
    temp_fhand.seek(0)
    props = {'Description': [textwrap.fill(request_data['description'], width=70)],
             'Reference': [request_data['reference']]}
    try:
        feature = add_feature(name=name, type_name=category_name,
                              vector=DOMESTICATED_VECTOR, genbank=temp_fhand,
                              props=props, owner=request.user, is_public=False)

    except IntegrityError as error:
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render_to_response('feature_exists.html',
                                      {},
                                      context_instance=RequestContext(request))
        else:
            return HttpResponseServerError()
    except Exception as error:
        return HttpResponseServerError()
    # if everithing os fine we show the just added feature
    return redirect(feature.url)

def fungal_domestication_view_protocol(request):
    def function(seq, category, prefix, suffix, enzymes, with_noncoding_seq,
                 reverse_orientation):
        pcrs, seq = domesticate(seq, category, prefix, suffix, enzymes,
                                with_intron=with_noncoding_seq, reverse_orientation=reverse_orientation)
                                
        protocol = write_domestication_protocol(pcrs)
        response = HttpResponse(protocol, content_type='text/plain')
        response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
        return response

    return _fungal_domestication_view_no_template(request, function)

def fungal_domestication_view_sbol(request):
    def function(seq, category, prefix, suffix,  enzymes, with_noncoding_seq, reverse_orientation):
        seq = domesticate(seq, category, prefix, suffix, enzymes,
                          with_intron=with_noncoding_seq, reverse_orientation=reverse_orientation)[1]
        response = HttpResponse(convert_to_sbol(seq),
                                content_type='text/xml')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}.xml"'.format(seq.id)
        return response
    return _fungal_domestication_view_no_template(request, function)


def synthesis_view_genbank(request):
    def function(record, category, prefix, suffix, enzymes, with_intron):
        parsed_record = SeqIO.read(StringIO(record), "genbank")
        response = HttpResponse(record.format('genbank'),
                                content_type='text/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}.gb"'.format(parsed_record.id)
        return response
    return _domestication_view_genbank_no_template(request, function)


def synthesis_view_sbol(request):
    def function(seq, category, prefix, suffix, enzymes, with_intron):
        seq = domesticate_for_synthesis(seq, category, prefix, suffix, enzymes,
                                        with_intron)[1]
        response = HttpResponse(convert_to_sbol(seq),
                                content_type='xml/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}.xml"'.format(seq.id)
        return response
    return _domestication_view_no_template(request, function)


def domestication_view_genbank(request):
    def function(seq, category, prefix, suffix, enzymes, with_intron):
        seq = domesticate(seq, category, prefix, suffix, enzymes,
                          with_intron=with_intron)[1]
        response = HttpResponse(seq.format('genbank'),
                                content_type='text/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}.gb"'.format(seq.id)
        return response
    return _domestication_view_no_template(request, function)


def  domestication_view_regular_genbank(request):
    def function(record, category, prefix, suffix, enzymes, with_intron):
        #seq = domesticate(seq, category, prefix, suffix, enzymes,
                           #with_intron=with_intron)[1]
        
        parsed_record = SeqIO.read(StringIO(record), "genbank")
        response = HttpResponse(record.format('genbank'),
                                content_type='text/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}.gb"'.format(parsed_record.id)
        return response
    return _domestication_view_genbank_no_template(request, function)


def domestication_view_sbol(request):
    def function(seq, category, prefix, suffix, enzymes, with_intron):
        seq = domesticate(seq, category, prefix, suffix, enzymes,
                          with_intron=with_intron)[1]
        response = HttpResponse(convert_to_sbol(seq),
                                content_type='text/xml')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}.xml"'.format(seq.id)
        return response
    return _domestication_view_no_template(request, function)


def domestication_view_protocol(request):
    def function(seq, category, prefix, suffix, enzymes, with_intron):
        pcrs, seq = domesticate(seq, category, prefix, suffix, enzymes,
                                with_intron=with_intron)
        protocol = write_domestication_protocol(pcrs)
        response = HttpResponse(protocol, content_type='text/plain')
        response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
        return response

    return _domestication_view_no_template(request, function)


def write_synthesis_protocol(seq_for_syn):
    protocol = '''Domestication Protocol for the Synthetic Strategy

Order the following sequence:
{0}

Once you have the synthetic sequence the domestication reaction should be performed as follows:
40 ng of the synthesized product
75 ng of pUPD2
5-10u BsmBI
3u T4 Ligase
1 microlitre Ligase Buffer

Final volume: 10 microlitres

We use Promega T4 DNA ligase(M180B) and fermentas BsmBI/Esp3I (ER0451). We haven't tried other enzymes suppliers but they will problably work as well.

Set your reaction in a thermocycler: 25 cycles x (37C 2', 16C 5').
One microlitre of the reaction is enough to be transform E.coli electrocompetent cells. Positive clones are selected in Chloramphenicol (25 microgram ml-1), IPTG (0.5mM) and Xgal (40 microgram ml-1) plates. You will distinguish between colonies carrying intact vectors (blue) and those transformed with your construction (white).

'''
    return protocol.format(fill(seq_for_syn, 80))

def write_cas12a_protocol(seq_for_syn):
    protocol = '''Domestication Protocol for the Cas12a Multiplex Strategy

Order the following sequence:
{0}

Once you have the synthetic sequence the domestication reaction should be performed as follows:
40 ng of the synthesized product
75 ng of pUPD2
5-10u BsmBI
3u T4 Ligase
1 microlitre Ligase Buffer

Final volume: 10 microlitres

We use Promega T4 DNA ligase(M180B) and fermentas BsmBI/Esp3I (ER0451). We haven't tried other enzymes suppliers but they will problably work as well.

Set your reaction in a thermocycler: 25 cycles x (37C 2', 16C 5').
One microlitre of the reaction is enough to be transform E.coli electrocompetent cells. Positive clones are selected in Chloramphenicol (25 microgram ml-1), IPTG (0.5mM) and Xgal (40 microgram ml-1) plates. You will distinguish between colonies carrying intact vectors (blue) and those transformed with your construction (white).

'''
    return protocol.format(fill(seq_for_syn, 80))


def write_domestication_protocol(pcrs):
    protocol = '''Domestication Protocol

Perform a PCR amplification for each patch with the given pair of oligos by using your DNA Polymerase manufacturer's protocol:

{0}


Once you have all your patches the domestication reaction should be performed as follows:
40 ng of each patch
75 ng of pUPD2
5-10u BsmBI
3u T4 Ligase
1 microlitre Ligase Buffer

Final volume: 10 microlitres

We use Promega T4 DNA ligase(M180B) and fermentas BsmBI/Esp3I (ER0451). We haven't tried other enzymes suppliers but they will problably work as well.

Set your reaction in a thermocycler: 25 cycles x (37C 2', 16C 5').
One microlitre of the reaction is enough to be transform E.coli electrocompetent cells. Positive clones are selected in Chloramphenicol (25 microgram ml-1), IPTG (0.5mM) and Xgal (40 microgram ml-1) plates. You will distinguish between colonies carrying intact vectors (blue) and those transformed with your construction (white).
'''
    pcr_str = ''
    for pcr in pcrs:
        pcr_str += '\tPCR product: {0}\n'.format(pcr['pcr_product'])
        pcr_str += '\tOligo forward: {0}\n'.format(pcr['oligo_forward'])
        pcr_str += '\tOligo reverse: {0}\n'.format(pcr['oligo_reverse'])
        pcr_str += '\n'

    return protocol.format(pcr_str)


@login_required
def crispr_view_add(request):
    context = {}
    context.update(csrf(request))
    request_data = request.POST
    # request_data
    seq = request_data['seq']
    category = request_data['category']
    prefix = request_data['prefix']
    suffix = request_data['suffix']
    seq_name = request_data['seq_name']
    name = request_data['name']
    if category == 'None':
        category_name = 'Other'
    else:
        category_name = CRYSPER_CATEGORIES[category][0]
    
    seq = Seq(seq, alphabet=generic_dna)
    annotated_seq = annotate_domesticated_crispr(seq, seq_name, prefix, 
                                                 suffix, category)
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(annotated_seq.id),
                                    suffix='.gb', dir=TMP_DIR)
    temp_fhand.write(annotated_seq.format('gb').encode())
    temp_fhand.flush()
    temp_fhand.seek(0)
    props = {'Description': [textwrap.fill(request_data['description'], width=70)],
             'Reference': [request_data['reference']]}
    try:
        # prefix has an added A nt for view comfort, but prefix must be 4 nts long
        if prefix == "GTGCA":
            prefix = "GTGC"
        feature = add_feature(name=name, type_name=category_name,
                              vector=None, genbank=temp_fhand, props=props,
                              owner=request.user, is_public=False,
                              prefix=prefix, suffix=suffix)

    except IntegrityError as error:
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render(request, 'feature_exists.html', context=context)
        else:
            return HttpResponseServerError()
    except Exception as error:
        raise
        return HttpResponseServerError(error)
    # if everithing os fine we show the just added feature
    return redirect(feature.url)


def _parse_unicode_into_list(string):
    # template return list of enzymes as unicode
    # this convert that in a list
    string = string.encode("ascii")
    list = [item for item in string.strip("[]").split(",")]
    list = [item.strip(" u'") for item in list]
    return list


@login_required
def domestication_view_add(request):
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST
    # request_data
    seq = request_data['seq']
    category = request_data['category']
    prefix = request_data['prefix']
    suffix = request_data['suffix']
    seq_name = request_data['seq_name']
    name = request_data['name']
    # Enzymes not working ???
    try:
        enzymes = request_data['enzymes']
        # it returns a unicode string, converting to list
        if isinstance(enzymes, unicode):
            enzymes = _parse_unicode_into_list(enzymes)

    except:
        enzymes = None

    with_intron = bool(int(request_data['with_intron']))
    if category == 'None':
        category_name = 'Other'
    else:
        category_name = CATEGORIES[category][0]
    seq = SeqRecord(Seq(seq), id=seq_name, name=seq_name)
    if 'seq_features' in request_data:
        merged_features = []
        features = eval(request_data['seq_features'])
        qualifiers = eval(request_data['features_qualifiers'])
        for feature in features:
            feature.qualifiers = qualifiers.pop(0)
            merged_features.append(feature)
            seq = domesticate(seq, category, prefix, suffix, enzymes=enzymes,
                              with_intron=with_intron, features=merged_features)[1]
    else:
        seq = domesticate(seq, category, prefix, suffix, enzymes=enzymes,
                          with_intron=with_intron)[1]
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(seq.id),
                                    suffix='.gb', dir=TMP_DIR)
    temp_fhand.write(seq.format('gb').encode())
    temp_fhand.flush()
    temp_fhand.seek(0)
    props = {'Description': [textwrap.fill(request_data['description'], width=70)],
             'Reference': [request_data['reference']]}
    try:
        feature = add_feature(name=name, type_name=category_name,
                              vector=DOMESTICATED_VECTOR, genbank=temp_fhand,
                              props=props, owner=request.user, is_public=False)

    except IntegrityError as error:
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render_to_response('feature_exists.html',
                                      {},
                                      context_instance=RequestContext(request))
        else:
            return HttpResponseServerError()
    except Exception as error:
        return HttpResponseServerError()
    # if everithing os fine we show the just added feature
    return redirect(feature.url)



def _get_features(record):
    parsed_record = SeqIO.read(StringIO(record), "genbank")
    return parsed_record.features


@login_required
def domestication_view_regular_add(request):
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST
    seq = request_data['seq']
    features = _get_features(request_data['record'])
    category = request_data['category']
    prefix = request_data['prefix']
    suffix = request_data['suffix']
    seq_name = request_data['seq_name']
    name = request_data['name']
    # Enzymes not working ???
    try:
        enzymes = request_data['enzymes']
        # it returns a unicode string, converting to list
        if isinstance(enzymes, unicode):
            enzymes = _parse_unicode_into_list(enzymes)

    except:
        enzymes = None

    with_intron = bool(int(request_data['with_intron']))
    if category == 'None':
        category_name = 'Other'
    else:
        category_name = CATEGORIES[category][0]
    seq = SeqRecord(Seq(seq), id=seq_name, name=seq_name)
    seq = domesticate(seq, category, prefix, suffix, enzymes=enzymes,
                          with_intron=with_intron)[1]
    for feature in features:
        seq.features.append(feature)
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(seq.id),
                                    suffix='.gb', dir=TMP_DIR)
    temp_fhand.write(seq.format('gb').encode())
    temp_fhand.flush()
    temp_fhand.seek(0)
    props = {'Description': [textwrap.fill(request_data['description'], width=70)],
             'Reference': [request_data['reference']]}
    try:
        feature = add_feature(name=name, type_name=category_name,
                              vector=DOMESTICATED_VECTOR, genbank=temp_fhand,
                              props=props, owner=request.user, is_public=False)

    except IntegrityError as error:
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render_to_response('feature_exists.html',
                                      {},
                                      context_instance=RequestContext(request))
        else:
            return HttpResponseServerError()
    except Exception as error:
        return HttpResponseServerError()
    # if everithing os fine we show the just added feature
    return redirect(feature.url)


CRYSPER_PROTOCOL = '''In order to use this sequence as a target in a guide RNA you have to order these primers:
Forward primer: 5' {forw_primer} 3'
Reverse primer: 5' {rev_primer} 3'
Dilute them to a final concentration of 1uM, mix 5ul of each primer and let them anneal for at least 30 minutes at room temperature before setting up the multipartite reaction.

'''


def crispr_view_protocol(request):
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST
    forw_primer = request_data['forw_primer']
    rev_primer = request_data['rev_primer']
    protocol = CRYSPER_PROTOCOL.format(forw_primer=forw_primer,
                                       rev_primer=rev_primer)
    response = HttpResponse(protocol, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
    return response


def synthesis_view_protocol(request):
    def function(seq, category, prefix, suffix, enzymes, with_intron):
        seq = domesticate_for_synthesis(seq, category, prefix, suffix, enzymes,
                                        with_intron)[0]
        protocol = write_synthesis_protocol(seq)
        response = HttpResponse(protocol, content_type='text/plain')
        response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
        return response

    return _domestication_view_no_template(request, function)


def cas12a_multiplexing_view_protocol(request):
    def function(targets, category):
        seq, _ = assemble_cas12_unit_zero(targets, category)
        protocol = write_cas12a_protocol(str(seq.seq))
        response = HttpResponse(protocol, content_type='text/plain')
        response['Content-Disposition'] = 'attachment; filename="protocol.txt"'
        return response

    return _cas12a_multiplex_domestication_view_no_template(request, function)


def _cas12a_multiplex_domestication_view_no_template(request, function):
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST
    targets = eval(request_data['targets'])
    category = request_data['category']
    #seq_name = request_data['seq_name']
    return function(targets, category)

def _domestication_view_no_template(request, function):
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST
    
    category = request_data['category']
    prefix = request_data['prefix']
    suffix = request_data['suffix']
    enzymes = json.loads(request_data['enzymes'])
    with_intron = request_data['with_intron']
    with_intron = bool(int(with_intron))
    seq_name = request_data['seq_name']
    seq = request_data['seq']
    seq = SeqRecord(Seq(seq), id=seq_name, name=seq_name)
    if 'seq_features' in request_data:
        features = eval(request_data['seq_features'])
        qualifiers = eval(request_data['features_qualifiers'])
        for feature in features:
            feature.qualifiers = qualifiers.pop(0)
            seq.features.append(feature)
        
    return function(seq, category, prefix, suffix, enzymes, with_intron)


def _domestication_view_genbank_no_template(request, function):
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST
    category = request_data['category']
    prefix = request_data['prefix']
    suffix = request_data['suffix']
    enzymes = json.loads(request_data['enzymes'])
    with_intron = request_data['with_intron']
    with_intron = bool(int(with_intron))
    seq_name = request_data['seq_name']
    record = request_data['record']
    return function(record, category, prefix, suffix, enzymes, with_intron)


def _fungal_domestication_view_no_template(request, function):
    if not request.POST:
        msg = "To show the protocol you need first to assemble parts"
        return HttpResponseBadRequest(msg)
    context = RequestContext(request)
    context.update(csrf(request))
    request_data = request.POST
    seq = request_data['seq']
    category = request_data['category']
    prefix = request_data['prefix']
    suffix = request_data['suffix']
    enzymes = json.loads(request_data['enzymes'])
    with_intron = request_data['with_noncoding_seq']
    if category == REVERSE_MARKER:
        reverse_orientation = True
    else:
        reverse_orientation = False
    with_intron = bool(int(with_intron))
    seq_name = request_data['seq_name']
    seq = SeqRecord(Seq(seq), id=seq_name, name=seq_name)
    return function(seq, category, prefix, suffix, enzymes, with_intron,
                    reverse_orientation)
