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

from tempfile import NamedTemporaryFile
from textwrap import fill
import json

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet

from django.core.context_processors import csrf
from django.template.context import RequestContext
from django.shortcuts import render_to_response, redirect
from django.http import HttpResponse, HttpResponseBadRequest
from django.contrib.auth.decorators import login_required
from django.http.response import HttpResponseServerError
from django.db.utils import IntegrityError

from goldenbraid.domestication import (domesticate, domesticate_for_synthesis,
                                       domestication_crispr)
from goldenbraid.settings import (CATEGORIES, CRYSPER_CATEGORIES,
                                  DOMESTICATED_VECTOR)

from goldenbraid.forms.domestication import (DomesticationForm,
                                             DomesticationCrisprForm)
from goldenbraid.views.feature import add_feature
from goldenbraid.sbol import convert_to_sbol


def synthesis_view(request):
    return _domestication_view(request, kind='synthesis')


def domestication_view(request):
    return _domestication_view(request, kind='domestication')


def crispr_view(request):
    context = RequestContext(request)
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

            return render_to_response('crispr_result.html',
                                      {'category': category,
                                       'prefix': prefix,
                                       'suffix': suffix,
                                       'seq': str(new_seq.seq),
                                       'seq_name': new_seq.name,
                                       'forw_dim': str(forw_dim.seq),
                                       'rev_dim': str(rev_dim.seq)},
                                  context_instance=RequestContext(request))
    else:
        form = DomesticationCrisprForm()

    context['form'] = form

    template = 'crispr_template.html'
    content_type = None
    return render_to_response(template, context, content_type=content_type)


def _domestication_view(request, kind):
    context = RequestContext(request)
    context.update(csrf(request))
    if request.method == 'POST':
        request_data = request.POST
    elif request.method == 'GET':
        request_data = request.GET
    else:
        request_data = None
    if request_data:
        form = DomesticationForm(request_data, request.FILES)
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
                    pcr = domesticate(seq, category, prefix, suffix, enzymes,
                                      with_intron)[0]
                except RuntimeError as error:
                    return render_to_response('goldenbraid_info.html',
                                              {'title': 'Can not domesticate sequence',
                                               'info': error},
                                              context_instance=RequestContext(request))

                return render_to_response('domestication_result.html',
                                          {'category': category,
                                           'prefix': prefix,
                                           'suffix': suffix,
                                           'pcrs': pcr,
                                           'seq': str(seq.seq),
                                           'seq_name': seq.name,
                                           'enzymes': enzymes,
                                           'with_intron': with_intron_str},
                                      context_instance=RequestContext(request))
            elif kind == 'synthesis':
                try:
                    seq_for_syn, prepared_seq = domesticate_for_synthesis(seq,
                                                                          category,
                                                                          prefix,
                                                                          suffix,
                                                                          enzymes,
                                                                          with_intron)
                except RuntimeError as error:
                    return render_to_response('goldenbraid_info.html',
                                              {'title': 'Can not domesticate sequence',
                                               'info': error},
                                              context_instance=RequestContext(request))
                return render_to_response('synthesis_result.html',
                                          {'category': category,
                                           'prefix': prefix,
                                           'suffix': suffix,
                                           'seq_syn': seq_for_syn,
                                           'seq': str(seq.seq),
                                           'seq_name': prepared_seq.name,
                                           'enzymes': enzymes,
                                           'with_intron': with_intron_str},
                                      context_instance=RequestContext(request))
    else:
        form = DomesticationForm()
    context['form'] = form
    context['kind'] = kind

    template = 'domestication_template.html'
    content_type = None
    return render_to_response(template, context, content_type=content_type)


def synthesis_view_genbank(request):
    def function(seq, category, prefix, suffix, enzymes, with_intron):
        seq = domesticate_for_synthesis(seq, category, prefix, suffix, enzymes,
                                        with_intron)[1]
        response = HttpResponse(seq.format('genbank'),
                                content_type='text/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}.gb"'.format(seq.id)
        return response
    return _domestication_view_no_template(request, function)


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
                          with_intron)[1]
        response = HttpResponse(seq.format('genbank'),
                                content_type='text/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}.gb"'.format(seq.id)
        return response
    return _domestication_view_no_template(request, function)


def domestication_view_sbol(request):
    def function(seq, category, prefix, suffix, enzymes, with_intron):
        seq = domesticate(seq, category, prefix, suffix, enzymes,
                          with_intron)[1]
        response = HttpResponse(convert_to_sbol(seq),
                                content_type='text/xml')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{0}.xml"'.format(seq.id)
        return response
    return _domestication_view_no_template(request, function)


def domestication_view_protocol(request):
    def function(seq, category, prefix, suffix, enzymes, with_intron):
        pcrs, seq = domesticate(seq, category, prefix, suffix, enzymes,
                                with_intron)
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
    if category == 'None':
        category_name = 'Other'
    else:
        category_name = CRYSPER_CATEGORIES[category][0]
    print 'seq', seq
    seq = SeqRecord(Seq(seq, DNAAlphabet()), id=seq_name, name=seq_name)
    print 'seq', seq
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(seq.id),
                                    suffix='.gb')
    temp_fhand.write(seq.format('gb'))
    temp_fhand.flush()
    temp_fhand.seek(0)
    props = {'Description': [request_data['description']],
             'Reference': [request_data['reference']]}
    try:
        feature = add_feature(name=name, type_name=category_name,
                              vector=None, genbank=temp_fhand, props=props,
                              owner=request.user, is_public=False,
                              prefix=prefix, suffix=suffix)

    except IntegrityError as error:
        print error
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render_to_response('feature_exists.html',
                                      {},
                                      context_instance=RequestContext(request))
        else:
            return HttpResponseServerError()
    except Exception as error:
        raise
        print error
        return HttpResponseServerError(error)
    # if everithing os fine we show the just added feature
    return redirect(feature.url)


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
    with_intron = bool(int(request_data['with_intron']))
    if category == 'None':
        category_name = 'Other'
    else:
        category_name = CATEGORIES[category][0]
    seq = SeqRecord(Seq(seq), id=seq_name, name=seq_name)
    seq = domesticate(seq, category, prefix, suffix, with_intron)[1]
    temp_fhand = NamedTemporaryFile(prefix='{0}.'.format(seq.id),
                                    suffix='.gb')
    temp_fhand.write(seq.format('gb'))
    temp_fhand.flush()
    temp_fhand.seek(0)
    props = {'Description': [request_data['description']],
             'Reference': [request_data['reference']]}
    try:
        feature = add_feature(name=name, type_name=category_name,
                              vector=DOMESTICATED_VECTOR, genbank=temp_fhand,
                              props=props, owner=request.user, is_public=False)

    except IntegrityError as error:
        print error
        if 'feature already in db' in str(error):
            # TODO choose a template
            return render_to_response('feature_exists.html',
                                      {},
                                      context_instance=RequestContext(request))
        else:
            return HttpResponseServerError()
    except Exception as error:
        print error
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




def _domestication_view_no_template(request, function):
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
    with_intron = request_data['with_intron']
    with_intron = bool(int(with_intron))
    seq_name = request_data['seq_name']
    seq = SeqRecord(Seq(seq), id=seq_name, name=seq_name)
    return function(seq, category, prefix, suffix, enzymes, with_intron)
