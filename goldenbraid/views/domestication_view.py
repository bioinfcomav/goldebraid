from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from django.core.context_processors import csrf
from django.template.context import RequestContext
from django.shortcuts import render_to_response
from django.http import HttpResponse, HttpResponseBadRequest

from goldenbraid.domestication import domesticate
from goldenbraid.forms import DomesticationForm
from goldenbraid.settings import CATEGORIES


def domestication_view(request):
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
        if form.is_valid():
            # do domestication
            seq = form.cleaned_data['seq']
            category = form.cleaned_data.get('category', None)
            if category is None:
                prefix = form.cleaned_data.get('prefix')
                suffix = form.cleaned_data.get('suffix')
            else:
                prefix = CATEGORIES[category][1]
                suffix = CATEGORIES[category][2]
            pcr = domesticate(seq, category, prefix, suffix)[0]
            return render_to_response('domestication_result.html',
                                      {'category': category,
                                       'prefix': prefix,
                                       'suffix': suffix,
                                       'pcrs': pcr,
                                       'seq': str(seq.seq),
                                       'seq_name': seq.name},
                                context_instance=RequestContext(request))
    else:
        form = DomesticationForm()
    context['form'] = form

    template = 'domestication_template.html'
    mimetype = None
    return render_to_response(template, context, mimetype=mimetype)


def domestication_view_genbank(request):
    def function(pcrs, seq):
        response = HttpResponse(seq.format('genbank'),
                                mimetype='text/plain')
        response['Content-Disposition'] = 'attachment; '
        response['Content-Disposition'] += 'filename="{}.gb"'.format(seq.id)
        return response
    return _domestication_view_no_template(request, function)


def write_domestication_protocol(pcrs):
    protocol = '''Domestication Protocol

Perform a PCR amplification for each patch with the given pair of oligos by using your DNA Polymerase manufacturer's protocol:

{0}


Once you have all your patches the domestication reaction should be performed as follows:
40 ng of each patch
75 ng of pUPD
3u BsmBI
3u T4 Ligase
1 microlitre Ligase Buffer

Final volume: 10 microlitres

Set your reaction in a thermocycler: 25 cycles x (37C 2', 16C 5').
One microlitre of the reaction is enough to be transform E.coli electrocompetent cells. Positive clones are selected in Ampicillin (50 microgram ml-1), IPTG (0.5mM) and Xgal (40 microgram ml-1) plates. You will distinguish between colonies carrying intact vectors (blue) and those transformed with your construction (white).
'''
    pcr_str = ''
    for pcr in pcrs:
        pcr_str += '\tPCR product: {0}\n'.format(pcr['pcr_product'])
        pcr_str += '\tOligo forward: {0}\n'.format(pcr['oligo_forward'])
        pcr_str += '\tOligo reverse: {0}\n'.format(pcr['oligo_reverse'])
        pcr_str += '\n'

    return protocol.format(pcr_str)


def domestication_view_protocol(request):
    def function(pcrs, seq):
        protocol = write_domestication_protocol(pcrs)
        response = HttpResponse(protocol, mimetype='text/plain')
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
    seq_name = request_data['seq_name']
    seq = SeqRecord(Seq(seq), id=seq_name, name=seq_name)
    pcrs, seq = domesticate(seq, category, prefix, suffix)
    return function(pcrs, seq)
