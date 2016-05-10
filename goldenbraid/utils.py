import re

from django.db.models import Q

from Bio.Seq import Seq

from goldenbraid.settings import REBASE_FILE, MANDATORY_DOMEST_ENZYMES


def filter_feature_by_user_perms(query, user):
    if user.is_staff:
        return query
    if user.is_authenticated():
        query = query.filter(Q(featureperm__owner__username=user) |
                             Q(featureperm__is_public=True))
    else:
        query = query.filter(featureperm__is_public=True)

    return query

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


def get_ret_sites(enzymes):
    'It returns the restriction site of the given enzymes'
    sites = []
    existing_enzymes = parse_rebase_file(REBASE_FILE)
    for enzyme in enzymes:
        site = existing_enzymes[enzyme]
        if '^' in site:
            raise ValueError('We only use type IIS enzymes')
        site = site.split('(')[0]
        rev_site = str(Seq(site).reverse_complement())
        sites.extend([site, rev_site])
    return sites


def has_rec_sites(seq, enzymes=None):
    if enzymes is None:
        enzymes = MANDATORY_DOMEST_ENZYMES
    rec_sites = get_ret_sites(enzymes)
    # regex with the sites to domesticate
    rec_sites_regex = '(' + '|'.join(rec_sites) + ')'
    rec_sites_regex = re.compile(rec_sites_regex, flags=re.IGNORECASE)
    match = rec_sites_regex.search(seq)
    return True if match else False


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


def _pref_suf_index_from_rec_sites(seq, forw_site, rev_site, rec_site,
                                   forw_cut_delta, rev_cut_delta):

    prefix_index = forw_site + len(rec_site) + forw_cut_delta
    if prefix_index >= len(seq):
        prefix_index = prefix_index - len(seq)

    suffix_index = rev_site - rev_cut_delta
    if suffix_index < 0:
        suffix_index = len(seq) - abs(suffix_index)
    return prefix_index, suffix_index


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
