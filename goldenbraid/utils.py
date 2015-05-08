'''
Created on 2015 mai. 6

@author: peio
'''
from goldenbraid.settings import REBASE_FILE, ENZYMES_USED_IN_GOLDENBRAID
import re
from Bio.Seq import Seq


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


def has_rec_sites(seq):
    rec_sites = get_ret_sites(ENZYMES_USED_IN_GOLDENBRAID)
    # regex with the sites to domesticate
    rec_sites_regex = '(' + '|'.join(rec_sites) + ')'
    rec_sites_regex = re.compile(rec_sites_regex, flags=re.IGNORECASE)
    match = rec_sites_regex.search(seq)
    return True if match else False
