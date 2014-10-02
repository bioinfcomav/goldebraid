'''
Created on 2014 uzt 1

@author: peio
'''

from django.db import transaction
from django.core.management.base import BaseCommand, CommandError

from Bio import SeqIO
from gb_genome_domestication.models import Db, Dbxref, Feature, FeatureDbxref
from gb_genome_domestication.management.commands.db_utils import get_or_load_db
from gb_genome_domestication import settings


DB = settings.DB
DB_URLPREFIX = settings.DB_URLPREFIX


class Command(BaseCommand):
    args = '<fasta_file> <spp> <database name> <database_urlprefix>'
    help = 'Adds the given cdss to the database'

    def handle(self, *args, **options):
        'Adds the given cdss to the database'
        if not args:
            raise CommandError('No file given')
        else:
            cds_fpath = args[0]
            spp = args[1]
            database_name = args[2]
            urlprefix = args[3]
        try:
            run_command(cds_fpath, spp, database_name, urlprefix)
        except Exception as error:
            raise CommandError(str(error))


def run_command(cds_fpath, spp, database_name, urlprefix):
    with transaction.commit_on_success():
        if not database_name or not urlprefix:
            ext_db = None
        else:
            ext_db = get_or_load_db(database_name, urlprefix=urlprefix)
        prim_db = get_or_load_db(DB, urlprefix=DB_URLPREFIX)

        for seq in SeqIO.parse(open(cds_fpath), format='fasta'):
            name = seq.name
            id_ = seq.id
            desc = seq.description
            residues = str(seq.seq)
            prim_dbxref = Dbxref.objects.using(DB).create(db=prim_db,
                                                          accession=id_)
            feat = Feature.objects.using(DB).create(uniquename=id_,
                                                    name=name, species=spp,
                                                    description=desc,
                                                    dbxref=prim_dbxref,
                                                    residues=residues)
            if ext_db is not None:
                ext_dbxref = Dbxref.objects.using(DB).create(db=ext_db,
                                                             accession=id_)
                FeatureDbxref.objects.using(DB).create(feature=feat,
                                                       dbxref=ext_dbxref)
