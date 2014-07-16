import re

from django.db import models

from gb_genome_domestication.settings import DB, ENZYMES_USED_IN_GOLDENBRAID

from goldenbraid.domestication import get_ret_sites


class Db(models.Model):
    db_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=255)
    description = models.CharField(max_length=255)
    urlprefix = models.CharField(max_length=255)
    url = models.CharField(max_length=255)

    class Meta:
        db_table = u'db'


class Dbxref(models.Model):
    dbxref_id = models.AutoField(primary_key=True)
    db = models.ForeignKey(Db)
    accession = models.CharField(max_length=255)

    class Meta:
        db_table = u'dbxref'

    @property
    def url(self):
        return  self.db.urlprefix + self.accession


class Feature(models.Model):
    feature_id = models.AutoField(primary_key=True)
    uniquename = models.CharField(max_length=255, unique=True)
    name = models.CharField(max_length=255)
    description = models.CharField(max_length=1024)
    species = models.CharField(max_length=255)
    residues = models.TextField()
    dbxref = models.ForeignKey(Dbxref)

    class Meta:
            db_table = u'feature'

    @property
    def primary_dbxref(self):
        return self.dbxref

    @property
    def num_rec_sites(self):
        rec_sites = get_ret_sites(ENZYMES_USED_IN_GOLDENBRAID)
        # regex with the sites to domesticate
        rec_sites_regex = '(' + '|'.join(rec_sites) + ')'
        rec_sites_regex = re.compile(rec_sites_regex)
        return len(re.findall(rec_sites_regex, self.residues))

    @property
    def secondary_dbxrefs(self):
        for feat_dbxref in FeatureDbxref.objects.using(DB).filter(feature=self):
            yield feat_dbxref.dbxref


class FeatureDbxref(models.Model):
    feature_dbxref_id = models.AutoField(primary_key=True)
    feature = models.ForeignKey(Feature)
    dbxref = models.ForeignKey(Dbxref)

    class Meta:
        db_table = u'feature_dbxref'

