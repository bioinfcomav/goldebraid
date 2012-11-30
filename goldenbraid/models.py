from django.db import models
from goldenbraid import settings
from goldenbraid.tags import DESCRIPTION_TYPE_NAME

DB = settings.DB


class Db(models.Model):
    db_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=255)
    description = models.CharField(max_length=255)
    urlprefix = models.CharField(max_length=255)
    url = models.CharField(max_length=255)


class Cv(models.Model):
    cv_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=255)
    definition = models.TextField()


class Dbxref(models.Model):
    dbxref_id = models.AutoField(primary_key=True)
    db = models.ForeignKey(Db)
    accession = models.CharField(max_length=255)


class Cvterm(models.Model):
    cvterm_id = models.AutoField(primary_key=True)
    cv = models.ForeignKey(Cv)
    name = models.CharField(max_length=1024)
    definition = models.TextField()


class Contact(models.Model):
    contact_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    email = models.EmailField(unique=True)


class Stockcollection(models.Model):
    stockcollection_id = models.AutoField(primary_key=True)
    # type = models.ForeignKey(Cvterm)
    contact = models.ForeignKey(Contact)
    name = models.CharField(max_length=255)
    uniquename = models.TextField()


class Stock(models.Model):
    stock_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    uniquename = models.TextField()
    description = models.TextField()
    stockcollection = models.ForeignKey(Stockcollection)
    feature = models.ForeignKey("Feature", related_name='stocks',
                                null=True)


class ReadOnlyDict(dict):
    'a dict with disabled setitem'
    def __init__(self, *args, **kwargs):
        'initiate as a regular dict'
        super(ReadOnlyDict, self).__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        'disabled setitem'
        raise ValueError('OnlyRead dictionary')


class Feature(models.Model):
    "The feature model"
    feature_id = models.AutoField(primary_key=True)
    uniquename = models.CharField(max_length=255, unique=True)
    name = models.CharField(max_length=255)
    type = models.ForeignKey(Cvterm)
    residues = models.TextField()
    dbxref = models.ForeignKey(Dbxref)
    vector = models.ForeignKey("Feature", null=True)
    timecreation = models.DateTimeField(auto_now_add=True)
    timelastmodified = models.DateTimeField(auto_now=True)

    @property
    def seq_len(self):
        return len(self.residues)

    @property
    def props(self):
        'It returns a list of cvterms'
        new_prop_dict = {}
        props = Featureprop.objects.using(DB).filter(feature=self)
        for prop in props:
            type_ = prop.type.name
            value = prop.value
            rank = prop.rank
            if type_ not in new_prop_dict:
                new_prop_dict[type_] = []
            new_prop_dict[type_].append((value, rank))
        prop_dict = {}
        for type_, values in new_prop_dict.items():

            values = sorted(values, key=lambda x: x[1])
            values = [value[0] for value in values]
            prop_dict[type_] = values

        return ReadOnlyDict(prop_dict)

    @property
    def description(self):
        'Get descrition if it has one'
        try:
            return self.props[DESCRIPTION_TYPE_NAME][0]
        except KeyError:
            return None


class Featureprop(models.Model):
    'Model to store the properties of the features'
    featureprop_id = models.AutoField(primary_key=True)
    feature = models.ForeignKey(Feature)
    type = models.ForeignKey(Cvterm)
    value = models.TextField()
    rank = models.IntegerField()

    class meta:
        unique_together = ('feature', 'type', 'rank')

