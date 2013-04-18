from django.db import models
from django.contrib.auth.models import User
from goldenbraid import settings
from goldenbraid.tags import (DESCRIPTION_TYPE_NAME, ENZYME_IN_TYPE_NAME,
                              VECTOR_TYPE_NAME, ENZYME_OUT_TYPE_NAME,
                              RESISTANCE_TYPE_NAME, REFERENCE_TYPE_NAME,
                              FORWARD, REVERSE)


class Db(models.Model):
    db_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=255)
    description = models.CharField(max_length=255)
    urlprefix = models.CharField(max_length=255)
    url = models.CharField(max_length=255)

    class Meta:
        db_table = u'db'


class Cv(models.Model):
    cv_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=255)
    definition = models.TextField()

    class Meta:
        db_table = u'cv'


class Dbxref(models.Model):
    dbxref_id = models.AutoField(primary_key=True)
    db = models.ForeignKey(Db)
    accession = models.CharField(max_length=255)

    class Meta:
        db_table = u'dbxref'

    @property
    def url(self):
        return  self.db.urlprefix + self.accession


class Cvterm(models.Model):
    cvterm_id = models.AutoField(primary_key=True)
    cv = models.ForeignKey(Cv)
    name = models.CharField(max_length=1024)
    definition = models.TextField()

    class Meta:
        db_table = u'cvterm'


class Count(models.Model):
    count_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    value = models.IntegerField()

    class Meta:
        db_table = u'count'

    @property
    def next(self):
        next_ = self.value
        self.value += 1
        self.save()
        return hex(next_)[2:].upper()


class Contact(models.Model):
    contact_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    email = models.EmailField(unique=True)

    class Meta:
        db_table = u'contact'


class Stockcollection(models.Model):
    stockcollection_id = models.AutoField(primary_key=True)
    # type = models.ForeignKey(Cvterm)
    contact = models.ForeignKey(Contact)
    name = models.CharField(max_length=255)
    uniquename = models.TextField()

    class Meta:
        db_table = u'stockcollection'


class Stock(models.Model):
    stock_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    uniquename = models.TextField()
    description = models.TextField()
    stockcollection = models.ForeignKey(Stockcollection)
    feature = models.ForeignKey("Feature", related_name='stocks', null=True)

    class Meta:
        db_table = u'stock'


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
    prefix = models.CharField(max_length=4)
    suffix = models.CharField(max_length=4)
    genbank_file = models.FileField(upload_to=settings.GENBANK_DIR)
    timecreation = models.DateTimeField(auto_now_add=True)
    timelastmodified = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = u'feature'

    @property
    def url(self):
        urlprefix = self.dbxref.db.urlprefix
        feat_dir = 'feature'
        accession = self.dbxref.accession
        return '{0}{1}/{2}'.format(urlprefix, feat_dir, accession)

    @property
    def seq_len(self):
        return len(self.residues)

    @property
    def props(self):
        'It returns a list of cvterms'
        new_prop_dict = {}
        props = Featureprop.objects.filter(feature=self)
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
    def enzyme_in(self):
        'It returns the enzyme in of the feature'
        return self.props.get(ENZYME_IN_TYPE_NAME, None)

    @property
    def enzyme_out(self):
        'It returns the enzyme out of the feature'
        if self.type.name == VECTOR_TYPE_NAME:
            return self.props[ENZYME_OUT_TYPE_NAME]
        else:
            return self.vector.enzyme_out

    @property
    def resistance(self):
        'It returns the resistance of the feature'
        if self.type.name == VECTOR_TYPE_NAME:
            return self.props[RESISTANCE_TYPE_NAME]
        else:
            return self.vector.resistance

    @property
    def description(self):
        'Get description if it has one'
        try:
            return self.props[DESCRIPTION_TYPE_NAME][0]
        except KeyError:
            return None

    @property
    def reference(self):
        'Get references of the feature if it has them'
        try:
            return self.props[REFERENCE_TYPE_NAME][0]
        except KeyError:
            return None

    @property
    def direction(self):
        'direction of the feature'
        if self.type.name != VECTOR_TYPE_NAME:
            vec_suffix = self.vector.suffix
            vec_prefix = self.vector.prefix
        else:
            vec_suffix = self.suffix
            vec_prefix = self.prefix
        if vec_prefix == 'CGCT' and vec_suffix == 'GGAG':
            direction = FORWARD
        elif vec_prefix == 'CTCC' and vec_suffix == 'AGCG':
            direction = REVERSE
        else:
            direction = None
        return direction

    @property
    def owner(self):
        'owner of the feat'
        return FeaturePerm.objects.get(feature=self).owner
        # return self.featureperm.owner

    @property
    def is_public(self):
        'owner of the feat'
        return FeaturePerm.objects.get(feature=self).is_public


class FeaturePerm(models.Model):
    'Model to store the perms of the features'
    feature = models.OneToOneField(Feature, primary_key=True)
    owner = models.ForeignKey(User)
    is_public = models.BooleanField()

    class Meta:
        db_table = u'featureperm'


class Featureprop(models.Model):
    'Model to store the properties of the features'
    featureprop_id = models.AutoField(primary_key=True)
    feature = models.ForeignKey(Feature)
    type = models.ForeignKey(Cvterm)
    value = models.TextField()
    rank = models.IntegerField()

    class Meta:
        db_table = u'featureprop'


