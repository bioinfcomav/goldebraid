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

from Bio import SeqIO

from django.db import models
from django.contrib.auth.models import User
from django.conf import settings as site_settings
from goldenbraid import settings
from goldenbraid.tags import (DESCRIPTION_TYPE_NAME, ENZYME_IN_TYPE_NAME,
                              VECTOR_TYPE_NAME, ENZYME_OUT_TYPE_NAME,
                              RESISTANCE_TYPE_NAME, REFERENCE_TYPE_NAME,
                              FORWARD, REVERSE, DERIVES_FROM, OTHER_TYPE_NAME,
                              TARGET_DICOT, TARGET_MONOCOT)
from goldenbraid.excel import plot_from_excel
from django.core.files.temp import NamedTemporaryFile
from django.core.urlresolvers import reverse
from goldenbraid.settings import (DOMESTICATION_VECTORS_IN_GB, CATEGORIES,
                                  SBOL_IMAGES, CRYSPER_CATEGORIES,
                                  MOCLO_INCOMPATIBLE_RESISTANCES)
from goldenbraid.utils import has_rec_sites, get_prefix_and_suffix_index

LEVEL_0 = '0'
LEVEL_1ALPHA = '1-alpha'
LEVEL_1_OMEGA = '1-omega'


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

    def __unicode__(self):
        return u'{}'.format(self.name)


class Dbxref(models.Model):
    dbxref_id = models.AutoField(primary_key=True)
    db = models.ForeignKey(Db)
    accession = models.CharField(max_length=255)

    class Meta:
        db_table = u'dbxref'

    @property
    def url(self):
        return self.db.urlprefix + self.accession


class Cvterm(models.Model):
    cvterm_id = models.AutoField(primary_key=True)
    cv = models.ForeignKey(Cv)
    name = models.CharField(max_length=1024)
    definition = models.TextField()

    class Meta:
        db_table = u'cvterm'

    def __unicode__(self):
        return u'{}'.format(self.name)


class Count(models.Model):
    "This clase is used to count how many entrys of type name has been created"
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

    def __unicode__(self):
        return u'{}'.format(self.uniquename)

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
        elif self.vector:
            return self.vector.enzyme_out

    @property
    def resistance(self):
        'It returns the resistance of the feature'
        if self.type.name == VECTOR_TYPE_NAME:
            return self.props[RESISTANCE_TYPE_NAME]
        else:
            if self.type.name in (TARGET_DICOT, TARGET_MONOCOT):
                return None
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
    def level(self):
        if self.type.name in (TARGET_DICOT, TARGET_MONOCOT):
            return LEVEL_0
        if not self.vector:
            return None
        vector_name = self.vector.name
        if vector_name in DOMESTICATION_VECTORS_IN_GB:
            return LEVEL_0
        elif 'alpha' in vector_name:
            return LEVEL_1ALPHA
        elif 'omega' in vector_name:
            return LEVEL_1_OMEGA

    def _get_sbol_image(self, direction=FORWARD):
        if self.level == LEVEL_0:
            if direction == FORWARD:
                return SBOL_IMAGES[FORWARD].get(self.gb_category, None)
            elif direction == REVERSE:
                return SBOL_IMAGES[REVERSE].get(self.gb_category, None)
            else:
                raise RuntimeError('Sequence direction must be forw or rev')

    @property
    def sbol_images(self):
        images = []
        level = self.level
        if level == LEVEL_0:
            images.append(self._get_sbol_image())
        elif level in (LEVEL_1_OMEGA, LEVEL_1ALPHA):
            sub_children = []
            for child in self.children:
                if child.level == LEVEL_0:
                    sub_children.append(child)
                elif level in (LEVEL_1_OMEGA, LEVEL_1ALPHA):
                    images.append(child.sbol_images)
            if sub_children:
                if self.direction == REVERSE:
                    sub_children = sub_children[::-1]
                images.extend([child._get_sbol_image(self.direction)
                               for child in sub_children])

        images = list(self._flatten_list(images))
        # Sometimes we have old pieces and we can not know the type
        if any([True if x is None else False for x in images]):
            return None
        return images

    def _flatten_list(self, list_):
        for item in list_:
            if hasattr(item, '__iter__'):
                for item2 in self._flatten_list(item):
                    yield item2
            else:
                yield item

    @property
    def gb_version(self):
        category = (self.type.name, self.prefix, self.suffix)
        v2_categories = [(u'CDS', u'AATG', u'GCAG'), (u'CT', u'GCAG', u'GCTT')]
        if category in v2_categories:
            return '2.0'
        else:
            return '3.0'

    @property
    def moclo_compatible(self):
        if not self.level or self.level != LEVEL_0 or self.type.name in (TARGET_DICOT, TARGET_MONOCOT):
            return 'not_evaluable'
        if self.vector.resistance not in (MOCLO_INCOMPATIBLE_RESISTANCES):
            enzyme = self.enzyme_out[0]
            residues = self.residues
            pref_idx, suf_idx = get_prefix_and_suffix_index(residues, enzyme)[:2]
            seq = residues[pref_idx:suf_idx + len(self.suffix)]

            # TODO: maybe we should look only to the part seq. not with the vector
            return not has_rec_sites(seq, enzymes=('BpiI', 'BsaI'))
        else:
            return False

    @property
    def gb_category(self):
        if (self.level != LEVEL_0 or  self.type.name == VECTOR_TYPE_NAME or
                self.type.name == OTHER_TYPE_NAME):
            return self.type.name
        type_ = self.type.name
        prefix = self.prefix
        suffix = self.suffix

        for category, values in CATEGORIES.items():
            if values == (type_, prefix, suffix):
                return category.strip()
        for category, values in CRYSPER_CATEGORIES.items():
            if values == (type_, prefix, suffix):
                return category.strip()

    @property
    def gb_category_sections(self):
        gb_category = self.gb_category
        if gb_category == OTHER_TYPE_NAME:
            return None

        if gb_category is not None:
            return gb_category.split(' ')[-1].strip()

    @property
    def gb_category_name(self):
        gb_category = self.gb_category
        if gb_category == OTHER_TYPE_NAME:
            return gb_category
        if gb_category is not None:
            gb_category_split = ' '.join(gb_category.split(' ')[:-1]).strip()
            if gb_category_split:
                return gb_category_split
            else:
                return gb_category

    @property
    def owner(self):
        'owner of the feat'
        return FeaturePerm.objects.get(feature=self).owner
        # return self.featureperm.owner

    @property
    def is_public(self):
        'owner of the feat'
        return FeaturePerm.objects.get(feature=self).is_public

    @property
    def children(self):
        children = []
        for frls in FeatureRelationship.objects.filter(object=self):
            children.append(frls.subject)
        # sort
        current_suffix = None
        ordered_children = []
        while children:
            for index, child in enumerate(children):
                if ((current_suffix is None and child.prefix == 'GGAG') or
                        (current_suffix and current_suffix == child.prefix)):
                    current_suffix = child.suffix
                    break

            ordered_children.append(children.pop(index))
        return ordered_children

    def add_relations(self, seq=None):
        """If seq object is not given, it will look in the genbank file
        associated with the feature"""
        if seq is None:
            path = self.genbank_file.path
            if os.path.exists(path):
                seq = SeqIO.read(path, 'gb')
        children = _parse_children_relations_from_gb(seq)
        if not children:
            return

        for child in children:
            try:
                child = Feature.objects.get(uniquename=child)
            except Feature.DoesNotExist:
                continue
            _get_or_create_feature_relationship(object_=self, subject=child)

    @property
    def ordered_experiments(self):
        experiments = []
        for featsubexpe in ExperimentSubFeature.objects.filter(feature=self):
            exp = featsubexpe.experiment
            if exp not in experiments:
                experiments.append(exp)
        for featexpe in ExperimentFeature.objects.filter(feature=self):
            exp = featexpe.experiment
            if exp not in experiments:
                experiments.append(exp)
        return experiments

    @property
    def experiment_images(self, user):
        experiments = self.ordered_experiments
        image_urls = []
        for experiment in experiments:
            urls = experiment.image_urls
            if urls:
                image_urls.append(urls[0])
            if len(image_urls) >= 2:
                break
        return image_urls


def _parse_children_relations_from_gb(seq):
    definition = seq.description
    if '(' in definition and ')' in definition:
        match = re.search('\((.+)\)', definition)
        if match:
            return match.group(1).split(',')
    else:
        return None


def _get_or_create_feature_relationship(object_, subject):
    derives_from = Cvterm.objects.get(name=DERIVES_FROM)
    try:
        FeatureRelationship.objects.get(subject=subject,
                                        type=derives_from, object=object_)
    except FeatureRelationship.DoesNotExist:
        FeatureRelationship.objects.create(subject=subject,
                                           type=derives_from, object=object_)


class FeaturePerm(models.Model):
    'Model to store the perms of the features'
    feature = models.OneToOneField(Feature, primary_key=True)
    owner = models.ForeignKey(User)
    is_public = models.BooleanField(default=False)

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


class FeatureRelationship(models.Model):
    'Store the relationships between parts'
    featurerelationship_id = models.AutoField(primary_key=True)
    type = models.ForeignKey(Cvterm)
    subject = models.ForeignKey(Feature, related_name='subject')
    object = models.ForeignKey(Feature, related_name='object')

    class Meta:
        db_table = u'feature_relationship'


class Experiment(models.Model):
    'Store the experiments associated to features'
    experiment_id = models.AutoField(primary_key=True)
    uniquename = models.CharField(max_length=255, unique=True)
    chasis_1 = models.CharField(max_length=255)
    chasis_2 = models.CharField(max_length=255)
    description = models.TextField(max_length=3000)
    type = models.ForeignKey(Cvterm)
    timecreation = models.DateTimeField(auto_now_add=True)
    dbxref = models.ForeignKey(Dbxref)

    class Meta:
        db_table = 'experiment'

    @property
    def features_used_in_experiment(self):
        try:
            exp_feats = ExperimentFeature.objects.filter(experiment=self)
        except ExperimentFeature.DoesNotExist:
            exp_feats = None
        if exp_feats:
            return [exp_feat.feature for exp_feat in exp_feats]
        else:
            return []

    @property
    def key_features(self):
        try:
            exp_subfeats = ExperimentSubFeature.objects.filter(experiment=self)
        except ExperimentSubFeature.DoesNotExist:
            exp_subfeats = None
        if exp_subfeats:
            return [exp_subfeat.feature for exp_subfeat in exp_subfeats]

    @property
    def owner(self):
        'owner of the feat'
        return ExperimentPerm.objects.get(experiment=self).owner
        # return self.featureperm.owner

    @property
    def is_public(self):
        'owner of the feat'
        return ExperimentPerm.objects.get(experiment=self).is_public

    @is_public.setter
    def is_public(self, value):
        exp_perm = ExperimentPerm.objects.get(experiment=self)
        exp_perm.is_public = value
        exp_perm.save()

    @property
    def url(self):
        urlprefix = self.dbxref.db.urlprefix
        feat_dir = 'experiment'
        accession = self.dbxref.accession
        return '{0}{1}/{2}'.format(urlprefix, feat_dir, accession)

    @property
    def numeric_props(self):
        prop_dict = {}
        props = ExperimentPropNumeric.objects.filter(experiment=self).order_by('type__name')
        for prop in props:
            type_ = prop.type.name
            value = prop.value
            if type_ not in prop_dict:
                prop_dict[type_] = []
            prop_dict[type_].append(value)
        return sorted(prop_dict.iteritems())
        # return ReadOnlyDict(prop_dict)

    @property
    def text_props(self):
        prop_dict = {}
        props = ExperimentPropText.objects.filter(experiment=self)
        for prop in props:
            title = prop.title
            value = prop.value
            if title not in prop_dict:
                prop_dict[title] = []
            prop_dict[title].append(value)

        return ReadOnlyDict(prop_dict)

    @property
    def image_props(self):
        return [(image_prop.description, image_prop.image)
                for image_prop in ExperimentPropImage.objects.filter(experiment=self)]

    @property
    def excel_props(self):
        return ExperimentPropExcel.objects.filter(experiment=self)

    @property
    def image_urls(self):
        urls = []
        for exp_excel in self.excel_props:
            url = exp_excel.image_url
            alt = exp_excel.description
            urls.append((url, alt))
        for image_desc, image in self.image_props:
            url = image.url
            alt = image_desc
            urls.append((url, alt))
        return urls

    @property
    def generic_file_props(self):
        return [(generic_file_prop.description, generic_file_prop.file)
                for generic_file_prop in ExperimentPropGenericFile.objects.filter(experiment=self)]

    @property
    def keywords(self):
        exp_keywords = ExperimentKeyword.objects.filter(experiment=self)
        keywords = [exp_keyword.keyword for exp_keyword in exp_keywords]
        return keywords


class ExperimentPerm(models.Model):
    experiment = models.OneToOneField(Experiment, primary_key=True)
    owner = models.ForeignKey(User)
    is_public = models.BooleanField(default=False)

    class Meta:
        db_table = u'experimentperm'


class ExperimentFeature(models.Model):
    experiment_feature_id = models.AutoField(primary_key=True)
    experiment = models.ForeignKey(Experiment)
    feature = models.ForeignKey(Feature)

    class Meta:
        db_table = u'experimentfeature'
        unique_together = ('experiment', 'feature')


class ExperimentSubFeature(models.Model):
    experiment_key_subfeature_id = models.AutoField(primary_key=True)
    experiment = models.ForeignKey(Experiment)
    feature = models.ForeignKey(Feature)

    class Meta:
        db_table = u'experimentkeysubfeature'
        unique_together = ('experiment', 'feature')


class ExperimentPropNumeric(models.Model):
    experiment_prop_numeric_id = models.AutoField(primary_key=True)
    experiment = models.ForeignKey(Experiment)
    type = models.ForeignKey(Cvterm)
    value = models.FloatField(null=True)

    class Meta:
        db_table = u'experimentpropnumeric'


class ExperimentPropText(models.Model):
    experiment_prop_text_id = models.AutoField(primary_key=True)
    experiment = models.ForeignKey(Experiment)
    title = models.CharField(max_length=255)
    value = models.TextField(max_length=3000)

    class Meta:
        db_table = u'experimentproptext'


# TODO: File name must be unique
class ExperimentPropImage(models.Model):
    experiment_prop_image_id = models.AutoField(primary_key=True)
    experiment = models.ForeignKey(Experiment)
    image = models.ImageField(upload_to=settings.RESULTS_DIR)
    description = models.CharField(max_length=255)

    class Meta:
        db_table = u'experimentpropimage'


class ExperimentPropExcel(models.Model):
    experiment_prop_excel_id = models.AutoField(primary_key=True)
    experiment = models.ForeignKey(Experiment)
    excel = models.FileField(upload_to=settings.RESULTS_DIR)
    description = models.CharField(max_length=255)

    class Meta:
        db_table = u'experimentpropexcel'

    @property
    def drawed_image(self):
        temp_fhand = NamedTemporaryFile()
        plot_from_excel(self.excel.path, temp_fhand)
        content_type = 'image/svg+xml'
        return open(temp_fhand.name).read(), content_type

    @property
    def image_url(self):
        return reverse('api_excel_image', args=[self.experiment_prop_excel_id])


class ExperimentPropGenericFile(models.Model):
    experiment_prop_generic_file_id = models.AutoField(primary_key=True)
    experiment = models.ForeignKey(Experiment)
    file = models.FileField(upload_to=settings.RESULTS_DIR)
    description = models.CharField(max_length=255)

    class Meta:
        db_table = u'experimentpropgenericfile'


class ExperimentKeyword(models.Model):
    experiment_keyword_id = models.AutoField(primary_key=True)
    experiment = models.ForeignKey(Experiment)
    keyword = models.CharField(max_length=50)

    class Meta:
        db_table = u'experimentkeyword'
