import re

from django import template
from django.utils.safestring import mark_safe
from django.utils.html import conditional_escape
from django.template.defaultfilters import stringfilter
from django.utils.encoding import force_unicode
from django.utils.functional import allow_lazy
from django.core.serializers import serialize
from django.db.models.query import QuerySet
import json
from django.contrib.auth.models import User, AnonymousUser

register = template.Library()


def link_if_url(identifier, object_, autoescape=True):
    '''Given a text and an object it returns either an str or an html <a>

    It will return an anchor if it is capable of building an url from the
    object.
    '''
    if autoescape:
        esc = conditional_escape
    else:
        esc = lambda x: x

    url = getattr(object_, 'url', None)
    identifier = esc(identifier)
    if url:
        url = esc(url)
        result = "<a href='%s'>%s</a>" % (url, identifier)
    else:
        result = identifier
    return mark_safe(result)

register.filter('link_if_url', link_if_url)


def replaceunderscore(value):
    """
    Replaces underscore with spaces
    """
    return mark_safe(re.sub('_', ' ', value))

slugify = stringfilter(replaceunderscore)
register.filter('replaceunderscore', replaceunderscore, is_safe=True)


def wrap(text, width):
    """
    A word-wrap function that preserves existing line breaks and most spaces in
    the text. Expects that existing line breaks are posix newlines.
    """
    text = force_unicode(text)

    def _lines(text):
        for index in range(0, len(text), width):
            yield text[index: index + width] + '\n'
    return u''.join(_lines(text))
wrap = allow_lazy(wrap, unicode)


def letterwrap(value, arg):
    """
    Wraps words at specified line length.
    Argument: number of characters to wrap the text at.
    """
    return wrap(value, int(arg))
letterwrap = stringfilter(letterwrap)
register.filter('letterwrap', letterwrap, is_safe=True)


def first_item(value, separator):
    """
    It returns the key part of a string.
    """
    return value.split(separator, 1)[0]
register.filter('first_item', first_item, is_safe=True)


def not_first_item(value, separator):
    """
    It returns the key part of a string.
    """
    return ' '.join(value.split(separator, 1)[1:])
register.filter('not_first_item', not_first_item, is_safe=True)


def jsonify(item):
    if isinstance(item, QuerySet):
        return serialize('json', item)
    return json.dumps(item)

register.filter('jsonify', jsonify)


def zip_lists(a, b):
    return zip(a, b)

register.filter('zip', zip_lists, is_safe=True)


def filter_private_exps(experiments, user):
    try:
        user = User.objects.get(username=user)
    except User.DoesNotExist:
        user = AnonymousUser()
    print(user)
    public_experiments = []
    print(experiments)
    for experiment in experiments:
        print(experiment.uniquename)
        if experiment.owner == user or user.is_staff or experiment.is_public:
            public_experiments.append(experiment)
    print public_experiments
    return public_experiments

register.filter('filter_private_exps', filter_private_exps, is_safe=True)

def filter_2best_images(experiments):
    urls = []
    for exp in experiments:
        exp_urls = exp.image_urls
        if exp_urls:
            urls.append(exp_urls[0])
        if len(urls) >= 2:
            break
    return urls
            
register.filter('filter_2best_images', filter_2best_images, is_safe=True)
