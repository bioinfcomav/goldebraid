import re

from django import template
from django.utils.safestring import mark_safe
from django.template.defaultfilters import stringfilter
from django.utils.encoding import force_unicode
from django.utils.functional import allow_lazy

register = template.Library()


def replaceunderscore(value):
    """
    Replaces underscore with spaces
    """
    return mark_safe(re.sub('_', ' ', value))

replaceunderscore.is_safe = True
slugify = stringfilter(replaceunderscore)
register.filter('replaceunderscore', replaceunderscore)


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
letterwrap.is_safe = True
letterwrap = stringfilter(letterwrap)
register.filter('letterwrap', letterwrap)


def first_item(value, separator):
    """
    It returns the key part of a string.
    """
    return value.split(separator, 1)[0]
first_item.is_safe = True
register.filter('first_item', first_item)


def not_first_item(value, separator):
    """
    It returns the key part of a string.
    """
    return ' '.join(value.split(separator, 1)[1:])
not_first_item.is_safe = True
register.filter('not_first_item', not_first_item)
