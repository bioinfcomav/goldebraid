#code taken from:
# http://djangosnippets.org/snippets/2332/

from django import template

register = template.Library()

class QuerystringNode(template.Node):
    def __init__(self, base_url, query_params, override=None, var_value=None):

        self.base_url = template.Variable(base_url)
        self.query_params = template.Variable(query_params)
        if var_value is not None: # a variable was passed in instead of a literal
            self.var_value = template.Variable(var_value)
        else:
            self.var_value = None
        self.override = override

    def render(self, context):
        base_url = self.base_url.resolve(context)
        url = base_url + '?'

        params = self.query_params.resolve(context).copy()

        if self.override is not None:
            if self.var_value is None:
                k,v = self.override.split('=')
            else:
                var_value = self.var_value.resolve(context)
                k = self.override
                v = var_value
            params[k] = v

        # TODO url encode. meh
        for k,v in params.items():
            url += "%s=%s&" % (k,v)
        return url[:-1]


@register.tag
def build_url(parser, token):
    """
    Entry point for the build_url template tag. This tag allows you to maintain
    a set of default querystring parameters and override an individual param.

    It was written to support list views that need to keep sort, filter, and page
    parameters, and this is the most common use case.

    Usage:

        {% build_url base_url query_params override_key=override_value_literal %}

                  base_url: string variable -- the URL's prefix
                            try {% url some-url as base_url %}
              query_params: dictionary of default querystring values.
                            {'k1':'v1', 'k2':'mountain'}
                            -> ?k1=v1&k2=mountain
              override_key: key to replace in query_params when building the url.
    override_value_literal: literal (string pls) value for the override
                  (output): (string) the url

    There is a second form which allows you to pass in a variable
    for the override value:

        {% build_url base_url query_params override_key override_value_variable %}
                                                       ^ that's the difference
    """
    try:
        args = token.split_contents()

        base_url = args[1]
        query_params = args[2]

        override = args[3] if len(args) >= 4 else None

        var_value = args[4] if len(args) == 5 else None

    except ValueError:
        raise template.TemplateSyntaxError("%r tag requires 3 or 4 arguments" % token.contents.split()[0])

    return QuerystringNode(base_url, query_params, override, var_value)
