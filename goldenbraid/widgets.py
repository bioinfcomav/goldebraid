'''
Created on 2015 mar. 26

@author: peio
'''

from os.path import join

from django.forms.util import flatatt
from django.forms.widgets import TextInput
from django.utils.safestring import mark_safe
from django.utils.encoding import force_unicode
from django.conf import settings


class AutocompleteTextInput(TextInput):
    '''A text input that autocompletes getting a json list'''
    class Media:
        'The css and javascript files required by this widget'
        base_url = settings.STATIC_URL
        # TODO add join url
        css = {'all': (join(base_url, 'style/base/jquery.ui.all.css'),
                       join(base_url, 'style/autocomplete_ui/autocomplete.css')
                       ,)}
#         js = (os.path.join(base_url, 'js/jquery-1.4.4.min.js'),
#               os.path.join(base_url + 'js/jquery-ui-1.8.10.custom.min.js'))

    def __init__(self, attrs=None, source=None, min_length=3,
                 result_limit=100):
        '''It inits the widget.

        A source url for the json list should be given.
        '''
        super(TextInput, self).__init__(attrs)
        if source is None:
            raise ValueError('A source url should be given')
        self.source = source
        self.min_length = int(min_length)
        self.result_limit = result_limit

    def render(self, name, value, attrs=None):
        'It renders the html and the javascript'
        if value is None:
            value = ''
        final_attrs = self.build_attrs(attrs, type=self.input_type, name=name)
        if value != '':
            # Only add the 'value' attribute if a value is non-empty.
            final_attrs['value'] = force_unicode(self._format_value(value))

        html = u'<input%s />' % flatatt(final_attrs)
        javascript = self.render_js(attrs['id'])
        return mark_safe(html + javascript)

    def render_js(self, field_id):
        'The javascript that does the autocomplete'

        javascript = u'''<script type="text/javascript">
$.expr[':'].textEquals = function (a, i, m) {
  return $(a).text().match("^" + m[3] + "$");
};
$(function() {
  $("#%(field_id)s").autocomplete({
    source: "%(source)s?limit=%(limit)s",
    minLength: %(min_length)i,
    change: function(event, ui) {
      // if the value of the textbox does not match a suggestion, clear the
      // textbox and the pk
      if ($(".ui-autocomplete li:textEquals('" + $(this).val() + "')").size() == 0) {
         $(this).val('');
      }
    },
  }).on('keydown', function (e) {
    var keyCode = e.keyCode || e.which;
    // if TAB or RETURN is pressed and the text in the textbox
    // does not match a suggestion clear the value of the pk and
    // the textbox
    if((keyCode == 9 || keyCode == 13) &&
       ($(".ui-autocomplete li:textEquals('" + $(this).val() + "')").size() == 0)) {
      $(this).val('');
    }
  });
});
</script>
'''
        javascript %= {'field_id': field_id, 'source': self.source,
                       'min_length': self.min_length,
                       'limit': self.result_limit}
        return javascript
