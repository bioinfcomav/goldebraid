'''
Created on 2015 mar. 26

@author: peio
'''

from os.path import join

from django.forms.utils import flatatt
from django.forms.widgets import TextInput, SelectMultiple
from django.utils.safestring import mark_safe
from django.utils.encoding import force_unicode
from django.conf import settings
from django.utils.html import format_html


class AutocompleteTextInput(TextInput):
    '''A text input that autocompletes getting a json list'''
    class Media:
        'The css and javascript files required by this widget'
        base_url = settings.STATIC_URL
        # TODO add join url
        css = {'all': (join(base_url, 'style/base/jquery.ui.all.css'),
                       join(base_url, 'style/autocomplete_ui/autocomplete.css')
                       ,)}

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
  })

});
</script>
'''
        javascript %= {'field_id': field_id, 'source': self.source,
                       'min_length': self.min_length,
                       'limit': self.result_limit}
        return javascript


class DinamicSelectMultiple(SelectMultiple):
    '''A multiple choice widgte that populates dinamically'''
    def __init__(self, attrs=None, source=None, parent_class=None, choices=()):
        '''It inits the widget.

        A source url for the json list should be given.
        '''
        super(SelectMultiple, self).__init__(attrs,)
        if source is None:
            raise ValueError('A source url should be given')
        self.source = source
        self._parent_class = parent_class
        self.choices = list(choices)

    def render(self, name, value, attrs=None, choices=()):
        if value is None:
            value = []
        final_attrs = self.build_attrs(attrs, name=name)
        output = [format_html('<select multiple="multiple"{0}>',
                              flatatt(final_attrs))]
        options = self.render_options(choices, value)
        if options:
            output.append(options)
        output.append('</select>')
        javascript = self.render_js(attrs['id'])
        return mark_safe('\n'.join(output) + javascript)

    def render_js_old(self, field_id):
        'The javascript that does the select options'

        javascript = u'''<script type="text/javascript">
$(document).ready(function() {
inputs = ssacar ldlde fdom inputys
for input in finputs:
    fillalnfsdlk(input)

    $("input").focusout(function(){
        if ($(this).hasClass('%(parent_class)s')) {
            //buscar todos los input con esa classe y mirar que no vacio
            var values = [];
            $('.ui-autocomplete-input').each(function(){
                var value = $(this).val();
                if(value != ''){ values.push(value);}
            })
            $.getJSON("%(source)s", {'features':values}, function(result){
                 var toAppend = '';
                $.each(result, function(i, val){
                    toAppend += '<option>'+val+'</option>';
                });
                $('#%(field_id)s').empty().append(toAppend);
            })
        }

    });
});
</script>'''
        javascript %= {'field_id': field_id, 'source': self.source,
                       'parent_class': self._parent_class}
        return javascript

    def render_js(self, field_id):
        'The javascript that does the select options'

        javascript = u'''<script type="text/javascript">
(function( $ ){
    $.fn.fill_children = function (index, input_element){
        if ($(input_element).hasClass('%(parent_class)s')) {
            //buscar todos los input con esa classe y mirar que no vacio
            var values = [];
            $('.ui-autocomplete-input').each(function(){
                var value = $(this).val();
                if(value != ''){ values.push(value);}
            })
            $.getJSON("%(source)s", {'features':values}, function(result){
                 var toAppend = '';
                $.each(result, function(i, val){
                    toAppend += '<option>'+val+'</option>';
                });
                $('#%(field_id)s').empty().append(toAppend);
            })
        }
    }
})( jQuery );

$(document).ready(function() {
    var inputs = $(this).find("input");
    $.each(inputs, function(index, val){
        $.fn.fill_children(index, val)});
    $("input").focusout(function(){
                    $.fn.fill_children('0', $(this))});
});
</script>'''
        javascript %= {'field_id': field_id, 'source': self.source,
                       'parent_class': self._parent_class}
        return javascript


#  .on('keydown', function (e) {
#     var keyCode = e.keyCode || e.which;
#     // if TAB or RETURN is pressed and the text in the textbox
#     // does not match a suggestion clear the value of the pk and
#     // the textbox
#     if((keyCode == 9 || keyCode == 13) &&
#        ($(".ui-autocomplete li:textEquals('" + $(this).val() + "')").size() == 0)) {
#       $(this).val('');
#     }
#   });
