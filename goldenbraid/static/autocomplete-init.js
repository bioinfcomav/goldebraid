$.expr[':'].textEquals = function (a, i, m) {
  return $(a).text().match("^" + m[3] + "$");
};
function enableAutocomplete(context, source) {
    $('input', context || null).autocomplete({
            source: source,
            minLength: 1,
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
};
function enableAutocomplete_nocheck(context, source) {
    $('input', context || null).autocomplete({
            source: source,
            minLength: 1,
    })
}