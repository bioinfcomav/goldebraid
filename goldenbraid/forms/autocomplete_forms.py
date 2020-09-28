from dal import autocomplete

from goldenbraid.models import Feature
from django import forms 
from goldenbraid.forms.feature import create_feature_validator 


class Bipart1Form(forms.ModelForm):
    class Meta:
        model = Feature
        fields = ('uniquename',)
        labels = {'uniquename': "Select Part 1"}
        widgets = {'uniquename': autocomplete.ListSelect2('bipartite1_autocomplete', 
                                                          attrs={'style':'width:1000px',
                                                                 }),
                  }

    def clean(self):
        pass

    def clean_uniquename(self):
        return create_feature_validator('uniquename')(self)


class Bipart2Form(forms.ModelForm):
    part1 = forms.CharField(max_length=100)
    part1.widget.attrs['readonly'] = True

    class Meta:
        model = Feature
        fields = ('uniquename',)
        labels = {'uniquename': "Select Part 2"}
        widgets = {'uniquename': autocomplete.ListSelect2('bipartite2_autocomplete', 
                                                          attrs={'style':'width:1000px',
                                                                 },
                                                          forward=['part1'])
                  }
    def clean(self):
        pass

    def clean_uniquename(self):
        return create_feature_validator('uniquename')(self)


class Bipart3Form(forms.ModelForm):
    part1 = forms.CharField(max_length=100)
    part1.widget.attrs['readonly'] = True
    part2 = forms.CharField(max_length=100)
    part2.widget.attrs['readonly'] = True
    class Meta:
        model = Feature
        fields = ('uniquename',)
        labels = {'uniquename': "Select Vector"}
        widgets = {'uniquename': autocomplete.ListSelect2('bipartite3_autocomplete', 
                                                          attrs={'style':'width:1000px',
                                                                 },
                                                          forward=['part1'])
                  }
    def clean(self):
        pass

    def clean_uniquename(self):
        return create_feature_validator('uniquename')(self)
