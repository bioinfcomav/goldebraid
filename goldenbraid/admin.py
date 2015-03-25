'''
Created on 2015 mar. 13

@author: peio
'''
from django.contrib import admin
from goldenbraid.models import Cvterm, Cv


class CvAdmin(admin.ModelAdmin):
    list_display = ('name', 'definition')
admin.site.register(Cv, CvAdmin)


class CvtermAdmin(admin.ModelAdmin):
    list_display = ('name', 'definition')
admin.site.register(Cvterm, CvtermAdmin)

