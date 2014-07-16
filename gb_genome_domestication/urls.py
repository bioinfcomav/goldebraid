'''
Created on 2014 uzt 1

@author: peio
'''
from django.conf.urls import patterns, url
from gb_genome_domestication.views import feature_view, search_view

urlpatterns = patterns('',
                       url(r'search/$', search_view, name='search_view'),
                       url(r'^feature/(?P<uniquename>.+)/$', feature_view,
                           name='feature_view'),
                       )
