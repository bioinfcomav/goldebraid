
from django.conf.urls import patterns, url
from goldenbraid.views import add_feature_view, feature_view

urlpatterns = patterns('',
            url(r'^add/feature', add_feature_view, name='add_feature'),
            url(r'^feature/(?P<uniquename>.+)/$', feature_view,
                name='feature_view')
                       )
