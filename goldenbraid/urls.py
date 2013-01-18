
from django.conf.urls import patterns, url
from goldenbraid.views.feature_views import add_feature_view, feature_view
from goldenbraid.views.multipartite_views import multipartite_view

urlpatterns = patterns('',
        url(r'^add/feature', add_feature_view, name='add_feature'),
        url(r'^feature/(?P<uniquename>.+)/$', feature_view,
                name='feature_view'),
        url(r'do/multipartite/(?P<multi_type>.+)?/?$', multipartite_view,
            name='multipartite_view')
                       )
