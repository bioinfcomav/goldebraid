
from django.conf.urls import patterns, url
from goldenbraid.views.feature_views import add_feature_view, feature_view
from goldenbraid.views.feature_search_view import search_features_view
from goldenbraid.views.multipartite_views import (multipartite_view,
                                                  multipartite_protocol_view,
                                                  multipartite_view_genbank)
from goldenbraid.views.bipartite_views import (bipartite_view,
                                               bipartite_view_genbank,
                                               bipartite_view_protocol)
from goldenbraid.views.domestication_view import (domestication_view,
                                                  domestication_view_genbank,
                                                  domestication_view_protocol)


urlpatterns = patterns('',
        url(r'^add/feature', add_feature_view, name='add_feature'),
        url(r'^search/features/$', search_features_view,
            name='search_features'),
        url(r'^feature/(?P<uniquename>.+)/$', feature_view,
                name='feature_view'),
        url(r'do/multipartite/protocol/$', multipartite_protocol_view,
            name="multipartite_view_protocol"),
        url(r'do/multipartite/(?P<multi_type>.+)?/genbank/$',
            multipartite_view_genbank,
            name="multipartite_view_genbank"),
         url(r'do/multipartite/(?P<multi_type>.+)?/?$', multipartite_view,
            name='multipartite_view'),
        url(r'do/bipartite/genbank/$', bipartite_view_genbank,
            name='bipartite_view_genbank'),
        url(r'do/bipartite/protocol/$', bipartite_view_protocol,
            name='bipartite_view_protocol'),
        url(r'do/bipartite/(?P<form_num>.+)?/?$', bipartite_view,
            name='bipartite_view'),
        url(r'do/domestication/genbank/$', domestication_view_genbank,
            name='domestication_view_genbank'),
        url(r'do/domestication/protocol/$', domestication_view_protocol,
            name='domestication_view_protocol'),
        url(r'do/domestication/$', domestication_view,
            name='domestication_view'),


                       )

