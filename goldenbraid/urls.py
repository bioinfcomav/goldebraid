# Copyright 2013 Diego Orzaez, Univ.Politecnica Valencia, Consejo Superior de
# Investigaciones Cientificas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from django.conf.urls import patterns, url
from goldenbraid.views.feature import (add_feature_view, feature_view,
                                       add_vector_view, search_features_view)

from goldenbraid.views.api import(feature_uniquenames, features_children,
                                  features_key_elements, excel_image,
                                  experiment_keywords, feature_sbol,
                                  combined_excel_image)

from goldenbraid.views.multipartite import (multipartite_view,
                                            multipartite_protocol_view,
                                            multipartite_view_genbank,
                                            multipartite_view_sbol,
                                            multipartite_view_free,
                                            multipartite_view_free_genbank,
                                            multipartite_view_free_sbol,
                                            multipartite_view_free_protocol,
                                            multipartite_view_add)
from goldenbraid.views.bipartite import (bipartite_view,
                                         bipartite_view_genbank,
                                         bipartite_view_sbol,
                                         bipartite_view_protocol,
                                         bipartite_view_add)
from goldenbraid.views.domestication import (domestication_view,
                                             domestication_view_genbank,
                                             domestication_view_protocol,
                                             domestication_view_add,
                                             domestication_view_sbol,
                                             synthesis_view,
                                             synthesis_view_genbank,
                                             synthesis_view_sbol,
                                             synthesis_view_protocol,
                                             crispr_view,
                                             crispr_view_protocol,
                                             crispr_view_add)
from goldenbraid.views.experiment import (search_experiment,
                                          add_experiment_view, experiment_view)
from goldenbraid.views.user import user_view


urlpatterns = patterns('',
    url(r'^user/(?P<username>.+)/$', user_view, name='user_view'),  # @IgnorePep8

    url(r'^do/domestication/$', domestication_view, name='domestication_view'),
    url(r'^do/domestication/add/$', domestication_view_add, name='domestication_view_add'),  # @IgnorePep8
    url(r'^do/domestication/genbank/$', domestication_view_genbank, name='domestication_view_genbank'),  # @IgnorePep8
    url(r'^do/domestication/sbol/$', domestication_view_sbol, name='domestication_view_sbol'),  # @IgnorePep8
    url(r'^do/domestication/protocol/$', domestication_view_protocol, name='domestication_view_protocol'),  # @IgnorePep8

    url(r'^do/synthesis/$', synthesis_view, name='synthesis_view'),
    url(r'^do/synthesis/genbank/$', synthesis_view_genbank, name='shinthesis_view_genbank'),  # @IgnorePep8
    url(r'^do/synthesis/sbol/$', synthesis_view_sbol, name='shinthesis_view_sbol'),  # @IgnorePep8
    url(r'^do/synthesis/protocol/$', synthesis_view_protocol, name='synthesis_view_protocol'),  # @IgnorePep8

    url(r'^do/crispr/$', crispr_view, name='crispr_view'),
    url(r'^do/crispr/protocol/$', crispr_view_protocol, name='crispr_view_add'),  # @IgnorePep8
    url(r'^do/crispr/add/$', crispr_view_add, name='crispr_view_add'),

    url(r'^add/vector/$', add_vector_view, name='add_vector'),
    url(r'^add/feature/$', add_feature_view, name='add_feature'),
    url(r'^add/experiment/(?P<exp_type>.+)/?$', add_experiment_view, name='add_experiment'),  # @IgnorePep8

    url(r'^search/experiment/$', search_experiment, name='search_experiments'),
    url(r'^search/features/$', search_features_view, name='search_features'),

    url(r'^experiment/(?P<uniquename>.+)/$', experiment_view, name='experiment_view'),  # @IgnorePep8
    url(r'^feature/(?P<uniquename>.+)/$', feature_view, name='feature_view'),

    url(r'^do/multipartite/add/?$', multipartite_view_add, name="multipartite_view_add"),  # @IgnorePep8
    url(r'^do/multipartite/free/protocol/?$', multipartite_view_free_protocol, name="multipartite_view_free_protocol"),  # @IgnorePep8
    url(r'^do/multipartite/free/genbank/?$', multipartite_view_free_genbank, name="multipartite_view_free_genbank"),  # @IgnorePep8
    url(r'^do/multipartite/free/sbol/?$', multipartite_view_free_sbol, name="multipartite_view_free_sbol"),  # @IgnorePep8
    url(r'^do/multipartite/free/(?P<form_num>.+)?/?$', multipartite_view_free, name="multipartite_view_free"),  # @IgnorePep8
    url(r'^do/multipartite/protocol/$', multipartite_protocol_view, name="multipartite_view_protocol"),  # @IgnorePep8
    url(r'^do/multipartite/(?P<multi_type>.+)?/genbank/$', multipartite_view_genbank, name="multipartite_view_genbank"),  # @IgnorePep8
    url(r'^do/multipartite/(?P<multi_type>.+)?/sbol/$', multipartite_view_sbol, name="multipartite_view_sbol"),  # @IgnorePep8
    url(r'^do/multipartite/(?P<multi_type>.+)?/?$', multipartite_view, name='multipartite_view'),  # @IgnorePep8

    url(r'^do/bipartite/add/?$', bipartite_view_add, name='bipartite_view_add'),  # @IgnorePep8
    url(r'^do/bipartite/genbank/?$', bipartite_view_genbank, name='bipartite_view_genbank'),  # @IgnorePep8
    url(r'^do/bipartite/sbol/?$', bipartite_view_sbol, name='bipartite_view_sbol'),  # @IgnorePep8
    url(r'^do/bipartite/protocol/?$', bipartite_view_protocol, name='bipartite_view_protocol'),  # @IgnorePep8
    url(r'^do/bipartite/(?P<form_num>.+)?/?$', bipartite_view, name='bipartite_view'),  # @IgnorePep8

    url('api/feature_uniquenames/$', feature_uniquenames, name='api_feature'),
    url('api/sbol/(?P<uniquename>.+)/$', feature_sbol, name='api_feature_sbol'),

    url('api/features_children/$', features_children, name='api_feature_children'),  # @IgnorePep8
    url('api/features_key_elements/$', features_key_elements, name='api_feature_key_elements'),  # @IgnorePep8
    url('api/excel_graph/(?P<excel_id>.+)?', excel_image, name='api_excel_image'),  # @IgnorePep8
    url('api/exp_keywords/$', experiment_keywords, name='api_exp_keywords'),  # @IgnorePep8
    url('api/excel_combined/(?P<uniquename>.+)?/(?P<exp_type>.+)?/.svg', combined_excel_image, name='api_combined_excel_images'),  # @IgnorePep8
                       )
