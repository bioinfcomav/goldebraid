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

from django.conf.urls import url, include
from goldenbraid.views.feature import (add_feature_view, feature_view,
                                       add_vector_view, search_features_view)

from goldenbraid.views.api import (feature_uniquenames, features_children,
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
                                            multipartite_view_add,
                                            level_0_genbank,
                                            level_0_sbol,
                                            level_0_add,
                                            level_0_protocol,
                                            cas12a_multiplexing_view_genbank,
                                            cas12a_multiplexing_view_sbol,
                                            cas12a_multiplexing_view_add,
                                            crispr_view_cas12_single_TU,
                                            crispr_view_cas12a_multiplexing_TU,
                                            multipartite_fungal_tu_view,
                                            fungal_gene_disruption_view)
from goldenbraid.views.bipartite import (bipartite_view,
                                         bipartite_view_genbank,
                                         bipartite_view_sbol,
                                         bipartite_view_protocol,
                                         bipartite_view_add,
                                         fungal_bipartite_view)
from goldenbraid.views.domestication import (domestication_view,
                                             domestication_view_genbank,
                                             domestication_view_protocol,
                                             domestication_view_add,
                                             domestication_view_sbol,
                                             synthesis_view,
                                             synthesis_view_genbank,
                                             synthesis_view_sbol,
                                             synthesis_view_protocol,
                                             crispr_view_cas9_single,
                                             crispr_view_cas9_multiplexing,
                                             crispr_view_cas12_single,
                                             crispr_view_protocol,
                                             crispr_view_add,
                                             crispr_view_cas12_multiplexing_2X,
                                             crispr_view_cas12_multiplexing_3X,
                                             crispr_view_cas12_multiplexing_4X,
                                             crispr_view_cas12_multiplexing_5X,
                                             crispr_view_cas12_multiplexing_6X,
                                             cas12a_multiplexing_view_protocol,
                                             fungal_domestication_view,
                                             fungal_domestication_view_add,
                                             fungal_domestication_view_genbank,
                                             fungal_domestication_view_protocol,
                                             fungal_domestication_view_sbol)

from goldenbraid.views.crispr_multiplexing import (level_0_editing,
                                                   level_0_regulation,
                                                   crispr_multiplexing,
                                                   crispr_view_cas9_multiplexing_auto,
                                                   crispr_view_cas9_multiplexing_auto_add_tu,
                                                   crispr_view_cas9_multiplexing_auto_protocol,
                                                   cripsr_view_cas9_auto_omega_add,
                                                   cripsr_view_cas9_auto_omega_protocol)
from goldenbraid.views.autocomplete_views import BipartitePart1Autocomplete, BipartitePart2Autocomplete, BipartitePart3Autocomplete
from goldenbraid.views.experiment import (search_experiment,
                                          add_experiment_view, experiment_view)
from goldenbraid.views.user import user_view
from goldenbraid.views.viewsets import search_feature_index, FeatureViewSet
from rest_framework import routers

router = routers.DefaultRouter()
router.register(r'search_features', FeatureViewSet, basename='features')

urlpatterns = [
    url(r'^search/features/api/', include(router.urls)),
    url(r'^user/(?P<username>.+)/$', user_view, name='user_view'),  # @IgnorePep8
    url(r'^do/domestication/$', domestication_view, name='domestication_view'),
    url(r'^do/domestication/add/$', domestication_view_add, name='domestication_view_add'),  # @IgnorePep8
    url(r'^do/domestication/genbank/$', domestication_view_genbank, name='domestication_view_genbank'),  # @IgnorePep8
    url(r'^do/domestication/sbol/$', domestication_view_sbol, name='domestication_view_sbol'),  # @IgnorePep8
    url(r'^do/domestication/protocol/$', domestication_view_protocol, name='domestication_view_protocol'),  # @IgnorePep8
    url(r'^fungal/do/domestication/$', fungal_domestication_view, name='domestication_view'),
    url(r'^fungal/do/domestication/add/$', fungal_domestication_view_add, name='domestication_view_add'),  # @IgnorePep8
    url(r'^fungal/do/domestication/genbank/$', fungal_domestication_view_genbank, name='domestication_view_genbank'),  # @IgnorePep8
    url(r'^fungal/do/domestication/sbol/$', fungal_domestication_view_sbol, name='domestication_view_sbol'),  # @IgnorePep8
    url(r'^fungal/do/domestication/protocol/$', fungal_domestication_view_protocol, name='domestication_view_protocol'),  # @IgnorePep8
    url(r'^do/synthesis/$', synthesis_view, name='synthesis_view'),
    url(r'^do/synthesis/genbank/$', synthesis_view_genbank, name='shinthesis_view_genbank'),  # @IgnorePep8
    url(r'^do/synthesis/sbol/$', synthesis_view_sbol, name='shinthesis_view_sbol'),  # @IgnorePep8
    url(r'^do/synthesis/protocol/$', synthesis_view_protocol, name='synthesis_view_protocol'),  # @IgnorePep8
    url(r'^do/crispr/Single_Cas9_gRNA_Domesticator$', crispr_view_cas9_single, name='cas9_single_domestication'),
    url(r'^do/crispr/multi_cas9_gRNA_domesticator_1$', crispr_view_cas9_multiplexing, name='cas9_multiplexing_domestication'),
    url(r'^do/crispr/cas9_multiplexing/crispr_for_dummies/(?P<section>.+)?/?$', crispr_view_cas9_multiplexing_auto, name='cas9_multiplexing_domestication_auto'),
    url(r'^do/crispr/cas9_multiplexing/auto_add_tu/?$', crispr_view_cas9_multiplexing_auto_add_tu, name='cas9_multiplexing_domestication_auto_add_tu'),
    url(r'^do/cripsr/cas9_auto_omega_add/?$', cripsr_view_cas9_auto_omega_add, name='cripsr_view_cas9_auto_omega_add'),
    url(r'^do/crispr_cas9_auto_omega/protocol/', cripsr_view_cas9_auto_omega_protocol, name='cripsr_view_cas9_auto_omega_protocol'),
    url(r'^do/crisper_cas9_auto/protocol/?$', crispr_view_cas9_multiplexing_auto_protocol, name='cas9_multiplexing_auto_protocol'),
    url(r'^do/crispr/Single_Cas12a_gRNA_Domesticator', crispr_view_cas12_single, name='cas12_single_domestication'),
    url(r'^do/crispr/cas12_single/Single_Cas12a_gRNA_assembler$', crispr_view_cas12_single_TU, name='cas12_single_tu'),
    url(r'^do/crispr/multiple_Cas12a_gRNA_domesticator_2X$', crispr_view_cas12_multiplexing_2X, name='cas12_multiplexing_domestication_2X'),
    url(r'^do/crispr/multiple_Cas12a_gRNA_domesticator_3X$', crispr_view_cas12_multiplexing_3X, name='cas12_multiplexing_domestication_3X'),
    url(r'^do/crispr/multiple_Cas12a_gRNA_domesticator_4X$', crispr_view_cas12_multiplexing_4X, name='cas12_multiplexing_domestication_4X'),
    url(r'^do/crispr/multiple_Cas12a_gRNA_domesticator_5X$', crispr_view_cas12_multiplexing_5X, name='cas12_multiplexing_domestication_5X'),
    url(r'^do/crispr/multiple_Cas12a_gRNA_domesticator_6X$', crispr_view_cas12_multiplexing_6X, name='cas12_multiplexing_domestication_6X'),
    url(r'^do/crispr/cas12a_multiplexing/genbank$', cas12a_multiplexing_view_genbank, name='cas12_multiplexing_genbank'),
    url(r'^do/crispr/cas12a_multiplexing/sbol$', cas12a_multiplexing_view_sbol, name='cas12_multiplexing_sbol'),
    url(r'^do/crispr/cas12a_multiplexing/add$', cas12a_multiplexing_view_add, name='cas12_multiplexing_add'),
    url(r'^do/crispr/cas12a_multiplexing_view_protocol', cas12a_multiplexing_view_protocol, name="cas12a_protocol"),
    url(r'^do/crispr/cas12a_multiplexing/multiple_Cas12a_gRNA_assembler', crispr_view_cas12a_multiplexing_TU, name="cas12a_multiplex_tu"),
    url(r'^do/crispr/level0/add/$', level_0_add, name='crispr_level0_add'),
    url(r'^do/crispr/multi_cas9_gRNA_domesticator_2/(?P<section>.+)?/?$', level_0_editing, name='crispr_level0_editing'),
    url(r'^do/crispr/multi_regulatory_cas9_gRNA_domesticator_2/(?P<section>.+)?/?$', level_0_regulation, name='crispr_level0_regulation'),
    url(r'^do/crispr/multiplexing/(?P<section>.+)?/?$', crispr_multiplexing, name='crispr_multiplexing'),
    url(r'^do/crispr/level0/genbank/$', level_0_genbank, name='crispr_level0_genbank'),
    url(r'^do/crispr/level0/protocol/$', level_0_protocol, name='crispr_level0_protocol'),
    url(r'^do/crispr/level0/sbol/$', level_0_sbol, name='crispr_level0_sbol'),
    url(r'^do/crispr/protocol/$', crispr_view_protocol, name='crispr_view_add'),
    url(r'^do/crispr/add/$', crispr_view_add, name='crispr_view_add'),
    url(r'^fungal/do/multipartite/basic/$', multipartite_fungal_tu_view, name='fungal multipartite'),
    url(r'^fungal/do/gene_disruption/$', fungal_gene_disruption_view, name='fungal_gene_disruption_view'),
    url(r'^fungal/do/bipartite/(?P<form_num>.+)?/?$', fungal_bipartite_view, name='fungal_bipartite_view'),
    url(r'^add/vector/$', add_vector_view, name='add_vector'),
    url(r'^add/feature/$', add_feature_view, name='add_feature'),
    url(r'^add/experiment/(?P<exp_type>.+)/?$', add_experiment_view, name='add_experiment'),  # @IgnorePep8
    url(r'^search/experiment/$', search_experiment, name='search_experiments'),
    #url(r'^search/features/$', search_features_view, name='search_features'),
    url(r'^search/features/$', search_feature_index, name='search_features'),
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
    url(r'api/feature_uniquenames/$', feature_uniquenames, name='api_feature'),
    url(r'api/sbol/(?P<uniquename>.+)/$', feature_sbol, name='api_feature_sbol'),
    url(r'api/features_children/$', features_children, name='api_feature_children'),  # @IgnorePep8
    url('api/features_key_elements/$', features_key_elements, name='api_feature_key_elements'),  # @IgnorePep8
    url('api/excel_graph/(?P<excel_id>.+)?', excel_image, name='api_excel_image'),  # @IgnorePep8
    url('api/exp_keywords/$', experiment_keywords, name='api_exp_keywords'),  # @IgnorePep8
    url('api/excel_combined/(?P<uniquename>.+)?/(?P<exp_type>.+)?/.svg', combined_excel_image, name='api_combined_excel_images'),  # @IgnorePep8
    url(r'^bipart_autocomplete_1/$', BipartitePart1Autocomplete.as_view(), name='bipartite1_autocomplete'),  # @IgnorePep8,
    url(r'^bipart_autocomplete_2/$', BipartitePart2Autocomplete.as_view(), name='bipartite2_autocomplete'),  # @IgnorePep8
    url(r'^bipart_autocomplete_3/$', BipartitePart3Autocomplete.as_view(), name='bipartite3_autocomplete'),  # @IgnorePep8

    ]