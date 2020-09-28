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

from django.conf.urls import include, url
from django.conf import settings
from django.conf.urls.static import static
from django.contrib.auth.views import LoginView
# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = [
    # Examples:
    # url(r'^$', 'dj_goldenbraid.views.home', name='home'),
    # url(r'^dj_goldenbraid/', include('dj_goldenbraid.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    url(r'^accounts/login/$', LoginView.as_view()),
    # Uncomment the next line to enable the admin:
    url(r'^admin/', admin.site.urls),
    url(r'^genome_domestication/search/', include('restcmd_client.urls')),
    url(r'^genome_domestication/', include('gb_genome_domestication.urls')),
    # (r'^cmd_client/', include('restcmd_client.urls')),

    url(r'', include('goldenbraid.urls')),
]
urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
#     static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT),
