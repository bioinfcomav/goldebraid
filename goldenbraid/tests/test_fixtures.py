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

from django.test import TestCase
from goldenbraid.models import Db

FIXTURES_TO_LOAD = ['initial_data2.json']


class TestFixtures(TestCase):
    'It tests that we can load the fixtures'
    fixtures = FIXTURES_TO_LOAD
    multi_db = True

    def test_initial(self):
        'It loads the fixtures.'
        assert Db.objects.all()
