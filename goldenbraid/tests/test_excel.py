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

import os.path
from django.test import TestCase
from unittest import main

from goldenbraid.tests.test_models import TEST_DATA
from goldenbraid.excel import (parse_xlsx, SCATTER, COLUMNS, YVALUES, XLABEL,
    draw_columns, XVALUES, draw_scatter)
from tempfile import NamedTemporaryFile


class TestExcel(TestCase):
    def test_parse_xlsx(self):
        plot_type, labels, data = parse_xlsx(os.path.join(TEST_DATA,
                                                          'columns.xlsx'))
        assert plot_type == COLUMNS
        assert data[YVALUES] == [1.0, 3.0, 4.0, 5.0, 5.0]
        assert labels[XLABEL] == 'xlabel2'
        plot_type, labels, data = parse_xlsx(os.path.join(TEST_DATA,
                                                          'scatter.xlsx'))
        plot_type = SCATTER
        assert data[YVALUES] == [1.0, 3.0, 4.0, 5.0, 5.0]

        assert labels[XLABEL] == 'xlabel1'

        try:
            parse_xlsx('/home/peio/tomato_gbs.xlsx')
            self.fail('runtimeerror expected')
        except RuntimeError:
            pass

    def test_draw_columns(self):
        _, labels, data = parse_xlsx(os.path.join(TEST_DATA, 'columns.xlsx'))
        out_fhand = NamedTemporaryFile(suffix='.png')
        draw_columns(labels, data, out_fhand)
        # raw_input(out_fhand.name)

    def test_draw_scatter(self):
        _, labels, data = parse_xlsx(os.path.join(TEST_DATA, 'scatter.xlsx'))
        out_fhand = NamedTemporaryFile(suffix='.png')
        draw_scatter(labels, data, out_fhand)
        # raw_input(out_fhand.name)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'TestExcel.test_draw_scatter']
    main()
