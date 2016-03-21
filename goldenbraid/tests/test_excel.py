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
                               draw_columns, draw_scatter,
                               draw_combined_graph)
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

#         try:
#             parse_xlsx('/home/peio/tomato_gbs.xlsx')
#             self.fail('runtimeerror expected')
#         except RuntimeError:
#             pass

    def test_draw_columns(self):
        _, labels, data = parse_xlsx(os.path.join(TEST_DATA, 'columns.xlsx'))
        out_fhand = NamedTemporaryFile(suffix='.svg')
        draw_columns(labels, data, out_fhand)
        raw_input(out_fhand.name)

    def test_draw_1column(self):
        _, labels, data = parse_xlsx(os.path.join(TEST_DATA, '1column.xlsx'))
        out_fhand = NamedTemporaryFile(suffix='.svg')
        draw_columns(labels, data, out_fhand)
        raw_input(out_fhand.name)

    def test_draw_3columns(self):
        _, labels, data = parse_xlsx(os.path.join(TEST_DATA, '3columns.xlsx'))
        out_fhand = NamedTemporaryFile(suffix='.svg')
        draw_columns(labels, data, out_fhand)
        raw_input(out_fhand.name)

    def test_draw_10columns(self):
        _, labels, data = parse_xlsx(os.path.join(TEST_DATA, '10columns.xlsx'))
        out_fhand = NamedTemporaryFile(suffix='.svg')
        draw_columns(labels, data, out_fhand)
        raw_input(out_fhand.name)

    def test_draw_scatter(self):
        _, labels, data = parse_xlsx(os.path.join(TEST_DATA, 'scatter.xlsx'))
        out_fhand = NamedTemporaryFile(suffix='.svg')
        draw_scatter(labels, data, out_fhand)
        raw_input(out_fhand.name)

    def test_empty_file_excel(self):
        _, labels, data = parse_xlsx(os.path.join(TEST_DATA, 'empty_file.xlsx'))


COMBINED_DATA = {u'GB_EXP_69': ({u'X-label': u'fill in x axis label',
                                 u'Y-label': u'fill in y axis label',
                                 u'title': u'title GB_EXP_69'},
                                {u'Y-values': [1.0, 3.0, 5.0, 5.0, 3.0, 7.0],
                                 u'Y-stdev': [0.1, 0.1, 0.1, 0.1, 0.1, 0.2],
                                 u'X-values': ['time1', 'time2', 'time4',
                                               'time5', 'time6', 'time7']}),
                 u'GB_EXP_68': ({u'X-label': u'fill in x axis label',
                                 u'Y-label': u'fill in y axis label',
                                 u'title': u'title GB_EXP_68'},
                                {u'Y-values': [4.0, 5.0, 5.0],
                                 u'Y-stdev': [0.1, 0.1, 0.1],
                                 u'X-values': ['time3', 'time4', 'time5']}),
                 u'GB_EXP_52': ({u'X-label': u'xlabel2',
                                 u'Y-label': u'ylabel2',
                                 u'title': u'title GB_EXP_52'},
                                {u'Y-values': [5.0],
                                 u'Y-stdev': [0.1],
                                 u'X-values': ['time3']}),
                 u'GB_EXP_64': ({u'X-label': u'fill in x axis label',
                                 u'Y-label': u'fill in y axis label',
                                 u'title': u'title GB_EXP_64'},
                                {u'Y-values': [2.0, 3.0, 4.0, 5.0, 5.0],
                                 u'Y-stdev': [0.1, 0.1, 0.1, 0.1, 0.1],
                                 u'X-values': ['time1', 'time2', 'time3',
                                               'time4', 'time5']})
                 }


class TestExcelCombined(TestCase):
    def test_combined(self):
        fhand = NamedTemporaryFile(suffix='.svg')
        fhand = open('/tmp/tmpf.svg', 'w')
        draw_combined_graph(COMBINED_DATA, fhand)
        # raw_input(fhand.name)
        fhand.close()


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'TestExcelCombined']
    main()
