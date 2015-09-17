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

from __future__ import division

from openpyxl import load_workbook
from openpyxl.utils import column_index_from_string

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

FIGURE_SIZE = (15.0, 11.0)  # inche

KEY_CELL_NAME = 'plot_type'
COLUMNS = 'columns'
SCATTER = 'scatter'
XLABEL = 'X-label'
YLABEL = 'Y-label'
TITLE = 'title'
YVALUES = 'Y-values'
XVALUES = 'X-values'
YSTDEV = 'Y-stdev'
XSTDEV = 'X-stdev'

LABELS = (XLABEL, YLABEL, TITLE)
COLUMNS_DATA_TYPE = (YVALUES, XVALUES, YSTDEV)
SCATTER_DATA_TYPE = COLUMNS_DATA_TYPE + (XSTDEV,)


def search_data(ws):
    row_start = None
    row_end = None
    column_start = None
    for row in ws.rows:
        for cell in row:
            cell_content = cell.value
            column_index = cell.column
            row_index = cell.row
            if cell_content == KEY_CELL_NAME:
                row_start = row_index
                column_start = column_index
            if (row_start is not None and column_start == cell.column and
                    cell_content is None):

                row_end = cell.row - 1
                if column_start is None or row_start is None:
                    raise RuntimeError('No data found in excel')
                return (row_start, row_end,
                        column_index_from_string(column_start))
    else:
        row_end = cell.row
        if column_start is None or row_start is None:
            raise RuntimeError('No data found in excel')
        return row_start, row_end, column_index_from_string(column_start)


def parse_xlsx(fpath_or_fhand):

    wb = load_workbook(filename=fpath_or_fhand)
    ws = wb.active
    row_start, row_end, column_start = search_data(ws)
    plot_type = ws.cell(row=row_start, column=column_start + 1).value
    labels = {}
    next_row = row_start + 1
    for row_index in range(next_row, next_row + len(LABELS)):
        key = ws.cell(row=row_index, column=column_start).value
        value = ws.cell(row=row_index, column=column_start + 1).value
        if key not in LABELS:
            raise RuntimeError('Given excel is malformed')
        labels[key] = value
    next_row = row_index + 1

    if plot_type == COLUMNS:
        header_items = COLUMNS_DATA_TYPE
    elif plot_type == SCATTER:
        header_items = SCATTER_DATA_TYPE
    else:
        raise RuntimeError('Can not detect plot type')

    data = {}
    ordered_header = []
    for row_index in range(next_row, row_end + 1):
        for column_index in range(column_start,
                                  column_start + len(header_items)):
            cell_value = ws.cell(row=row_index, column=column_index).value
            if row_index == next_row:
                if cell_value not in header_items:
                    msg = '{} not a valid data type name'
                    raise RuntimeError(msg.format(cell_value))
                data[cell_value] = []
                ordered_header.append(cell_value)
            else:
                data_type = ordered_header[column_index - column_start]
                if data_type in (XSTDEV, YSTDEV, YVALUES):
                    value = float(cell_value)
                elif plot_type == SCATTER and data_type == XVALUES:
                    value = float(cell_value)
                elif plot_type == COLUMNS and data_type == XVALUES:
                    value = str(cell_value)

                data[data_type].append(value)
    # print plot_type, labels, data
    return plot_type, labels, data


def draw_columns(labels, data, out_fhand):
    canvas, axes = get_canvas_and_axes()
    x_vals = data[XVALUES]
    y_vals = data[YVALUES]
    y_stdev = data.get(YSTDEV, None)
    bar_width = 0.8

    half_width = bar_width / 2
    xlabel_pos = list(range(len(x_vals)))
    left = [i - half_width for i in xlabel_pos]
    bottom = [0] * len(x_vals)
    height = y_vals
    kwargs = {}
    if y_stdev:
        kwargs['yerr'] = y_stdev
        kwargs['ecolor'] = 'red'
        # kwargs['capsize'] = 8
    axes.bar(left=left, width=bar_width, bottom=bottom, height=height,
             **kwargs)
    axes.set_title(labels[TITLE])
    axes.set_ylabel(labels[YLABEL])
    axes.set_xlabel(labels[XLABEL], clip_on=False)
    axes.set_xticklabels([''] + data[XVALUES])
    canvas.print_figure(out_fhand, plot_format='png')


def draw_scatter(labels, data, out_fhand):
    canvas, axes = get_canvas_and_axes()
    x_vals = data[XVALUES]
    y_vals = data[YVALUES]
    y_stdev = data.get(YSTDEV, None)
    x_stdev = data.get(XSTDEV, None)
    # axes.scatter(x_vals, y_vals)
    axes.errorbar(x_vals, y_vals, xerr=x_stdev, yerr=y_stdev, fmt='o')

    axes.set_title(labels[TITLE])
    axes.set_ylabel(labels[YLABEL])
    axes.set_xlabel(labels[XLABEL], clip_on=False)
    canvas.print_figure(out_fhand, plot_format='png')


def plot_from_excel(fpath_or_fhand, out_fhand):
    plot_type, labels, data = parse_xlsx(fpath_or_fhand)
    if plot_type == 'columns':
        draw_columns(labels, data, out_fhand)
    elif plot_type == 'scatter':
        draw_scatter(labels, data, out_fhand)
    else:
        raise RuntimeError()


def get_canvas_and_axes():
    'It returns a matplotlib canvas and axes instance'
    fig = Figure()
    canvas = FigureCanvas(fig)
    axes = fig.add_subplot(111)
    # fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    return canvas, axes
