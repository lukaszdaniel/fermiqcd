import csv
from optparse import OptionParser
import re
import os
import sys
import numpy as np
from scipy import stats
import wx
from matplotlib.axes import Subplot
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_ps import FigureCanvasPS
from matplotlib.backends.backend_wxagg import (
    FigureCanvasWxAgg as FigureCanvasWx,
    NavigationToolbar2WxAgg as NavigationToolbar2Wx,
)
from matplotlib.colors import ColorConverter
from matplotlib.colors import rgb2hex

import os.path
import zipfile
import datetime
import traceback

import ibootstrap
import ifit

# Script information
USAGE = "python iplot.py\n"
VERSION = (
    "iplot v1.0"
    "\n  Copyright (c) 2007 Massimo Di Pierro"
    "\n  All rights reserved"
    "\n  License: GPL 3.0"
    "\n\n  Written by Massimo Di Pierro <mdipierro@cs.depaul.edu>"
    "\n    and Vincent Harvey <vincent@vincentharvey.com>"
)
DESCRIPTION = "Plot the output of ibootstrap.py"

DEFAULT_OUTPUT_PREFIX = ".data/ibootstrap"


def clean(text):
    return text


def get_name(filename, prefix):
    return filename[:-4].replace(prefix, "").replace("_", " ")


def gen_plot(figure, plot_args):
    figure.clear()
    plot = plot_args[0](figure, *plot_args[1])
    figure.add_subplot(plot, 111)
    return plot


def csv_items(filename):
    return csv.reader(open(filename, "r"), delimiter=",", quoting=csv.QUOTE_NONNUMERIC)


# Notice that I defer the actual generation of the plots, and just pass
# along the function used to make them, along with some arguments. By convention
# it is assumed the function takes a matplotlib figure, then those arguments.
# Using a lambda would have been nicer, however due to the semantics of python
# loops, specifically writing over the loop variables, by the time the lambda
# was called many of its closure variables were overwritten


class IPlot:
    def __init__(self, filename, plot_type="ps", items=[], output_prefix=""):
        self.type = plot_type
        self.filename = filename
        self.output_prefix = output_prefix
        self.figure = Figure((8.5, 11))
        self.extension = ".png" if plot_type == "png" else ".ps"
        self.canvas = FigureCanvasAgg(self.figure) if plot_type == "png" else FigureCanvasPS(self.figure)
        self.plots = {}

        try:
            self.plot_raw_data(f"{filename}_raw_data.csv")
            self.plot_trails(f"{filename}_autocorrelations.csv")
            self.plot_trails(f"{filename}_trails.csv")
            self.plot_samples(f"{filename}_samples.csv")
            self.plot_min_mean_max(f"{filename}_min_mean_max.csv", items)
        except Exception as e:
            print(e)
        self.figure.clear()

    def print_all(self):
        def print_dict(d):
            for name, val in d.items():
                if isinstance(val, dict):
                    print_dict(val)
                else:
                    self.print_plot(name, val)

        print_dict(self.plots)

    def print_plot(self, name, plot_args):
        plot = gen_plot(self.figure, plot_args)
        fname = f"{self.output_prefix}{name}{self.extension}"
        self.canvas.print_figure(fname, 72)
        return fname

    def setup_plot(self, figure, name, xlabel, ylabel):
        plot = Subplot(figure, 111, axisbg="white", frameon="False")
        plot.set_title(name)
        plot.set_xlabel(xlabel)
        plot.set_ylabel(ylabel)
        return plot

    def make_error_bar(self, figure, name, xdata, ydata, errors, xlabel, ylabel, filename):
        plot = self.setup_plot(figure, name, xlabel, ylabel)
        plot.errorbar(xdata, ydata, errors, fmt="ko", markerfacecolor=None)
        plot.iplot_errorbar = True
        plot.filename = filename
        return plot

    def make_hist(self, figure, name, data, length, xlabel="x", ylabel="y"):
        plot = self.setup_plot(figure, name, xlabel, ylabel)
        plot.hist(data, length, facecolor="w", edgecolor="k")
        return plot

    def make_plot(self, figure, name, xdata, ydata, xlabel="x", ylabel="y"):
        plot = self.setup_plot(figure, name, xlabel, ylabel)
        plot.plot(xdata, ydata, "k-", markerfacecolor=None)
        return plot

    def plot_raw_data(self, filename):
        quants = {}
        prefix = get_name(filename, self.filename)
        self.plots[prefix] = quants
        for items in csv_items(filename):
            plots = {}
            tag = items[0]
            quants[f"{prefix} {tag}"] = plots

            data = items[1:]
            name = f"{prefix} {clean(tag)}"
            plots[name] = self.make_plot, (name, list(range(len(data))), data, "step", tag)
            name = f"{prefix} {clean(tag)} hist"
            plots[name] = self.make_hist, (name, data, len(data) // 20, tag, "frequency")

            mu = np.mean(data)
            sd = np.std(data)
            probs = [min(x, 1 - x) for x in [stats.norm.cdf((x - mu) / sd) for x in data]]
            name = f"{prefix} {clean(tag)} probability"
            plots[name] = self.make_plot, (name, list(range(len(probs))), probs, "step", f"probability {tag}")

    def plot_trails(self, filename):
        quants = {}
        prefix = get_name(filename, self.filename)
        self.plots[prefix] = quants
        for items in csv_items(filename):
            tag = items[0]
            data = items[1:]
            name = f"{prefix} {clean(tag)}"
            quants[name] = self.make_plot, (name, list(range(len(data))), data, "step", tag)

    def plot_samples(self, filename):
        quants = {}
        prefix = get_name(filename, self.filename)
        self.plots[prefix] = quants
        for items in csv_items(filename):
            tag = items[0]
            data = items[1:]
            name = f"{prefix} {clean(tag)} hist"
            quants[name] = self.make_hist, (name, data, len(data) // 10, tag, "frequency")

    def plot_min_mean_max(self, filename, xlab=None):
        if not xlab:
            xlab = [x.strip() for x in IPlot.indices.split(",")]
        lines = list(csv_items(filename))
        tags = lines[0]
        if not xlab or xlab[0] == "":
            xlab = [tags[1]]

        index = next((i - 1 for i in range(1, len(tags) - 3) if tags[i] == xlab[0]), -1)
        if index < 0:
            print("Error: Unable to find xlab tag:", xlab)
            raise Exception

        sets = {}
        for items in lines[1:]:
            tag = items[0]
            data = items[1:]
            legend = " ".join([f"{tags[i]}={data[i-1]}" for i in range(1, len(tags) - 3) if tags[i] not in xlab])
            sets.setdefault(legend, ([], [], [], []))
            x, y, yminus, yplus = sets[legend]
            t = data[index]
            x.append(t)
            y.append(data[-2])
            yminus.append(data[-3])
            yplus.append(data[-1])

        quants = {}
        prefix = get_name(filename, self.filename)
        self.plots[prefix] = quants
        for legend, (x, y, yminus, yplus) in sets.items():
            error_low = np.subtract(y, yminus)
            error_high = np.subtract(yplus, y)
            name = f"{prefix} {clean(legend)}"
            quants[name] = self.make_error_bar, (name, x, y, [error_low, error_high], tags[index + 1], tags[0], filename)



def make_gui():
    app = wx.App(False)
    root = wx.Frame(None, -1, "Iplot")

    root_sizer = wx.BoxSizer(wx.VERTICAL)
    root.SetSizer(root_sizer)

    log_widget = wx.TextCtrl(root, -1, style=wx.TE_MULTILINE, size=(600, 600))
    root_sizer.Add(log_widget, 1, wx.ALL | wx.EXPAND)

    clear_id = wx.Window.NewControlId()
    clear_button = wx.Button(root, clear_id, "Clear Log Window")
    root.Bind(wx.EVT_BUTTON, lambda e: log_widget.Clear(), id=clear_id)
    root_sizer.Add(clear_button, 0, wx.RIGHT | wx.EXPAND)

    menubar = wx.MenuBar(style=wx.MENU_TEAROFF)
    file_menu = wx.Menu(style=wx.MENU_TEAROFF)

    file_new_id = wx.Window.NewControlId()
    file_menu.Append(file_new_id, "&New")
    root.Bind(wx.EVT_MENU, lambda e: show_bootstrap_init(), id=file_new_id)

    file_save_id = wx.Window.NewControlId()
    file_menu.Append(file_save_id, "&Save")
    root.Bind(wx.EVT_MENU, lambda e: save_data(), id=file_save_id)

    file_open_id = wx.Window.NewControlId()
    file_menu.Append(file_open_id, "Open")
    root.Bind(wx.EVT_MENU, lambda e: open_data(), id=file_open_id)

    file_exit_id = wx.Window.NewControlId()
    file_menu.Append(file_exit_id, "&Exit")
    root.Bind(wx.EVT_MENU, lambda e: sys.exit(0), id=file_exit_id)

    menubar.Append(file_menu, "&File")
    root.SetMenuBar(menubar)

    def save_data():
        if not os.access(".data", os.R_OK):
            dlg = wx.MessageDialog(None, "No data to save", style=wx.ICON_ERROR | wx.OK)
            dlg.ShowModal()
            return
        fdialog = wx.FileDialog(None, "Choose save file", os.getcwd(), style=wx.SAVE, wildcard="*.zip")
        if fdialog.ShowModal() == wx.ID_OK:
            zipfilename = fdialog.GetFilename()
            with zipfile.ZipFile(f"{zipfilename.rstrip('.zip')}.zip", "w") as zip_file:
                for filename in os.listdir(".data"):
                    zip_file.write(filename)
        fdialog.Destroy()

    def open_data():
        log_msg("SORRY NOT IMPLEMENTED")
        return
        setup_data_dir()
        fdialog = wx.FileDialog(
            None,
            "Choose zipped file to open",
            os.getcwd(),
            style=wx.OPEN,
            wildcard="*.zip",
        )
        if fdialog.ShowModal() == wx.ID_OK:
            zipfilename = fdialog.GetFilename()
            zip = zipfile.ZipFile(str(zipfilename).rstrip(".zip") + ".zip", "r")
            for filename in zip.namelist():
                # os.makedirs(os.path.dirname(filename))
                file(".data/" + os.path.basename(filename), "wb").write(
                    zip.read(filename)
                )
            zip.close()
            load_iplot(DEFAULT_OUTPUT_PREFIX)
        fdialog.Destroy()

    def log_msg(msg):
        log_widget.AppendText(f"{datetime.datetime.now()}: {msg}\n")

    def to_wx_color(string_color):
        rgb = ColorConverter().to_rgb(string_color)
        return wx.Color(*[x * 255 for x in rgb])

    def to_mpl_color(wxcolor):
        return rgb2hex([float(x) / 255.0 for x in wxcolor.Get()])

    def show_plot(name, plot_args):
        pwin = wx.Frame(root, -1, name)

        # iplot.figure.clear() #delaxes(plot)
        figure = Figure(facecolor="#ffffff")
        plot = [gen_plot(figure, plot_args)]

        # iplot.figure.clear()
        # iplot.figure.add_axes(plot)

        frame_box = wx.BoxSizer(wx.VERTICAL)
        pwin.SetSizer(frame_box)

        graph_panel = wx.Panel(pwin, -1, style=wx.SUNKEN_BORDER)
        graph_vbox = wx.BoxSizer(wx.VERTICAL)
        graph_panel.SetSizer(graph_vbox)

        canvas = FigureCanvasWx(graph_panel, -1, figure)
        canvas.draw()
        graph_vbox.Add(canvas, 1, wx.ALL | wx.EXPAND)

        toolbar = NavigationToolbar2Wx(canvas)
        toolbar.Realize()
        graph_vbox.Add(toolbar, 0, wx.LEFT | wx.EXPAND)
        toolbar.update()

        edit_panel = wx.Panel(pwin, -1, style=wx.SUNKEN_BORDER)
        edit_box = wx.GridSizer(4, 2)
        edit_panel.SetSizer(edit_box)

        grid_items = []

        def make_entry(name, default=""):
            grid_items.append((wx.StaticText(edit_panel, -1, name), 0, 0))
            id = wx.Window.NewControlId()
            entry = wx.TextCtrl(edit_panel, id, default)  # , size=(150, -1))
            grid_items.append((entry, 1, wx.RIGHT | wx.EXPAND))
            return entry

        def set_attributes(*args):
            try:
                if hasattr(plot[0], "iplot_errorbar"):
                    loc = {}
                    try:
                        expression, initial = expr_entry.GetValue().split("@")
                        symbols, points = ifit.read_min_mean_max_file(plot[0].filename)
                        fit = ifit.IFit(
                            expression, points, symbols, condition=cond_entry.GetValue()
                        )
                        variables = ifit.restricted_eval("dict(%s)" % initial, loc)
                        variables, least_squares, hessian = fit.fit(**variables)
                        fit.scatter_points = 0
                        if scatter_entry.GetValue():
                            fit.scatter_points = int(scatter_entry.GetValue())
                        if fit.scatter_points:
                            fit.iterative_fit(
                                DEFAULT_OUTPUT_PREFIX + "_scatter.csv", **variables
                            )
                        for key, value in list(variables.items()):
                            log_msg("%s = %g" % (key, value))
                        log_msg("least_squares=%s" % least_squares)
                        log_msg("hessian=%s" % str(hessian))
                        bound_min, bound_max = (
                            plot[0].get_xaxis().get_data_interval().get_bounds()
                        )
                        if bound_min == bound_max:
                            return
                        new_ydata = []
                        new_xdata = numpy.arange(
                            bound_min,
                            bound_max,
                            (bound_max - bound_min) / 200.0,
                            numpy.float,
                        )
                        # str() to convert from unicode
                        # extrapolate_variable = str(extrap_entry.GetValue())
                        # if fit.scatter_points:
                        #    extrap = lambda f, c: fit.extrapolate_with_errors(**c)
                        # else:
                        #    extrap = lambda f, c: fit.extrapolate(**c)
                        # for item in new_xdata:
                        #    coordinates = {extrapolate_variable:item}
                        #    #coordinates=ifit.restricted_eval('dict(%s)' % item,loc)
                        #    e = extrap(fit, coordinates)
                        #    new_ydata.append(e)
                        #    print 'extrapolation %s -> %s' % (item,str(e))
                        #### replaced by
                        for item in new_xdata:
                            new_ydata.append(fit.f(**{IPlot.indices[0]: item}))
                        plot[0] = gen_plot(figure, plot_args)
                        if fit.scatter_points:
                            low, mid, high = [], [], []
                            for h, m, l in new_ydata:
                                low.append(l)
                                mid.append(m)
                                high.append(h)
                            plot[0].plot(new_xdata, low, color="#ff5555")
                            plot[0].plot(new_xdata, mid, color="#ff0000")
                            plot[0].plot(new_xdata, high, color="#bb0000")
                        else:
                            plot[0].plot(new_xdata, new_ydata, color="r")
                    except Exception as e:
                        log_msg(e)
                        log_msg("NO FIT")
                plot[0].set_xlabel(xlabel_entry.GetValue())
                plot[0].set_ylabel(ylabel_entry.GetValue())
                plot[0].set_title(title_entry.GetValue())
                plot[0].set_axis_bgcolor(background_entry.GetValue())
                figure.set_facecolor(background_entry.GetValue())
                canvas.draw()
            except Exception as e:
                traceback.print_tb(sys.exc_info()[2])
                dlg = wx.MessageDialog(
                    None,
                    "Error setting preferences:\t" + str(e),
                    style=wx.ICON_ERROR | wx.OK,
                )
                dlg.ShowModal()

        expr_entry = cond_entry = extrap_entry = scatter_entry = None
        if hasattr(plot[0], "iplot_errorbar"):
            expr_entry = make_entry("Expression", "a*t*t + b*t + c@a=1.0,b=1.0,c=1.0")
            cond_entry = make_entry("Condition", "True")
            # extrap_entry = make_entry('Extrapolate to', 't')
            scatter_entry = make_entry("Number of scatter points", "0")

        xlabel_entry = make_entry("X-Label", plot[0].xaxis.label.get_text())
        ylabel_entry = make_entry("Y-Label", plot[0].yaxis.label.get_text())
        title_entry = make_entry("Title", plot[0].title.get_text())
        background_entry = make_entry("Background Color", plot[0].get_axis_bgcolor())

        def color_dialog(*args):
            try:
                wxcolor = to_wx_color(background_entry.GetValue())
            except:
                wxcolor = to_wx_color(plot[0].get_axis_bgcolor())
            cdialog = wx.ColourDialog(None)
            cdialog.GetColourData().SetColour(wxcolor)
            if cdialog.ShowModal() == wx.ID_OK:
                rgb = to_mpl_color(cdialog.GetColourData().GetColour())
                background_entry.SetValue(rgb)

        id = wx.Window.NewControlId()
        button = wx.Button(edit_panel, id, "Color Select")
        pwin.Bind(wx.EVT_BUTTON, color_dialog, id=id)
        grid_items.append(button)

        id = wx.Window.NewControlId()
        button = wx.Button(edit_panel, id, "Set")
        pwin.Bind(wx.EVT_BUTTON, set_attributes, id=id)
        grid_items.append(button)
        edit_box.AddMany(grid_items)

        frame_box.Add(edit_panel, 0, wx.RIGHT | wx.EXPAND, 5)
        frame_box.Add(graph_panel, 1, wx.ALL | wx.EXPAND, 5)

        pwin.Show(True)
        pwin.Fit()

    def bind_plot_menu(name, value, id):
        root.Bind(wx.EVT_MENU, lambda e: show_plot(name, value), id=id)

    def sorted(d):
        items = list(d.items())
        items.sort()
        return items

    def set_menu(menu, d):
        for name, val in sorted(d):
            id = wx.Window.NewControlId()
            if isinstance(val, wx.Menu):
                menu.AppendMenu(id, name, val)
            else:
                menu.Append(id, name)
                bind_plot_menu(name, val, id)
                # root.Bind(wx.EVT_MENU, lambda e: show_plot(name, val), id=id)

    def make_dict_menu(d):
        items = {}
        for name, val in list(d.items()):
            if isinstance(val, dict):
                submenu_items = make_dict_menu(val)
                submenu = wx.Menu(style=wx.MENU_TEAROFF)
                set_menu(submenu, submenu_items)
                items[name] = submenu
            else:
                items[name] = val
        return items

    def show_bootstrap_init():
        frame = wx.Frame(root, -1, "Select Bootstrap parameters")
        frame_box = wx.BoxSizer(wx.VERTICAL)
        frame.SetSizer(frame_box)

        param_panel = wx.Panel(frame, -1, style=wx.SUNKEN_BORDER)
        parameters = [
            ("filepattern", "File Pattern", str, "*.log"),
            ("expression", "Expression", str, '"a[<t>]"/"b[<t>]"'),
            ("condition", "Condition", str, "True"),
            ("indices", "Indices", str, "t"),
            ("min_index", "Minimum Index", int, 0),
            ("max_index", "Maximum Index", int, 0),
            ("nsamples", "Number of samples", int, 100),
            ("percent", "Bootstrap percent", float, 0.158),
            ("raw", "Raw", bool, False),
            ("import_module", "Import module", str, None)
        ]

        def empty_if_none(e):
            return "" if e is None else str(e)

        psizer = wx.GridBagSizer(5, 10)
        entries = []
        for index, (pname, name, ptype, default) in enumerate(parameters):
            psizer.Add(wx.StaticText(param_panel, -1, name), (index, 0), wx.DefaultSpan, wx.EXPAND)
            entry = wx.CheckBox(param_panel, -1) if ptype == bool else wx.TextCtrl(param_panel, -1, empty_if_none(default), size=wx.Size(500, 20))
            entries.append(entry)
            psizer.Add(entry, (index, 1), wx.DefaultSpan, wx.RIGHT | wx.EXPAND)

        def bootstrap_run():
            param_dict = {pname: ptype(entry.GetValue()) for pname, entry, ptype in zip([p[0] for p in parameters], entries, [p[2] for p in parameters])}
            param_dict["output_prefix"] = "DEFAULT_OUTPUT_PREFIX"
            try:
                print(param_dict)
                bootstrap = ibootstrap.IBootstrap(**param_dict)
                for msg in bootstrap.report:
                    log_msg("IBootstrap MSG:\t" + msg)
                log_msg("IBootstrap status:\t" + bootstrap.status)
                load_iplot(DEFAULT_OUTPUT_PREFIX)
                frame.Destroy()
            except Exception as e:
                log_msg(str(e))
                traceback.print_exc(file=sys.stdout)

        param_panel.SetSizer(psizer)
        param_panel.Fit()
        frame_box.Add(param_panel, 1, wx.ALL | wx.EXPAND)

        button_panel = wx.Panel(frame, -1)
        button_panel_sizer = wx.BoxSizer(wx.HORIZONTAL)

        button = wx.Button(button_panel, wx.ID_CANCEL)
        frame.Bind(wx.EVT_BUTTON, lambda e: frame.Destroy(), id=wx.ID_CANCEL)
        button_panel_sizer.Add(button, 1, wx.LEFT | wx.EXPAND)

        button = wx.Button(button_panel, wx.ID_OK)
        frame.Bind(wx.EVT_BUTTON, lambda e: bootstrap_run(), id=wx.ID_OK)
        button_panel_sizer.Add(button, 1, wx.LEFT | wx.EXPAND)

        button_panel.SetSizer(button_panel_sizer)
        frame_box.Add(button_panel, 0, wx.ALIGN_BOTTOM | wx.EXPAND)

        frame.Fit()
        frame.Centre()
        frame.Show(True)

    def setup_data_dir():
        try:
            os.mkdir(".data")
        except Exception as e:
            log_msg(str(e))

    def load_iplot(path_prefix):
        iplot = IPlot(path_prefix)
        set_iplot_menu(iplot)

    def set_iplot_menu(iplot):
        if menubar.GetMenuCount() > 1 and menubar.GetMenu(1) != None:
            menubar.Remove(1).Destroy()
        iplot_items = make_dict_menu(iplot.plots)
        iplot_menu = wx.Menu(style=wx.MENU_TEAROFF)
        set_menu(iplot_menu, iplot_items)
        menubar.Append(iplot_menu, "&Plots")

    root.Fit()
    root.Show(True)

    # load_iplot(DEFAULT_OUTPUT_PREFIX)

    app.MainLoop()


def shell_iplot():
    parser = OptionParser(usage=USAGE, version=VERSION)
    parser.description = DESCRIPTION
    parser.add_option(
        "-d",
        "--dest_prefix",
        default="",
        dest="dest_prefix",
        help="the prefix used to build output filenames",
    )
    parser.add_option(
        "-o",
        "--origin_prefix",
        default="ibootstrap",
        dest="origin_prefix",
        help="the prefix used to build input filenames",
    )
    parser.add_option(
        "-p", "--plot_type", default="ps", dest="plot_type", help="ps or png"
    )
    parser.add_option(
        "-v",
        "--plot_variables",
        default="",
        dest="plot_variables",
        help="plotting variables",
    )
    parser.add_option(
        "-b",
        "--batch",
        default=False,
        action="store_true",
        dest="batch",
        help="Batch print all plots",
    )

    (options, args) = parser.parse_args()
    if options.batch:
        iplot = IPlot(
            options.origin_prefix,
            options.plot_type,
            options.plot_variables.split(","),
            output_prefix=options.dest_prefix,
        )
        iplot.print_all()
    else:
        make_gui()


if __name__ == "__main__":
    shell_iplot()
