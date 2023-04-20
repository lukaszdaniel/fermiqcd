from optparse import *
import numpy
from scipy import stats
import tkinter as Tk
from matplotlib.axes import Subplot
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_ps import FigureCanvasPS
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg as FigureCanvasTk,
    NavigationToolbar2TkAgg as NavigationToolbar2Tk,
)

# from rpy import *
import re

usage = "python iplot.py\n"
version = ("iplotv1.0"
           "\n  Copyright (c) 2007 Massimo Di Pierro"
           "\n  All rights reserved"
           "\n  License: GPL 3.0"
           "\n\n  Written by Massimo Di Pierro <mdipierro@cs.depaul.edu>")

description = "plot the output of ibootstrap.py"

# r.library('Hmisc')


def clean(text):
    return re.sub("\W", "", text)


def gen_plot(figure, plot_args):
    figure.clear()
    plot = plot_args[0](figure, *plot_args[1])
    figure.add_subplot(plot, 111)


# Notice that I defer the actual generation of the plots, and just pass
# along the function used to make them, along with some arguments. By convention
# it is assumed the function takes a matplotlib figure, then those arguments.
# Using a lambda would have been nicer, however due to the semantics of python
# loops, specifically writing over the loop variables, by the time the lambda
# was called many of its closure variables were overwritten


class IPlot:

    def __init__(self, filename, plot_type, items=[], output_prefix=""):
        self.type = plot_type

        self.output_prefix = output_prefix
        self.figure = Figure((8.5, 11))
        # self.figure = Figure((8.3,11.7))
        if plot_type == "png":
            self.canvas = FigureCanvasAgg(self.figure)
            self.extension = ".png"
        else:
            self.canvas = FigureCanvasPS(self.figure)
            self.extension = ".ps"
        self.plots = {}

        # if self.type=='quartz': r.quartz()
        self.plot_raw_data(filename + "_raw_data.csv")
        self.plot_trails(filename + "_trails.csv")
        self.plot_samples(filename + "_samples.csv")
        self.plot_min_mean_max(filename + "_min_mean_max.csv", items)
        self.figure.clear()

    def print_all(self):

        def print_dict(d):
            for name, val in list(d.items()):
                if isinstance(val, dict):
                    print_dict(val)
                else:
                    self.print_plot(name, val)

        print_dict(self.plots)

    def print_plot(self, name, plot_args):
        plot = gen_plot(self.figure, plot_args)
        fname = self.output_prefix + name + self.extension
        self.canvas.print_figure(fname, 72)
        return fname

    def setup_plot(self, figure, name, xlabel, ylabel):
        plot = Subplot(figure, 111)
        plot.set_title(name)
        plot.set_xlabel(xlabel)
        plot.set_ylabel(ylabel)
        return plot

    def make_error_bar(self, figure, name, xdata, ydata, errors, xlabel,
                       ylabel):
        plot = self.setup_plot(figure, name, xlabel, ylabel)
        plot.errorbar(xdata, ydata, errors, fmt="ro", markerfacecolor=None)
        return plot

    def make_hist(self,
                  figure,
                  name,
                  data,
                  length,
                  xlabel="x",
                  ylabel="y"):  # prob='T'
        # be sure to add a 'rug', whatever that is
        plot = self.setup_plot(figure, name, xlabel, ylabel)
        plot.hist(data, length)
        return plot
        # plot.set_xticks(data) # was not a good 'rug'

    def make_plot(self,
                  figure,
                  name,
                  xdata,
                  ydata,
                  xlabel="x",
                  ylabel="y"):  # type='p'
        plot = self.setup_plot(figure, name, xlabel, ylabel)
        plot.plot(xdata, ydata, "ko", markerfacecolor=None)
        return plot

    def plot_raw_data(self, filename):
        quants = {}
        self.plots[filename[:-4]] = quants
        for line in open(filename, "r"):
            plots = {}
            items = line.split(",")
            tag = items[0][1:-1]

            quants[filename[:-4] + "_" + tag] = plots

            data = [float(x) for x in items[1:]]
            name = filename[:-4] + "_%s" % clean(tag)
            plots[name] = self.make_plot, (
                name,
                list(range(len(data))),
                data,
                "step",
                tag,
            )
            name = filename[:-4] + "_%s_hist" % clean(tag)
            plots[name] = self.make_hist, (name, data, len(data) / 20, tag,
                                           "frequency")
            # self.begin(filename[:-4]+'_%s_qq.ps' % clean(tag))
            # r.qqnorm(data,xlab=tag+' quantiles',main='')
            # r.qqline(data)
            # self.end()
            probs = []
            mu = numpy.mean(data)
            sd = numpy.std(data)
            probs = [
                min(x, 1 - x)
                for x in [stats.norm.cdf((x - mu) / sd) for x in data]
            ]
            name = filename[:-4] + "_%s_probability" % clean(tag)
            plots[name] = self.make_plot, (
                name,
                list(range(len(probs))),
                probs,
                "step",
                "probability " + tag,
            )

    def plot_trails(self, filename):
        quants = {}
        self.plots[filename[:-4]] = quants
        for line in open(filename, "r"):
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(x) for x in items[1:]]
            name = filename[:-4] + "_%s" % clean(tag)
            quants[name] = self.make_plot, (
                name,
                list(range(len(data))),
                data,
                "step",
                tag,
            )

    def plot_samples(self, filename):
        quants = {}
        self.plots[filename[:-4]] = quants
        for line in open(filename, "r"):
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(x) for x in items[1:]]
            name = filename[:-4] + "_%s_hist" % clean(tag)
            quants[name] = self.make_hist, (
                name,
                data,
                len(data) / 10,
                tag,
                "frequency",
            )
            # self.begin(filename[:-4]+'_%s_qq.ps' % clean(tag))
            # r.qqnorm(data,xlab=tag+' quantiles',main='')
            # r.qqline(data)
            # self.end()

    def plot_min_mean_max(self, filename, xlab=["t"]):
        lines = open(filename, "r").readlines()
        tags = [x[1:-1] for x in re.compile('"[^"]*"').findall(lines[0])]
        if not xlab or xlab[0] == "":
            xlab = [tags[1]]
        index = -1
        for i in range(1, len(tags) - 3):
            if tags[i] == xlab[0]:
                index = i - 1
        if index < 0:
            print("error", xlab)
            raise Exception
        sets = {}
        for line in lines[1:]:
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(p) for p in items[1:]]
            legend = ""
            for i in range(1, len(tags) - 3):
                if not tags[i] in xlab:
                    legend += "%s=%g " % (tags[i], data[i - 1])
            if legend not in sets:
                x, y, yminus, yplus = [], [], [], []
                sets[legend] = (x, y, yminus, yplus)
            else:
                x, y, yminus, yplus = sets[legend]
            t = data[index]
            x.append(t)
            y.append(data[-2])
            yminus.append(data[-3])
            yplus.append(data[-1])
        # v=r.FALSE
        quants = {}
        self.plots[filename[:-4]] = quants
        for legend in list(sets.keys()):
            x, y, yminus, yplus = sets[legend]
            # matplotlib takes offsets for errors, not absolute y positions
            error_low = numpy.subtract(y, yminus)
            error_high = numpy.subtract(yplus, y)
            name = filename[:-4] + "_%s" % clean(legend)
            quants[name] = self.make_error_bar, (
                name,
                x,
                y,
                [error_low, error_high],
                tags[index + 1],
                tags[0],
            )
            # v=r.TRUE


def make_menu(iplot):
    # win = Tk.Tk()
    # win.wm_title("IPlot")

    main_window = [None]

    def show_plot(name, plot_args):
        pwin = Tk.Toplevel()
        pwin.wm_title(name)

        # iplot.figure.clear() #delaxes(plot)
        figure = Figure()
        plot = gen_plot(figure, plot_args)

        # iplot.figure.clear()
        # iplot.figure.add_axes(plot)

        canvas = FigureCanvasTk(figure, master=pwin)
        canvas.draw()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, pwin)
        toolbar.update()
        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        # def print_plot(*args):
        #    fname = iplot.print_plot(name, plot_args)
        #    top = Tk.Toplevel(pwin)
        #    Tk.Label(top, text="Printed to file '%s'" % fname).pack()
        #    Tk.Button(top, text="OK", command=lambda: top.destroy()).pack()
        #    pwin.wait_window(top)

        # button = Tk.Button(pwin, text="Print", command=print_plot)
        # button.pack(side=Tk.TOP, expand=0, fill=Tk.X)

    def make_dict_list(plot_dict, name):
        if not main_window[0]:
            main_window[0] = Tk.Tk()
            lwin = main_window[0]
        else:
            lwin = Tk.Toplevel()
        lwin.wm_title(name)

        scrollbar = Tk.Scrollbar(lwin, orient=Tk.VERTICAL)
        listbox = Tk.Listbox(lwin, width=40, yscrollcommand=scrollbar.set)
        scrollbar.config(command=listbox.yview)
        sorted_items = list(plot_dict.items())
        sorted_items.sort()
        for item in sorted_items:
            listbox.insert(Tk.END, item[0])
        listbox.pack(side=Tk.LEFT, expand=1, fill=Tk.BOTH)
        scrollbar.pack(side=Tk.LEFT, expand=0, fill=Tk.Y)

        def get_selected():
            items = listbox.curselection()
            if not items:
                return None, None
            return sorted_items[int(items[0])]

        def show(*args):
            key, plot = get_selected()
            if not plot:
                return
            if isinstance(plot, dict):
                make_dict_list(plot, key)
            else:
                show_plot(key, plot)

        listbox.bind("<Double-Button-1>", show)
        return lwin

    def make_dict_menu(d, root):
        menu = Tk.Menu(root)
        for name, val in list(d.items()):
            if isinstance(val, dict):
                m = make_dict_menu(val, menu)
                menu.add_cascade(label=name, menu=m)
            else:
                menu.add_command(label=name,
                                 command=lambda: show_plot(name, val))
        return menu

    lmenu = make_dict_list(iplot.plots, "IPlot")
    # lmenu.pack(side.Tk.TOP, expand=1, fill=Tk.BOTH)
    # menu = make_dict_menu(iplot.plots, win)
    # win.config(menu=menu)
    # menu.pack(side=Tk.TOP, expand=1, fill=Tk.BOTH)

    Tk.mainloop()


def shell_iplot():
    parser = OptionParser(usage, None, Option, version)
    parser.description = description
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
    parser.add_option("-p",
                      "--plot_type",
                      default="ps",
                      dest="plot_type",
                      help="ps or png")
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
    iplot = IPlot(
        options.origin_prefix,
        options.plot_type,
        options.plot_variables.split(","),
        output_prefix=options.dest_prefix,
    )
    if options.batch:
        iplot.print_all()
    else:
        make_menu(iplot)


if __name__ == "__main__":
    shell_iplot()
