# XXX - No class, just generate all plots into a list

from optparse import OptionParser
import numpy
from scipy import stats
import tkinter as Tk
from matplotlib.axes import Subplot
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_ps import FigureCanvasPS
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigureCanvasTk

# from rpy import *
import re

USAGE = "python iplot.py\n"
VERSION = ("iplotv1.0"
           "\n  Copyright (c) 2007 Massimo Di Pierro"
           "\n  All rights reserved"
           "\n  License: GPL 3.0"
           "\n\n  Written by Massimo Di Pierro <mdipierro@cs.depaul.edu>")

DESCRIPTION = "plot the output of ibootstrap.py"

# r.library('Hmisc')


def clean(text):
    return re.sub("\W", "", text)


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
        for name, plot in list(self.plots.items()):
            self.print_plot(name, plot)

    def print_plot(self, name, plot):
        self.figure.clear()
        self.figure.add_subplot(plot)
        fname = self.output_prefix + name + self.extension
        self.canvas.print_figure(fname, 72)
        return fname

    def setup_plot(self, name, xlabel, ylabel):
        plot = Subplot(self.figure, 111)
        plot.set_title(name)
        plot.set_xlabel(xlabel)
        plot.set_ylabel(ylabel)
        self.plots[name] = plot
        return plot

    def make_error_bar(self, name, xdata, ydata, errors, xlabel, ylabel):
        plot = self.setup_plot(name, xlabel, ylabel)
        plot.errorbar(xdata, ydata, errors, fmt="ro", markerfacecolor=None)

    def make_hist(self,
                  name,
                  data,
                  length,
                  xlabel="x",
                  ylabel="y"):  # prob='T'
        # be sure to add a 'rug', whatever that is
        plot = self.setup_plot(name, xlabel, ylabel)
        plot.hist(data, length)
        # plot.set_xticks(data) # was not a good 'rug'

    def make_plot(self,
                  name,
                  xdata,
                  ydata,
                  xlabel="x",
                  ylabel="y"):  # type='p'
        plot = self.setup_plot(name, xlabel, ylabel)
        plot.plot(xdata, ydata, "ko", markerfacecolor=None)

    def plot_raw_data(self, filename):
        for line in open(filename, "r"):
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(x) for x in items[1:]]
            name = filename[:-4] + "_%s" % clean(tag)
            self.make_plot(
                name,
                xdata=list(range(len(data))),
                ydata=data,
                xlabel="step",
                ylabel=tag,
            )
            name = filename[:-4] + "_%s_hist" % clean(tag)
            self.make_hist(name,
                           data,
                           length=len(data) / 20,
                           xlabel=tag,
                           ylabel="frequency")
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
            self.make_plot(
                name,
                list(range(len(probs))),
                probs,
                xlabel="step",
                ylabel="probability " + tag,
            )

    def plot_trails(self, filename):
        for line in open(filename, "r"):
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(x) for x in items[1:]]
            name = filename[:-4] + "_%s" % clean(tag)
            self.make_plot(
                name,
                xdata=list(range(len(data))),
                ydata=data,
                xlabel="step",
                ylabel=tag,
            )

    def plot_samples(self, filename):
        for line in open(filename, "r"):
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(x) for x in items[1:]]
            name = filename[:-4] + "_%s_hist" % clean(tag)
            self.make_hist(name,
                           data,
                           length=len(data) / 10,
                           xlabel=tag,
                           ylabel="frequency")
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
        for legend in list(sets.keys()):
            x, y, yminus, yplus = sets[legend]
            # matplotlib takes offsets for errors, not absolute y positions
            error_low = numpy.subtract(y, yminus)
            error_high = numpy.subtract(yplus, y)
            name = filename[:-4] + "_%s" % clean(legend)
            self.make_error_bar(
                name,
                x,
                y,
                [error_low, error_high],
                xlabel=tags[index + 1],
                ylabel=tags[0],
            )
            # v=r.TRUE


def make_window(iplot):
    win = Tk.Tk()
    win.wm_title("IPlot")

    canvas = FigureCanvasTk(iplot.figure, master=win)
    canvas.draw()
    canvas.get_tk_widget().pack(side=Tk.LEFT, expand=1, fill=Tk.BOTH)

    frame = Tk.Frame(win)
    frame.pack(side=Tk.LEFT, expand=1, fill=Tk.BOTH)

    lframe = Tk.Frame(frame)
    lframe.pack(side=Tk.TOP, expand=1, fill=Tk.BOTH)

    scrollbar = Tk.Scrollbar(lframe, orient=Tk.VERTICAL)
    listbox = Tk.Listbox(lframe, width=40, yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    sorted_keys = list(iplot.plots.keys())
    sorted_keys.sort()
    for name in sorted_keys:
        listbox.insert(Tk.END, name)
    listbox.pack(side=Tk.LEFT, expand=1, fill=Tk.BOTH)
    scrollbar.pack(side=Tk.LEFT, expand=0, fill=Tk.Y)

    displayed = [None, None]

    def key_plot():
        items = listbox.curselection()
        if not items:
            return None, None
        iplot.figure.clear()
        key = sorted_keys[int(items[0])]
        plot = iplot.plots[key]
        return key, plot

    def show(*args):
        key, plot = key_plot()
        if not plot:
            return
        displayed[0], displayed[1] = key, plot
        iplot.figure.clear()
        iplot.figure.add_axes(plot)
        canvas.draw()

    def print_plot(*args):
        key, plot = displayed
        if not plot:
            return
        name = iplot.print_plot(key, plot)
        top = Tk.Toplevel(win)
        Tk.Label(top, text="Printed to file '%s'" % name).pack()
        Tk.Button(top, text="OK", command=lambda: top.destroy()).pack()
        win.wait_window(top)

    listbox.bind("<Double-Button-1>", show)
    button = Tk.Button(frame, text="Print", command=print_plot)
    button.pack(side=Tk.TOP, expand=0, fill=Tk.X)

    Tk.mainloop()


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
        make_window(iplot)


if __name__ == "__main__":
    shell_iplot()
