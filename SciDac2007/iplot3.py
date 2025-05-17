import numpy
from scipy import stats
import tkinter as Tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_ps import FigureCanvasPS
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg as FigureCanvasTk,
    NavigationToolbar2Tk as NavigationToolbar2Tk,
)

import csv
from optparse import OptionParser
import re
from rpy import r

# Script information
USAGE = "python iplot.py\n"
VERSION = ("iplot v1.0"
           "\n  Copyright (c) 2007 Massimo Di Pierro"
           "\n  All rights reserved"
           "\n  License: GPL 3.0"
           "\n\n  Written by Massimo Di Pierro <mdipierro@cs.depaul.edu>")
DESCRIPTION = "Plot the output of ibootstrap.py"


def clean(text):
    """Cleans up text by replacing spaces with empty and slashes with '_div_'."""
    return re.sub(r"\s+", "", text.replace("/", "_div_"))


def gen_plot(figure, plot_args):
    """Generates a plot using the provided figure and arguments."""
    figure.clear()
    plot_func = plot_args[0]
    plot = plot_func(figure, *plot_args[1])
    figure.add_subplot(plot, 111)
    return plot


# Notice that I defer the actual generation of the plots, and just pass
# along the function used to make them, along with some arguments. By convention
# it is assumed the function takes a matplotlib figure, then those arguments.
# Using a lambda would have been nicer, however due to the semantics of python
# loops, specifically writing over the loop variables, by the time the lambda
# was called many of its closure variables were overwritten


class IPlot:
    def __init__(self, filename, plot_type, items=None, output_prefix=""):
        if items is None:
            items = []
        self.type = plot_type
        self.output_prefix = output_prefix
        self.figure = Figure((8.5, 11))
        
        # Choose the canvas based on plot type
        if plot_type == "png":
            self.canvas = FigureCanvasAgg(self.figure)
            self.extension = ".png"
        else:
            self.canvas = FigureCanvasPS(self.figure)
            self.extension = ".ps"
        
        self.plots = {}

        self.plot_raw_data(filename + "_raw_data.csv")
        self.plot_autocorrelations(filename + "_autocorrelations.csv")
        self.plot_trails(filename + "_trails.csv")
        self.plot_samples(filename + "_samples.csv")
        self.plot_min_mean_max(filename + "_min_mean_max.csv", items)
        self.figure.clear()

    def print_all(self):
        """Prints all plots."""
        def print_dict(d):
            for name, val in list(d.items()):
                if isinstance(val, dict):
                    print_dict(val)
                else:
                    self.print_plot(name, val)

        print_dict(self.plots)

    def print_plot(self, name, plot_args):
        """Generates and saves a plot to a file."""
        plot = gen_plot(self.figure, plot_args)
        fname = self.output_prefix + name + self.extension
        self.canvas.print_figure(fname, 72)
        return fname

    def setup_plot(self, figure, name, xlabel, ylabel):
        """Sets up a basic plot with title and axis labels."""
        plot = figure.add_subplot(111)
        plot.set_title(name)
        plot.set_xlabel(xlabel)
        plot.set_ylabel(ylabel)
        return plot

    def make_error_bar(self, figure, name, xdata, ydata, errors, xlabel, ylabel):
        """Generates an error bar plot."""
        plot = self.setup_plot(figure, name, xlabel, ylabel)
        plot.errorbar(xdata, ydata, errors, fmt="ko", markerfacecolor=None)
        return plot

    def make_hist(self, figure, name, data, length, xlabel="x", ylabel="y"):
        """Generates a histogram plot."""
        plot = self.setup_plot(figure, name, xlabel, ylabel)
        plot.hist(data, length, facecolor="w", edgecolor="k")
        return plot

    def make_plot(self, figure, name, xdata, ydata, xlabel="x", ylabel="y"):
        """Generates a simple line plot."""
        plot = self.setup_plot(figure, name, xlabel, ylabel)
        plot.plot(xdata, ydata, "k-", markerfacecolor=None)
        return plot

    def plot_raw_data(self, filename):
        """Plots raw data from the given CSV file."""
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
            plots[name] = self.make_hist, (name, data, 10, tag, "frequency")
            
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

    def plot_autocorrelations(self, filename):
        """Plots autocorrelations from the given CSV file."""
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

    def plot_trails(self, filename):
        """Plots trails from the given CSV file."""
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
        """Plots samples from the given CSV file."""
        quants = {}
        self.plots[filename[:-4]] = quants
        for line in open(filename, "r"):
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(x) for x in items[1:]]
            name = filename[:-4] + "_%s_hist" % clean(tag)
            quants[name] = self.make_hist, (name, data, 10, tag, "frequency")

    def plot_min_mean_max(self, filename, xlab=None):
        """Plots min, mean, and max from the given CSV file."""
        lines = list(
            csv.reader(open(filename, "r"), delimiter=",", quoting=csv.QUOTE_NONNUMERIC)
        )
        tags = lines[0]
        if not xlab or xlab[0] == "":
            xlab = [tags[1]]
        index = next(
            (k for k, tag in enumerate(tags) if tag == xlab[0]), None
        )
        if index is None:
            raise ValueError(f"Error: {xlab} not found in the tags.")
        sets = {}
        for items in lines[1:]:
            tag = items[0]
            data = items[1:]
            legend = "".join(f"{tags[i]}={data[i-1]} " for i in range(1, len(tags)) if tags[i] not in xlab)
            if legend not in sets:
                sets[legend] = ([], [], [], [])
            x, y, yminus, yplus = sets[legend]
            t = data[index]
            x.append(t)
            y.append(data[-2])
            yminus.append(data[-3])
            yplus.append(data[-1])

        for legend, (x, y, yminus, yplus) in sets.items():
            error_low = numpy.subtract(y, yminus)
            error_high = numpy.subtract(yplus, y)
            name = filename[:-4] + "_%s" % clean(legend)
            self.plots[name] = self.make_error_bar, (
                name,
                x,
                y,
                [error_low, error_high],
                tags[index + 1],
                tags[0],
            )


def make_menu(iplot):
    """Creates the menu interface for interacting with the plot."""
    main_window = Tk.Tk()

    def show_plot(name, plot_args):
        pwin = Tk.Toplevel()
        pwin.wm_title(name)
        figure = Figure(facecolor="#ffffff")
        plot = gen_plot(figure, plot_args)

        frame_mpl = Tk.Frame(master=pwin)
        frame_mpl.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

        canvas = FigureCanvasTk(figure, master=frame_mpl)
        canvas.draw()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, frame_mpl)
        toolbar.update()
        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        frame_edit = Tk.Frame(master=pwin)
        frame_edit.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)

        def make_entry(name, default, row):
            Tk.Label(master=frame_edit, text=name).grid(row=row, column=1)
            entry = Tk.Entry(master=frame_edit)
            entry.insert(0, default)
            entry.grid(row=row, column=2)
            return entry

        def set_attributes(*args):
            plot.set_xlabel(xlabel_entry.get())
            plot.set_ylabel(ylabel_entry.get())
            plot.set_title(title_entry.get())
            plot.set_axis_bgcolor(background_entry.get())
            figure.set_facecolor(background_entry.get())
            canvas.draw()

        xlabel_entry = make_entry("X-Label", plot.xaxis.label.get_text(), 1)
        ylabel_entry = make_entry("Y-Label", plot.yaxis.label.get_text(), 2)
        title_entry = make_entry("Title", plot.title.get_text(), 3)
        background_entry = make_entry("Background", plot.get_axis_bgcolor(), 4)
        button = Tk.Button(frame_edit, text="Set", command=set_attributes)
        button.grid(row=5, column=1, columnspan=2)

        # 5) Optional. For each plot it should allow to edit plot parameters such as labels, title, colors etc. by default background should be white.

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

    # main_window[0] = Tk.Tk()
    # menu = pmw.MenuBar(main_window[0])
    # menu.pack(fill='y')
    # lmenu = make_dict_menu(iplot.plots, menu)
    lmenu = make_dict_list(iplot.plots, "IPlot")
    # lmenu.pack(side.Tk.TOP, expand=1, fill=Tk.BOTH)
    # menu = make_dict_menu(iplot.plots, win)
    # win.config(menu=menu)
    # menu.pack(side=Tk.TOP, expand=1, fill=Tk.BOTH)

    Tk.mainloop()


def shell_iplot():
    """Main function to handle command-line arguments and execute plotting."""
    parser = OptionParser(USAGE, None, OptionParser, VERSION)
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
        make_menu(iplot)


if __name__ == "__main__":
    shell_iplot()
