import csv
from optparse import OptionParser
import re
from rpy2.robjects import r

# Script information
USAGE = "python iplot.py\n"
VERSION = ("iplot v1.0"
           "\n  Copyright (c) 2007 Massimo Di Pierro"
           "\n  All rights reserved"
           "\n  License: GPL 3.0"
           "\n\n  Written by Massimo Di Pierro <mdipierro@cs.depaul.edu>")
DESCRIPTION = "Plot the output of ibootstrap.py"

# R library initialization
r.library("Hmisc")


def clean(text):
    """Clean a string by removing non-alphanumeric characters."""
    return re.sub("\W", "", text)


class IPlot:
    def __init__(self, filename, plot_type, items=None):
        if items is None:
            items = []
        self.type = plot_type
        if self.type == "quartz":
            r.quartz()

        # Plot different data categories
        self.plot_raw_data(f"{filename}_raw_data.csv")
        self.plot_trails(f"{filename}_trails.csv")
        self.plot_samples(f"{filename}_samples.csv")
        self.plot_min_mean_max(f"{filename}_min_mean_max.csv", items)

    def begin(self, filename):
        """Starts the plotting device (ps or png)."""
        if self.type == "ps":
            r.postscript(f"{filename}.ps")
        elif self.type == "png":
            r.png(f"{filename}.png")

    def end(self):
        """Ends the plotting device."""
        if self.type == "ps":
            r.dev_off()
        else:
            input("Press Enter to continue...")

    def plot_raw_data(self, filename):
        """Plot raw data and its histogram."""
        for line in open(filename, "r"):
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(x) for x in items[1:]]

            # Plot raw data
            self.begin(f"{filename[:-4]}_{clean(tag)}")
            r.plot(x=list(range(len(data))), y=data, xlab="step", ylab=tag)
            self.end()

            # Plot histogram
            self.begin(f"{filename[:-4]}_{clean(tag)}_hist")
            r.hist(data, n=len(data) // 20, xlab=tag, ylab="frequency", main="", prob="T")
            r.rug(data)
            self.end()

            # Plot probability data
            mu = r.mean(data)
            sd = r.sd(data)
            probs = [min(x, 1 - x) for x in [r.pnorm((x - mu) / sd) for x in data]]
            self.begin(f"{filename[:-4]}_{clean(tag)}_probability")
            r.plot(list(range(len(probs))), probs, xlab="step", ylab=f"probability {tag}", type="p")
            self.end()

    def plot_trails(self, filename):
        """Plot data trails."""
        for line in open(filename, "r"):
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(x) for x in items[1:]]

            # Plot trails data
            self.begin(f"{filename[:-4]}_{clean(tag)}")
            r.plot(x=list(range(len(data))), y=data, xlab="step", ylab=tag, type="p")
            self.end()

    def plot_samples(self, filename):
        """Plot sample data and its histogram."""
        for line in open(filename, "r"):
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(x) for x in items[1:]]

            # Plot histogram for samples
            self.begin(f"{filename[:-4]}_{clean(tag)}_hist")
            r.hist(data, n=len(data) / 10, xlab=tag, ylab="frequency", prob="T")
            r.rug(data)
            self.end()

    def plot_min_mean_max(self, filename, xlab=["t"]):
        """Plot min, mean, max data."""
        lines = open(filename, "r").readlines()
        tags = [x[1:-1] for x in re.compile('"[^"]*"').findall(lines[0])]

        if not xlab or xlab[0] == "":
            xlab = [tags[1]]

        index = next((i - 1 for i in range(1, len(tags) - 3) if tags[i] == xlab[0]), -1)
        if index < 0:
            print("Error: Unable to find xlab tag:", xlab)
            raise Exception

        sets = {}
        for line in lines[1:]:
            items = line.split(",")
            tag = items[0][1:-1]
            data = [float(p) for p in items[1:]]
            legend = "".join([f"{tags[i]}={data[i - 1]} " for i in range(1, len(tags) - 3) if tags[i] not in xlab])

            if legend not in sets:
                sets[legend] = ([], [], [], [])  # x, y, yminus, yplus

            x, y, yminus, yplus = sets[legend]
            t = data[index]
            x.append(t)
            y.append(data[-2])
            yminus.append(data[-3])
            yplus.append(data[-1])

        for legend, (x, y, yminus, yplus) in sets.items():
            self.begin(f"{filename[:-4]}_{clean(legend)}")
            r.errbar(x, y, yminus, yplus, xlab=tags[index + 1], ylab=tags[0])
            self.end()

def shell_iplot():
    """Main function to handle command-line arguments and execute plotting."""
    parser = OptionParser(USAGE, None, OptionParser, VERSION)
    parser.description = DESCRIPTION
    parser.add_option("-o", "--origin_prefix", default="ibootstrap", dest="origin_prefix", help="Prefix for input filenames")
    parser.add_option("-p", "--plot_type", default="ps", dest="plot_type", help="Plot type: ps or quartz")
    parser.add_option("-v", "--plot_variables", default="", dest="plot_variables", help="Comma-separated list of variables to plot")

    options, _ = parser.parse_args()
    plot = IPlot(options.origin_prefix, options.plot_type, options.plot_variables.split(","))


if __name__ == "__main__":
    shell_iplot()
