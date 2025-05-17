import csv
import re
import argparse
from rpy2.robjects import r
from rpy2.robjects.packages import importr

# Script information
USAGE = "python iplot.py\n"
VERSION = ("iplot v1.0"
           "\n  Copyright (c) 2007 Massimo Di Pierro"
           "\n  All rights reserved"
           "\n  License: GPL 3.0"
           "\n\n  Written by Massimo Di Pierro <mdipierro@cs.depaul.edu>")
DESCRIPTION = "Plot the output of ibootstrap.py"

# R library initialization
Hmisc = importr("Hmisc")

def clean(text):
    """Cleans up text by replacing spaces and slashes."""
    return re.sub(r"\s+", "", text.strip().replace("/", "_div_"))

class IPlot:
    def __init__(self, filename, plot_type, items=None):
        if items is None:
            items = []
        self.type = plot_type
        if self.type == "quartz":
            r.quartz()

        # Plot different data categories
        self.plot_raw_data(f"{filename}_raw_data.csv")
        self.plot_autocorrelations(f"{filename}_autocorrelations.csv")
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

    def plot_data(self, filename, plot_function):
        """Handles reading data from a CSV and plotting it using a specific plot function."""
        with open(filename, "r", newline='') as file:
            reader = csv.reader(file, delimiter=",", quoting=csv.QUOTE_NONNUMERIC)
            for items in reader:
                tag = items[0]
                data = items[1:]
                plot_function(tag, data, filename)

    def plot_raw_data(self, filename):
        """Plots raw data and its corresponding histogram and probability plot."""
        def plot_function(tag, data, filename):
            base = f"{filename[:-4]}_{clean(tag)}"

            # Raw data plot
            self.begin(base)
            r.plot(x=list(range(len(data))), y=data, xlab="step", ylab=tag, main="")
            self.end()

            # Plot histogram
            self.begin(f"{base}_hist")
            r.hist(data, n=max(1, len(data) // 20), xlab=tag, ylab="frequency", main="", prob=True)
            r.rug(data)
            self.end()

            # Plot probability data
            mu = r.mean(data)[0]
            sd = r.sd(data)[0]

            def pnorm_val(x):
                return float(r.pnorm((x - mu) / sd)[0])

            probs = [min(pnorm_val(x), 1 - pnorm_val(x)) for x in data]
            self.begin(f"{base}_probability")
            r.plot(list(range(len(probs))), probs, xlab="step", ylab=f"probability {tag}", main="", type="p")
            self.end()

        self.plot_data(filename, plot_function)

    def plot_autocorrelations(self, filename):
        """Plots autocorrelations."""
        def plot_function(tag, data, filename):
            self.begin(f"{filename[:-4]}_{clean(tag)}")
            r.plot(x=list(range(len(data))), y=data, xlab="step", ylab=tag, main="")
            self.end()

        self.plot_data(filename, plot_function)

    def plot_trails(self, filename):
        """Plots trails from the given CSV file."""
        def plot_function(tag, data, filename):
            self.begin(f"{filename[:-4]}_{clean(tag)}")
            r.plot(x=list(range(len(data))), y=data, xlab="step", ylab=tag, main="", type="p")
            self.end()

        self.plot_data(filename, plot_function)

    def plot_samples(self, filename):
        """Plots samples and histograms from the given CSV file."""
        def plot_function(tag, data, filename):
            self.begin(f"{filename[:-4]}_{clean(tag)}_hist")
            r.hist(data, n=max(1, len(data) // 10), xlab=tag, ylab="frequency", main="", prob=True)
            r.rug(data)
            self.end()

        self.plot_data(filename, plot_function)

    def plot_min_mean_max(self, filename, xlab=None):
        """Plot min/mean/max with error bars from the given CSV file."""
        if not xlab:
            xlab = ["t"]

        with open(filename, "r", newline='') as file:
            lines = list(csv.reader(file, delimiter=",", quoting=csv.QUOTE_NONNUMERIC))

        tags = lines[0]
        index = next((i for i, tag in enumerate(tags) if tag.startswith("[")), len(tags) - 3)

        sets = {}
        for items in lines[1:]:
            tag = items[0]
            data = items[1:]

            legend = " ".join(
                f"{tags[i]}={data[i - 1]}"
                for i in range(1, len(tags) - 3)
                if tags[i] not in xlab
            )

            if legend not in sets:
                sets[legend] = ([], [], [], [])

            x, y, yminus, yplus = sets[legend]
            t = data[index]
            x.append(t)
            y.append(data[-2])
            yminus.append(data[-3])
            yplus.append(data[-1])

        for legend, (x, y, yminus, yplus) in sets.items():
            self.begin(f"{filename[:-4]}_{clean(legend)}")
            r.errbar(x, y, yminus, yplus, xlab=tags[index + 1], ylab=tags[0], main="")
            self.end()

def shell_iplot():
    """Handles command-line options and invokes IPlot."""
    parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
    parser.add_argument("-o", "--origin_prefix", default="ibootstrap", help="Prefix for input filenames")
    parser.add_argument("-p", "--plot_type", default="ps", help="Plot type: 'ps', 'png', or 'quartz'")
    parser.add_argument("-v", "--plot_variables", default="", help="Comma-separated list of variables to plot")
    parser.add_argument("-f", "--fit", default=[], action="append", help="Fits to be performed (not implemented)")

    options = parser.parse_args()
    items = options.plot_variables.split(",") if options.plot_variables else []

    if options.fit:
        print("Note: -f option is not implemented yet.")

    IPlot(options.origin_prefix, options.plot_type, items)

if __name__ == "__main__":
    shell_iplot()

