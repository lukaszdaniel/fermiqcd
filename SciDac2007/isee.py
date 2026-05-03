# Author: Antonio X. Garcia <antoniox.garcia@gmail.com>
# Modified by: Massimo Di Pierro <mdipierro@cs.depaul.edu>
# Copyright (c) 2007, DePaul University
# License: BSD Style

import os
import sys
import glob
import time
import threading
import pyvista as pv

# ================================================================================
# Show.py — Simple interactive VTK viewer based on PyVista
# ================================================================================
# 
# DESCRIPTION
# -----------
# This script provides a lightweight interactive viewer for VTK-compatible files
# (using PyVista). It supports both:
# 
#   • Animation over a sequence of files (e.g. simulation timesteps)
#   • Live monitoring of a file that is being updated during a simulation
# 
# The script is intended mainly for quick visualization of simulation outputs,
# especially those produced in time series (e.g. field snapshots).
# 
# INPUT
# -----
# The script takes a single command-line argument: a list of file patterns.
# 
# Example:
#     python show.py "*.vtk"
#     python show.py "out_*.vtk,data_*.vtk"
# 
# Patterns can be comma-separated and are expanded using glob.
# 
# FEATURES
# --------
# • Interactive 3D visualization using PyVista
# • Automatic loading of VTK (and other PyVista-supported formats)
# • Two animation modes:
#   
#   1. LOOP MODE ("l")
#      - Iterates over all input files in order
#      - Loops indefinitely
#      - Useful for visualizing time evolution
# 
#   2. WATCH MODE ("w")
#      - Monitors the first file in the list
#      - Reloads it when it changes on disk
#      - Useful for real-time simulation output
# 
# • Snapshot support (saving screenshots)
# • Keyboard-controlled interaction
# • Multithreaded animation (non-blocking UI)
# 
# KEYBOARD CONTROLS
# -----------------
# l : Start looping over all files (animation mode)
# w : Watch first file for changes (live mode)
# s : Stop animation (loop/watch)
# p : Save screenshot (if output directory is configured)
# 
# VISUALIZATION DETAILS
# ---------------------
# • Displays mesh with:
#     - semi-transparent surface (opacity = 0.5)
#     - optional edges disabled
# • Shows:
#     - axes
#     - bounding box
#     - filename index (during loop mode)
# 
# INTERNAL DESIGN
# ---------------
# • Uses PyVista Plotter for rendering
# • Background thread handles animation (loop/watch)
# • Mesh is reloaded completely on each update
# • File change detection based on modification time (mtime)
# 
# LIMITATIONS
# -----------
# • No CLI options beyond file patterns
# • WATCH mode only monitors the first file
# • No direct control over colormaps or scalar fields
# • Reloading clears and rebuilds the scene (not incremental)
# • Output directory for snapshots must be set manually in code
#
#
# ================================================================================


VERSION_INFO = (
    "Show.py v2.0\n"
    "Copyright (c) 2007, DePaul University\n"
    "All rights reserved.\n"
    "License: BSD Style\n"
    "Written by Antonio X. Garcia <antonioxgarcia@gmail.com>\n"
    "and Massimo Di Pierro <mdipierro@cs.depaul.edu>"
)


def parse():
    """Parse file patterns from CLI."""
    if len(sys.argv) < 2:
        print("Usage: show.py FILE(s)")
        sys.exit(1)

    filepatterns = sys.argv[1]
    files = []
    for pattern in filepatterns.split(","):
        files.extend(sorted(glob.glob(pattern)))

    return files


class PyVistaShow:
    def __init__(self, files, title="PyVista Viewer"):
        self.files = files
        self.title = title

        self.plotter = pv.Plotter()
        self.mesh_actor = None

        self.stop = True
        self.mode = None  # "loop" or "watch"
        self.index = 0
        self.file_mtime = None

        self.outputdir = None
        self.imagecount = 0

    # --------------------------------------------------------
    # Visualization setup
    # --------------------------------------------------------
    def setup_scene(self):
        self.plotter.add_text(self.title, position="upper_edge", font_size=12)
        self.plotter.add_axes()

        if not self.files:
            return

        mesh = pv.read(self.files[0])
        self.mesh_actor = self.plotter.add_mesh(
            mesh,
            opacity=0.5,
            show_edges=False
        )

        self.plotter.add_bounding_box()

    def update_mesh(self, filename):
        mesh = pv.read(filename)
        self.plotter.clear_actors()
        self.setup_scene()
        self.mesh_actor = self.plotter.add_mesh(mesh, opacity=0.5)

    # --------------------------------------------------------
    # Animation: LOOP
    # --------------------------------------------------------
    def loop(self):
        self.stop = False
        self.mode = "loop"

        def run():
            while not self.stop:
                if self.index >= len(self.files):
                    self.index = 0

                self.update_mesh(self.files[self.index])
                self.plotter.add_text(f"{self.index}", position="lower_right", font_size=10)

                self.index += 1
                time.sleep(1)

        threading.Thread(target=run, daemon=True).start()

    # --------------------------------------------------------
    # Animation: WATCH
    # --------------------------------------------------------
    def watch(self):
        self.stop = False
        self.mode = "watch"

        if not self.files:
            return

        self.file_mtime = os.stat(self.files[0]).st_mtime

        def run():
            while not self.stop:
                mtime = os.stat(self.files[0]).st_mtime
                if mtime != self.file_mtime:
                    self.file_mtime = mtime
                    self.update_mesh(self.files[0])
                time.sleep(1)

        threading.Thread(target=run, daemon=True).start()

    # --------------------------------------------------------
    # Stop animation
    # --------------------------------------------------------
    def stop_anim(self):
        self.stop = True

    # --------------------------------------------------------
    # Snapshot
    # --------------------------------------------------------
    def snapshot(self):
        if self.outputdir:
            path = os.path.join(self.outputdir, f"output{self.imagecount:05d}.png")
            self.plotter.screenshot(path)
            self.imagecount += 1

    def bind_keys(self):
        self.plotter.add_key_event("l", self.loop)
        self.plotter.add_key_event("s", self.stop_anim)
        self.plotter.add_key_event("w", self.watch)
        self.plotter.add_key_event("p", self.snapshot)

    def run(self):
        self.setup_scene()
        self.bind_keys()
        self.plotter.show(title=self.title)


if __name__ == "__main__":
    files = parse()
    app = PyVistaShow(files, "VTK Viewer")
    app.run()
