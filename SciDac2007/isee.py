# Author: Antonio X. Garcia <antoniox.garcia@gmail.com>
# Modified by: Massimo Di Pierro <mdipierro@cs.depaul.edu>
# Copyright (c) 2007, DePaul University
# License: BSD Style

import os
import sys
import glob
import wx
from enthought.mayavi.app import Mayavi
from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
from enthought.mayavi.modules.outline import Outline
from enthought.mayavi.modules.text import Text
from enthought.mayavi.modules.surface import Surface
from enthought.mayavi.modules.orientation_axes import OrientationAxes

# usage: show.py FILE(s) [OPTION]

VERSION_INFO = (
    "Show.py v1.0\n"
    "Copyright (c) 2007, DePaul University\n"
    "All rights reserved.\n"
    "License: BSD Style\n"
    "Written by Antonio X. Garcia <antonioxgarcia@gmail.com>\n"
    "and Massimo Di Pierro <mdipierro@cs.depaul.edu>"
)


def parse():
    """Parses file input from command-line arguments."""
    if len(sys.argv) < 2:
        print("Usage: show.py FILE(s)")
        sys.exit(1)
    filepatterns = sys.argv[1]
    return [glob.glob(p) for p in filepatterns.split(",")]


class Watcher(wx.Timer):
    def __init__(self, interval, callback, *args, **kwargs):
        super().__init__()
        self.callback = callback
        self.args = args
        self.kwargs = kwargs
        self.Start(interval)

    def Notify(self):
        self.callback(*self.args, **self.kwargs)


class MayaViShow(Mayavi):
    def __init__(self, filelists, title):
        super().__init__()
        self.filelists = filelists
        self.title = title
        self.data = []
        self.watcher = None
        self.stop = True
        self.imagecount = 0
        self.sceneTextCount = None
        self.filestatus = []

    def run(self):
        win = self.script.engine.application.gui.GetTopWindow()
        win.CenterOnScreen()
        self.script.new_scene()

        for idx, filegroup in enumerate(self.filelists):
            f = filegroup[0]
            reader = VTKFileReader()
            reader.initialize(f)
            self.data.append(reader)
            self.filestatus.append(os.stat(f).st_mtime)
            self.script.add_source(reader)
            self.add_modules(idx, add_text=(idx == 0))

        self.add_animation_menu()

    def add_animation_menu(self):
        """Adds a custom 'Animations' menu to the GUI."""
        window = self.script.engine.application.gui.GetTopWindow()
        menubar = window.GetMenuBar()
        animation_menu = wx.Menu()
        pos = menubar.GetMenuCount() - 1
        menubar.Insert(pos, animation_menu, "Animations")

        items = {
            "Loop": self.AnimLoop,
            "Stop": self.AnimStop,
            "Watch": self.AnimWatch,
            # "Save as MPEG": self.AnimMPEG  # Placeholder
        }

        for label, handler in items.items():
            menu_item = wx.MenuItem(animation_menu, wx.ID_ANY, label)
            animation_menu.AppendItem(menu_item)
            window.Bind(wx.EVT_MENU, handler, menu_item)

    def AnimStop(self, evt):
        self.stop = True

    def AnimLoop(self, evt):
        self.stop = False
        self.timestep = 0
        self.watcher = Watcher(1000, self.AnimLoopRec)

    def AnimLoopRec(self):
        if self.stop:
            self.watcher = None
            return

        filelist = self.filelists[0]
        if self.timestep < len(filelist):
            for idx, _ in enumerate(self.filelists):
                self.data[idx].initialize(filelist[self.timestep])

            self.sceneTextCount.text = str(self.timestep)
            self.timestep += 1
        else:
            self.watcher = None

    def AnimWatch(self, evt):
        self.stop = False
        self.frame = 0
        self.watcher = Watcher(1000, self.AnimWatchRec)

    def AnimWatchRec(self):
        if self.stop:
            self.watcher = None
            return

        changed = False
        for idx, filegroup in enumerate(self.filelists):
            f = filegroup[0]
            if self.filestatus[idx] != os.stat(f).st_mtime:
                self.data[idx].reader.modified()
                self.data[idx].update()
                self.data[idx].data_changed = True
                changed = True

        if changed:
            self.sceneTextCount.text = str(self.frame)
            self.frame += 1

    def add_modules(self, idx, add_text=False):
        """Adds default modules for rendering."""
        self.script.add_module(Outline())
        self.script.add_module(OrientationAxes())

        surface = Surface()
        surface.enable_contours = True
        surface.actor.property.opacity = 0.5
        self.script.add_module(surface)

        # Title Text
        title_text = Text()
        title_text.text = self.title
        title_text.actor.scaled_text = False
        title_text.actor.text_property.font_size = 18
        self.script.add_module(title_text)

        width = title_text.actor.mapper.get_width(title_text.scene.renderer) / title_text.scene.renderer.size[0]
        height = title_text.actor.mapper.get_height(title_text.scene.renderer) / title_text.scene.renderer.size[1]
        title_text.x_position = 0.5 - width / 2
        title_text.y_position = 1 - height

        if add_text:
            self.sceneTextCount = Text()
            self.sceneTextCount.text = "0"
            self.sceneTextCount.actor.scaled_text = False
            self.sceneTextCount.actor.text_property.font_size = 24
            self.sceneTextCount.x_position = 0.95
            self.sceneTextCount.y_position = 0.05
            self.script.add_module(self.sceneTextCount)

    def snapshot(self):
        """Saves a PNG snapshot of the current scene."""
        if hasattr(self, 'outputdir') and self.outputdir:
            scene = self.script.engine.current_scene.scene
            path = os.path.join(self.outputdir, f"output{self.imagecount:05d}.png")
            scene.save_png(path)
            self.imagecount += 1


if __name__ == "__main__":
    files = parse()
    MayaViShow(files, "TEST").main()
