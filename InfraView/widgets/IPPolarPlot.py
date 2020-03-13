import numpy as np

import pyqtgraph as pg
from pyqtgraph import LinearRegionItem

from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import Qt, pyqtSignal


class IPPolarPlot(pg.PlotItem):

    def __init__(self):
        super().__init__()
        self.drawPlot()

    def drawPlot(self):
        self.hideAxis('top')
        self.hideAxis('bottom')
        self.hideAxis('right')
        self.hideAxis('left')
        self.vb.setAspectLocked(lock=True, ratio=1)
        self.addLine(x=0, pen=0.8, z=10)
        self.addLine(y=0, pen=0.8, z=10)
        self.setXRange(-.004, .004)
        self.setYRange(-.004, .004)
        self.hideButtons()
        for r in range(1, 10):
            circle = pg.QtGui.QGraphicsEllipseItem(-r / 2500., -r / 2500., r * 2 / 2500., r * 2 / 2500.)
            circle.setPen(pg.mkPen(0.8))
            circle.setZValue(10)
            self.addItem(circle)
        self.addLabels()

    def addLabels(self):
        N_text = pg.TextItem(html='<b>N</b>')
        self.addItem(N_text)
        N_text.setPos(0., 10.5 / 2500.)
