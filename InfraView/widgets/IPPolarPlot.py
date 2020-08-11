import numpy as np

import pyqtgraph as pg
from pyqtgraph import LinearRegionItem

import PyQt5
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import Qt, pyqtSignal


class IPPolarPlot(pg.PlotItem):

    def __init__(self):
        super().__init__()
        self.drawPlot(250.0)

    def drawPlot(self, min_trace_vel):
        self.hideAxis('top')
        self.hideAxis('bottom')
        self.hideAxis('right')
        self.hideAxis('left')
        self.vb.setAspectLocked(lock=True, ratio=1)
        self.addLine(x=0, pen=0.8, z=10)
        self.addLine(y=0, pen=0.8, z=10)
        self.setRange(min_trace_vel)
        self.hideButtons()
        self.addLabels(min_trace_vel)

    def addLabels(self, min_trace_vel):
        N_text = pg.TextItem(html='<b>N</b>')
        self.addItem(N_text)
        N_text.setPos(0., 10.5 / (min_trace_vel*10))

    def setRange(self, min_trace_vel):
        # min trace vel is in m/s
        slowness = 1/min_trace_vel
        self.setXRange(-slowness, slowness)
        self.setYRange(-slowness, slowness)

        # need to draw new circles.  First clear out old ones.
        for item in self.items:
            if type(item) is PyQt5.QtWidgets.QGraphicsEllipseItem or type(item) is pg.graphicsItems.TextItem.TextItem:
                self.removeItem(item)

        for r in range(1, 10):
            circle = pg.QtGui.QGraphicsEllipseItem(-r / (min_trace_vel*10), -r / (min_trace_vel*10), r * 2 / (min_trace_vel*10), r * 2 / (min_trace_vel*10))
            circle.setPen(pg.mkPen(0.8))
            circle.setZValue(10)
            self.addItem(circle)

        # if the range changes, the location of the N label needs to change too...
        self.addLabels(min_trace_vel)

