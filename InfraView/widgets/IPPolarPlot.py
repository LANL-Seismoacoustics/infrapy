import numpy as np
import math

import pyqtgraph as pg


class IPPolarPlot(pg.PlotItem):

    def __init__(self):
        super().__init__()
        pg.setConfigOptions(antialias=True)
        self.drawPlot(250.0)

    def drawPlot(self, min_trace_vel):
        self.hideAxis('top')
        self.hideAxis('bottom')
        self.hideAxis('right')
        self.hideAxis('left')
        self.vb.setAspectLocked(lock=True, ratio=1)
        #self.addLine(x=0, pen=0.5, z=10)
        #self.addLine(y=0, pen=0.5, z=10)
        self.setRange(min_trace_vel*12/10)
        self.hideButtons()

    def drawRadials(self, min_trace_vel, count=16):
        # Draws the radial "spokes" in the plot
        # count is the number of lines to draw.  the spacing of lines will be 360/count degrees
        radius = 1/min_trace_vel
        angles = np.arange(0, 360, 360./count)

        for a in angles:
            x = radius * math.sin(math.radians(a))
            y = radius * math.cos(math.radians(a))
            self.plot([0,x],[0,y], pen=(0.75))

    def addLabels(self, min_trace_vel):
        N_text = pg.TextItem(html='<b>N</b>')
        self.addItem(N_text)
        N_text.setPos(0., 11.5 / (min_trace_vel*10))

    def setRange(self, min_trace_vel, circle_count=5):
        # min trace vel is in m/s
        slowness = 1.0/min_trace_vel
        self.setXRange(-slowness, slowness)
        self.setYRange(-slowness, slowness)

        # need to draw new circles.  First clear out old ones.
        for item in reversed(self.items):
            if type(item) is pg.QtWidgets.QGraphicsEllipseItem or type(item) is pg.TextItem:
                self.removeItem(item)

        for r in range(1, circle_count+1):
            circle = pg.QtWidgets.QGraphicsEllipseItem(-r / (min_trace_vel*circle_count), -r / (min_trace_vel*circle_count), r * 2 / (min_trace_vel*circle_count), r * 2 / (min_trace_vel*circle_count))
            circle.setPen(pg.mkPen(0.75))
            circle.setZValue(10)
            self.addItem(circle)

        # plot the "spokes"
        self.drawRadials(min_trace_vel)

        # if the range changes, the location of the N label needs to change too...
        self.addLabels(min_trace_vel)

