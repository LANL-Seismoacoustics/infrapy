import pyqtgraph as pg
from pyqtgraph.graphicsItems.GraphicsWidget import GraphicsWidget
from pyqtgraph.graphicsItems.LabelItem import LabelItem
from PyQt5 import QtCore, QtGui
from pyqtgraph import functions as fn
from pyqtgraph.Point import Point
from pyqtgraph.graphicsItems.ScatterPlotItem import drawSymbol
from pyqtgraph.graphicsItems.GraphicsWidgetAnchor import GraphicsWidgetAnchor


class IPSimpleLegend(GraphicsWidget, GraphicsWidgetAnchor):

    def __init__(self, size=None, offset=None):

        GraphicsWidget.__init__(self)
        GraphicsWidgetAnchor.__init__(self)
        self.setFlag(self.ItemIgnoresTransformations)
        self.layout = QtGui.QGraphicsGridLayout()
        self.setLayout(self.layout)
        self.size = size
        self.offset = offset
        if size is not None:
            self.setGeometry(QtCore.QRectF(0, 0, self.size[0], self.size[1]))

    def setParentItem(self, p):
        ret = GraphicsWidget.setParentItem(self, p)
        if self.offset is not None:
            offset = Point(self.offset)
            anchorx = 1 if offset[0] <= 0 else 0
            anchory = 1 if offset[1] <= 0 else 0
            anchor = (anchorx, anchory)
            self.anchor(itemPos=anchor, parentPos=anchor, offset=offset)
        return ret

    def addItem(self, name):
        """
        Add a new entry to the legend.

        ==============  ========================================================
        **Arguments:**
        item            A PlotDataItem from which the line and point style
                        of the item will be determined or an instance of
                        ItemSample (or a subclass), allowing the item display
                        to be customized.
        title           The title to display for this item. Simple HTML allowed.
        ==============  ========================================================
        """
        self.label = LabelItem(name)
        self.layout.addItem(self.label, 0, 0)
        self.updateSize()

    def removeItem(self, item):
        self.removeItem(self.label)                       # redraq box

    def updateSize(self):
        if self.size is not None:
            return

        height = self.label.height()
        width = self.label.width()

        self.setGeometry(0, 0, width, height)

    def boundingRect(self):
        return QtCore.QRectF(0, 0, self.width(), self.height())

    def paint(self, p, *args):
        p.setPen(fn.mkPen(255, 255, 255, 100))
        p.setBrush(fn.mkBrush(100, 100, 100, 80))
        # p.drawRect(self.boundingRect())

    def hoverEvent(self, ev):
        ev.acceptDrags(QtCore.Qt.LeftButton)

    def mouseDragEvent(self, ev):
        if ev.button() == QtCore.Qt.LeftButton:
            dpos = ev.pos() - ev.lastPos()
            self.autoAnchor(self.pos() + dpos)
