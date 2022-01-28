import pyqtgraph as pg
from pyqtgraph import Point

from pyqtgraph.graphicsItems.InfiniteLine import InfLineLabel

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal, Qt, QPoint
from PyQt5.QtGui import QPen
from PyQt5.QtWidgets import QMenu

import pdb

class IPEventLine(pg.InfiniteLine):

    sigEventLineMoving = pyqtSignal(float)
    sigEventLineMoved  = pyqtSignal(float)

    def __init__(self, eventPosition, eventID=None, parent=None):
        super().__init__(angle=90, movable=True, pen='b', label='')
        self.parent = parent
        self.setPos(eventPosition)
        self.setZValue(20)
        self.label = InfLineLabel(self, text=eventID, movable=True, color=(0,0,0))
        self.label.fill = pg.mkBrush(255,255,255,0)
        self.label.setPosition(.95)
        self.label.setZValue(20)
        self.setHoverPen('b',width=3)

    def mouseClickEvent(self, evt):
        if evt.button() == Qt.RightButton:
            evt.accept()
            pos = evt.screenPos()
            # does nothing right now

    def mouseDragEvent(self, ev):

        if self.movable and ev.button() == QtCore.Qt.LeftButton:
            pos = self.pos().x()
            self.sigEventLineMoving.emit(pos)           

        super().mouseDragEvent(ev)          

    def setColor(self, color):
        self.setPen(QPen(color))
        self.setHoverPen(color, width=3)

    def setID(self, label):
        self.label.setText(label)

            

            