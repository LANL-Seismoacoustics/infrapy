import pyqtgraph as pg


from pyqtgraph.graphicsItems.InfiniteLine import InfLineLabel

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal, pyqtSlot, Qt
from PyQt5.QtGui import QPen


class IPEventLine(pg.InfiniteLine):

    sigEventLineMoving = pyqtSignal(float)
    sigEventLineMoved = pyqtSignal(float)

    def __init__(self, eventPosition, eventID=None, parent=None):
        super().__init__(angle=90, movable=False, pen=(22,43,72), label='')
        self.parent = parent
        self.setPos(eventPosition)
        self.setZValue(20)
        self.label = InfLineLabel(self, text=eventID, movable=False, color=(22,43,72))
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

    @pyqtSlot(float)
    def updatePosition(self, newPosition):
        self.setPos(newPosition)
            

class IPArrivalLine(pg.InfiniteLine):
    def __init__(self, position, arrival_type, parent=None):
        '''
        valid values for type are 'stratospheric', 'tropospheric', or 'thermospheric'
        '''
        super().__init__(angle=90, movable=False, pen=(50,140,75), label='')

        if arrival_type not in ['Stratospheric', 'Tropospheric', 'Thermospheric']:
            # complain and bail
            raise ValueError
            return

        self.parent = parent
        self.setPos(position)

        self.label = InfLineLabel(self, text=arrival_type, movable=False, color=(50,140,75))
        self.label.fill = pg.mkBrush(255, 255, 255, 0)
        if arrival_type == 'Stratospheric': 
            self.label.setPosition(0.85)
        elif arrival_type == 'Thermospheric':
            self.label.setPosition(0.8)
        elif arrival_type == 'Tropospheric':
            self.label.setPosition(0.9)
        self.label.setZValue(20)

    @pyqtSlot(float)
    def updatePosition(self, new_position):
        self.setPos(new_position)

    