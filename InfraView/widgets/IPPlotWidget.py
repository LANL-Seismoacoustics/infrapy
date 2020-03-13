import numpy as np

import pyqtgraph as pg
from pyqtgraph import LinearRegionItem

from PyQt5 import QtCore, QtGui
from PyQt5.QtGui import QCursor
from PyQt5.QtWidgets import QMenu
from PyQt5.QtCore import pyqtSignal, QPoint


class NonScientific(pg.AxisItem):
    def __init__(self, *args, **kwargs):
        super(NonScientific, self).__init__(*args, **kwargs)

    def tickStrings(self, values, scale, spacing):
        # This line return the NonScientific notation value
        return ["%0.1f" % x for x in np.array(values).astype(float)]


class IPCustomViewBox(pg.ViewBox):

    def __init__(self, parent=None):
        super(IPCustomViewBox, self).__init__(parent)

        self.menu = None    # override pyqtgraph viewboxmenu
        self.menu = self.getMenu()

    def raiseContextMenu(self, ev):
        if not self.menuEnabled():
            return
        menu = self.getMenu()
        pos = ev.screenPos()
        menu.popup(QPoint(pos.x(), pos.y()))

    def getMenu(self):
        if self.menu is None:
            self.menu = QMenu()
            # self.viewAll = QtGui.QAction("View All", self.menu)
            self.exportImage = QtGui.QAction("Export Image", self.menu)
            # self.exportImage.triggered.connect(self.export)
            self.menu.addAction(self.exportImage)
        return self.menu


class IPPlotWidget(pg.PlotItem):

    sigNoiseRegionChanged = pyqtSignal(tuple)
    sigSignalRegionChanged = pyqtSignal(tuple)
    sigFreqRegionChanged = pyqtSignal(tuple)

    _noise_region = None
    _signal_region = None
    _freq_region = None

    _pickable = False

    def __init__(self, mode='plain', y_label_format=None, pickable=False):
        '''
        mode can currently be 'Waveform' or 'PSD' or 'Plain'

        label_format can be 'scientific' or nonscientific
        '''
        # if y_label_format == 'nonscientific':
        #     super().__init__(axisItems={'left': NonScientific(orientation='left')}, viewbox=IPCustomViewBox())
        # else:
        #     super().__init__(viewBox=IPCustomViewBox())

        if y_label_format == 'nonscientific':
            super().__init__(axisItems={'left': NonScientific(orientation='left')})
        else:
            super().__init__()

        # this will tell the widget if you can click on it and generate a 'pick'
        self._pickable = pickable

        self.showAxis('right')
        self.getAxis('right').setTicks('')
        self.showAxis('top')
        self.getAxis('top').setTicks('')

        self.getAxis('left').setWidth(80)

        if mode == 'waveform':
            self._noise_region = IPLinearRegionItem_Noise()
            self._signal_region = IPLinearRegionItem_Signal()

            self.addItem(self._noise_region)
            self.addItem(self._signal_region)

        elif mode == 'PSD':
            self._freq_region = IPFreqLinearRegionItem()

            self.addItem(self._freq_region)

        elif mode == 'plain':
            pass

    def setBackgroundColor(self, r, g, b):
        self.vb.setBackgroundColor(QtGui.QColor(r, g, b))

    def xaxis(self):
        return self.vb.XAxis

    def yaxis(self):
        return self.vb.YAxis

    def mouseClickEvent(self, evt):
        if evt.button() == QtCore.Qt.RightButton:
            if self._pickable:
                if evt.button() == QtCore.Qt.LeftButton:
                    p = QCursor.pos()   # This is the global coordinate of the mouse
                    scene_pos = evt.scenePos()

                    mousepoint = self.vb.mapSceneToView(scene_pos)
                    newx = mousepoint.x()
                    evt.accept()
        else:
            evt.ignore()

    def mouseDragEvent(self, evt):
        if evt.button() == QtCore.Qt.RightButton:
            evt.ignore()
        else:
            pass

    def getNoiseRegion(self):
        return self._noise_region

    def getSignalRegion(self):
        return self._signal_region

    def getFreqRegion(self):
        return self._freq_region

    def getNoiseRegionRange(self):
        if self._noise_region is not None:
            return self._noise_region.getRegion()
        else:
            return None

    def getSignalRegionRange(self):
        if self._signal_region is not None:
            return self._signal_region.getRegion()
        else:
            return None

    def getFreqRegionRange(self):
        if self._freq_region is not None:
            return self._freq_region.getRegion()
        else:
            return None

    def setNoiseRegionRange(self, range):
        if self._noise_region is not None:
            self._noise_region.setRegion(range)

    def setSignalRegionRange(self, range):
        if self._signal_region is not None:
            self._signal_region.setRegion(range)

    def setFreqRegionRange(self, range):
        if self._freq_region is not None:
            self._freq_region.setRegion(range)

    def copySignalRange(self, sourceSigRegion):
        reg = sourceSigRegion.getRegion()
        self._signal_region.setRegion(reg)

    def copyNoiseRange(self, sourceNoiseRegion):
        reg = sourceNoiseRegion.getRegion()
        self._noise_region.setRegion(reg)


class IPLinearRegionItem_Noise(LinearRegionItem):

    def __init__(self, values=[0, 1], orientation=LinearRegionItem.Vertical, brush=None, movable=True, bounds=None):
        super().__init__(values=values, orientation=orientation, brush=brush, movable=movable, bounds=bounds)
        self.setZValue(15)
        brush = QtGui.QBrush(QtGui.QColor(200, 100, 100, 50))
        self.setBrush(brush)

    def mouseClickEvent(self, ev):

        if ev.button() == QtCore.Qt.RightButton:
            ev.accept()
            pos = ev.screenPos()
            self.showMenu(pos)

    def showMenu(self, position):
        menu = QMenu()
        deleteSelection = menu.addAction("Remove Selection", self.hideMe)
        ret = menu.exec_(QPoint(position.x(), position.y()))

    def hideMe(self):
        self.setVisible(False)

    def showMe(self):
        self.setVisible(True)


class IPLinearRegionItem_Signal(LinearRegionItem):

    def __init__(self, values=[0, 1], orientation=LinearRegionItem.Vertical, brush=None, movable=True, bounds=None):
        super().__init__(values=values, orientation=orientation, brush=brush, movable=movable, bounds=bounds)
        self.setZValue(15)
        brush = QtGui.QBrush(QtGui.QColor(176, 224, 230, 100))
        self.setBrush(brush)

    def mouseClickEvent(self, ev):

        if ev.button() == QtCore.Qt.RightButton:
            ev.accept()
            pos = ev.screenPos()
            self.showMenu(pos)

    def showMenu(self, position):
        menu = QMenu()
        deleteSelection = menu.addAction("Remove Selection", self.hideMe)
        ret = menu.exec_(QPoint(position.x(), position.y()))

    def hideMe(self):
        self.setVisible(False)

    def showMe(self):
        self.setVisible(True)


class IPFreqLinearRegionItem(LinearRegionItem):

    def __init__(self, values=[np.log10(0.5), np.log10(5)], orientation=LinearRegionItem.Vertical, brush=None, movable=True, bounds=None):
        super().__init__(values=values, orientation=orientation, brush=brush, movable=movable, bounds=bounds)
        self.setZValue(15)
        brush = QtGui.QBrush(QtGui.QColor(100, 100, 100, 50))
        self.setBrush(brush)

    def mouseClickEvent(self, ev):

        if ev.button() == QtCore.Qt.RightButton:
            ev.accept()
            pos = ev.screenPos()
            self.showMenu(pos)

    def showMenu(self, position):
        menu = QMenu()
        # delete_selection = menu.addAction("Remove Selection", self.hideMe)
        ret = menu.exec_(QPoint(position.x(), position.y()))

    def hideMe(self):
        self.setVisible(False)

    def showMe(self):
        self.setVisible(True)
