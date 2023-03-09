import numpy as np

import pyqtgraph as pg
from pyqtgraph import LinearRegionItem

from PyQt5.QtGui import QCursor, QColor, QBrush, QFont
from PyQt5.QtWidgets import QMenu, QAction
from PyQt5.QtCore import pyqtSignal, QPoint, Qt

from obspy.core import UTCDateTime


class NonScientific(pg.AxisItem):
    #def __init__(self, *args, **kwargs):
    #    super(NonScientific, self).__init__(*args, **kwargs)

    def tickStrings(self, values, scale, spacing):
        # This line return the NonScientific notation value
        return ["%0.1f" % x for x in np.array(values).astype(float)]

class IPWaveformTimeAxis(pg.AxisItem):
    # subclass the basic axis item, mainly to make custom time axis
    def __init__(self, est, *args, **kwargs):
        super().__init__(orientation='bottom', *args, **kwargs)
        # est is the "earliest_start_time"
        self.set_earliest_start_time(est)

        # make font size smaller
        # font = QFont()
        # font.setPointSize(12)
        # self.setTickFont(font)

    def tickStrings(self, values, scale, spacing):
        return [(self.earliest_start_time + value).strftime("%H:%M:%S") for value in values]

    def set_earliest_start_time(self, est):
        self.earliest_start_time = est

    def get_start_time(self):
        return self.earliest_start_time

class IPSpectrogramTimeAxis(pg.AxisItem):
    # subclass the basic axis item, mainly to make custom time axis
    def __init__(self, *args, **kwargs):
        super().__init__(orientation='bottom', *args, **kwargs)
        # st: (UTCDateTime) is the start_time of the window which will be the earliest start 
        # time of the waveforms plus the offset seconds of the signal/noise window
        self.set_start_time(UTCDateTime(0))

    def tickStrings(self, values, scale, spacing):
        return [(self.start_time + value).strftime("%H:%M:%S") for value in values]

    def set_start_time(self, st):
        self.start_time = st


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
            # self.viewAll = QAction("View All", self.menu)
            self.exportImage = QAction("Export Image", self.menu)
            # self.exportImage.triggered.connect(self.export)
            self.menu.addAction(self.exportImage)
        return self.menu


class IPPlotItem(pg.PlotItem):

    sigNoiseRegionChanged = pyqtSignal(tuple)
    sigSignalRegionChanged = pyqtSignal(tuple)
    sigFreqRegionChanged = pyqtSignal(tuple)

    noise_region = None
    signal_region = None
    freq_region = None

    pickable = False

    labi = None

    def __init__(self, mode='plain', y_label_format=None, pickable=False, est=None):
        '''
        mode: (str) can currently be 'waveform' or 'PSD' or 'Plain' or 'Spectrogram'
        est: (UTCDateTime) Earliest Start Time
        pickable: (bool) Can you click on plot to make a pick
        '''

        if y_label_format == 'nonscientific':
            super().__init__(axisItems={'left': NonScientific(orientation='left')})
        else:
            if mode == 'waveform':
                if est is None:
                    est = UTCDateTime(0)
                super().__init__(axisItems={'bottom': IPWaveformTimeAxis(est=est)})
            elif mode == 'spectrogram':
                if est is None:
                    est = UTCDateTime(0)
                super().__init__(axisItems={'bottom': IPSpectrogramTimeAxis(est=est)})
            else:
                super().__init__()

        # this will tell the widget if you can click on it and generate a 'pick'
        self.pickable = pickable

        self.showAxis('right')
        self.getAxis('right').setTicks('')
        self.showAxis('top')
        self.getAxis('top').setTicks('')

        self.getAxis('left').setWidth(80)
        # font = QFont()
        # font.setPointSize(10)
        # self.getAxis('left').setTickFont(font)
        # self.getAxis('bottom').setTickFont(font)

        if mode == 'waveform':
            self.noise_region = IPLinearRegionItem_Noise()
            self.signal_region = IPLinearRegionItem_Signal()

            self.addItem(self.noise_region)
            self.addItem(self.signal_region)

        elif mode == 'PSD':
            self.freq_region = IPFreqLinearRegionItem()

            self.addItem(self.freq_region)

        elif mode == 'plain':
            pass

    def setBackgroundColor(self, r, g, b):
        self.vb.setBackgroundColor(QColor(r, g, b))

    def setEarliestStartTime(self, est):
        self.getAxis('bottom').set_earliest_start_time(est)

    def get_start_time(self):
        return self.getAxis('bottom').get_start_time()

    def setPlotLabel(self, text):
        if self.labi is not None:
            self.vb.removeItem(self.labi)
        self.labi = pg.LabelItem(text=text)
        self.labi.setParentItem(self.vb)
        self.labi.anchor(itemPos=(0.0,0.0), parentPos=(0.0,0.0))

    def clearPlotLabel(self):
        if self.labi is not None:
            self.vb.removeItem(self.labi)
        self.labi = None

    def xaxis(self):
        return self.vb.XAxis

    def yaxis(self):
        return self.vb.YAxis

    def mouseClickEvent(self, evt):
        if evt.button() == Qt.RightButton:
            if self.pickable:
                if evt.button() == Qt.LeftButton:
                    p = QCursor.pos()   # This is the global coordinate of the mouse
                    scene_pos = evt.scenePos()

                    mousepoint = self.vb.mapSceneToView(scene_pos)
                    newx = mousepoint.x()
                    evt.accept()
        else:
            evt.accept()

    def mouseDragEvent(self, evt):
        if evt.button() == Qt.RightButton:
            evt.ignore()
        else:
            pass

    def getNoiseRegion(self):
        return self.noise_region

    def getSignalRegion(self):
        return self.signal_region

    def getFreqRegion(self):
        return self.freq_region

    def getNoiseRegionRange(self):
        if self.noise_region is not None:
            return self.noise_region.getRegion()
        else:
            return None

    def getSignalRegionRange(self):
        if self.signal_region is not None:
            return self.signal_region.getRegion()
        else:
            return None

    def getFreqRegionRange(self):
        if self.freq_region is not None:
            return self.freq_region.getRegion()
        else:
            return None

    def setNoiseRegionRange(self, range):
        if self.noise_region is not None:
            self.noise_region.setRegion(range)

    def setSignalRegionRange(self, range):
        if self.signal_region is not None:
            self.signal_region.setRegion(range)

    def setFreqRegionRange(self, range):
        if self.freq_region is not None:
            self.freq_region.setRegion(range)

    def copySignalRange(self, sourceSigRegion):
        reg = sourceSigRegion.getRegion()
        self.signal_region.setRegion(reg)

    def copyNoiseRange(self, sourceNoiseRegion):
        reg = sourceNoiseRegion.getRegion()
        self.noise_region.setRegion(reg)


class IPLinearRegionItem_Noise(LinearRegionItem):

    sig_IPRegion_Change_finished = pyqtSignal(tuple)

    def __init__(self, values=[0, 1], orientation=LinearRegionItem.Vertical, brush=None, movable=True, bounds=None, swapMode='block'):
        super().__init__(values=values, orientation=orientation, brush=brush, movable=movable, bounds=bounds, swapMode=swapMode)
        self.setZValue(15)
        brush = QBrush(QColor(255, 71, 71, 50))
        self.setBrush(brush)

    def mouseClickEvent(self, ev):
        if ev.button() == Qt.RightButton:
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

    def lineMovedFinished(self):
        self.sig_IPRegion_Change_finished.emit(self.getRegion())
        super().lineMovedFinished()



class IPLinearRegionItem_Signal(LinearRegionItem):

    sig_IPRegion_Change_finished = pyqtSignal(tuple)

    def __init__(self, values=[0, 1], orientation=LinearRegionItem.Vertical, brush=None, movable=True, bounds=None, swapMode='block'):
        super().__init__(values=values, orientation=orientation, brush=brush, movable=movable, bounds=bounds, swapMode=swapMode)
        self.setZValue(15)
        brush = QBrush(QColor(80, 159, 250, 100))
        self.setBrush(brush)

    def mouseClickEvent(self, ev):

        if ev.button() == Qt.RightButton:
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

    def lineMovedFinished(self):
        print("move finished")
        self.sig_IPRegion_Change_finished.emit(self.getRegion())
        super().lineMovedFinished()


class IPFreqLinearRegionItem(LinearRegionItem):

    def __init__(self, values=[np.log10(0.5), np.log10(5)], orientation=LinearRegionItem.Vertical, brush=None, movable=True, bounds=None, swapMode='block'):
        super().__init__(values=values, orientation=orientation, brush=brush, movable=movable, bounds=bounds, swapMode=swapMode)
        self.setZValue(15)
        brush = QBrush(QColor(200, 200, 200, 50))
        self.setBrush(brush)

    def mouseClickEvent(self, ev):

        if ev.button() == Qt.RightButton:
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
