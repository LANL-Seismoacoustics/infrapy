from PyQt5.QtWidgets import (QWidget, QVBoxLayout)

from PyQt5.QtCore import pyqtSignal, pyqtSlot
from PyQt5 import QtCore, QtGui

import pyqtgraph as pg
import numpy as np
from scipy import signal

from InfraView.widgets import IPPlotItem

class IPSingleSensorWidget(pg.GraphicsLayoutWidget):

    waveform_data_item = None
    noise_data_item = None

    signal_fs = 1.0
    noise_fs = 1.0
    
    def __init__(self, parent=None):
        super().__init__(parent)

        self.parent = parent
        self.buildUI()

    def buildUI(self):
        self.lhWidget = pg.GraphicsLayoutWidget()

        self.waveformPlot = IPPlotItem.IPPlotItem(mode='waveform', est=None)
        self.waveformPlot.setLabel('left', 'Amplitude')
        self.waveformPlot.hideButtons()

        self.noisePlot = IPPlotItem.IPPlotItem(mode='waveform', est=None)
        self.noisePlot.setLabel('left', 'Amplitude')
        self.noisePlot.hideButtons()

        self.signalSpecWidget = IPSpectrogramWidget()
        self.noiseSpecWidget = IPSpectrogramWidget()

        self.addItem(self.waveformPlot)
        self.nextRow()
        self.addItem(self.noisePlot)
        self.nextRow()
        self.addItem(self.signalSpecWidget)
        self.nextRow()
        self.addItem(self.noiseSpecWidget)

    def get_earliest_start_time(self):
        return self.parent.waveformWidget.plotViewer.pl_widget.earliest_start_time

    #@pyqtSlot(pg.LinearRegionItem)
    def signal_region_changed(self):
        print('signal region changed')

   # @pyqtSlot(pg.LinearRegionItem)
    def noise_region_changed(self):
        print('noise region changed')

    @pyqtSlot(pg.PlotDataItem, tuple, str)
    def setSignalWaveform(self, plotLine, region, plot_label=None):
        # pretty much the same as the setWaveform in IPBeamformingWidget

        initial = False
        if self.waveform_data_item is not None:
            self.waveform_data_item.clear()
        else:
            self.waveform_data_item = pg.PlotDataItem()
            initial = True

        # bringing in a new waveform, we might have a new earliest_start_time, so update that in the 
        # plots so that the x-axes will be correct
        self.waveformPlot.setEarliestStartTime(self.get_earliest_start_time())

        # need to make a copy of the currently active plot and give it to the beamformingwidget for display
        self.waveform_data_item.setData(plotLine.xData, plotLine.yData)
        self.waveform_data_item.setPen(pg.mkPen(color=(100, 100, 100), width=1))
        self.waveformPlot.enableAutoRange(axis=pg.ViewBox.YAxis)

        # calculate the sampling frequency
        self.signal_fs = 1.0/(plotLine.xData[1] - plotLine.xData[0])

        if initial:
            # only need to add the item if it wasn't already added
            self.waveformPlot.addItem(self.waveform_data_item)
        if plot_label is not None:
            self.waveformPlot.setPlotLabel(plot_label)
        self.waveformPlot.setXRange(region[0], region[1], padding=0)
        
        # generate spectrogram
        if self.waveform_data_item is not None:
            self.signalSpecWidget.calc_spectrogram(self.waveform_data_item.getData(), Fs=self.signal_fs )

    @pyqtSlot(tuple)
    def updateSignalWaveformRange(self, new_range):
        self.waveformPlot.setXRange(new_range[0], new_range[1], padding=0)
        # we want to set the title of the plot to reflect the current start time of the view
        self.start_time = self.get_earliest_start_time() + new_range[0]
        self.waveformPlot.setTitle(str(self.start_time) + " (Signal)")

        # generate spectrogram
        if self.waveform_data_item is not None:
            self.signalSpecWidget.calc_spectrogram(self.waveform_data_item.getData(), Fs=self.signal_fs)

    @pyqtSlot(pg.PlotDataItem, tuple, str)
    def setNoiseWaveform(self, plotLine, region, plot_label=None):
        initial = False
        if self.noise_data_item is not None:
            self.noise_data_item.clear()
        else:
            self.noise_data_item = pg.PlotDataItem()
            initial = True

        # bringing in a new waveform, we might have a new earliest_start_time, so update that in the 
        # plots so that the x-axes will be correct
        self.noisePlot.setEarliestStartTime(self.get_earliest_start_time())

        # need to make a copy of the currently active plot and give it to the beamformingwidget for display
        self.noise_data_item.setData(plotLine.xData, plotLine.yData)
        self.noise_data_item.setPen(pg.mkPen(color=(100, 100, 100), width=1))
        self.noisePlot.enableAutoRange(axis=pg.ViewBox.YAxis)

        # calculate the sampling frequency
        self.noise_fs = 1.0/(plotLine.xData[1] - plotLine.xData[0])

        if initial:
            # only need to add the item if it wasn't already added
            self.noisePlot.addItem(self.noise_data_item)
        if plot_label is not None:
            self.noisePlot.setPlotLabel(plot_label)
        self.noisePlot.setXRange(region[0], region[1], padding=0)

        # generate spectrogram
        if self.noise_data_item is not None:
            self.noiseSpecWidget.calc_spectrogram(self.noise_data_item.getData(), Fs=self.noise_fs)

    @pyqtSlot(tuple)
    def updateNoiseWaveformRange(self, new_range):
        self.noisePlot.setXRange(new_range[0], new_range[1], padding=0)
        # we want to set the title of the plot to reflect the current start time of the view
        self.start_time = self.get_earliest_start_time() + new_range[0]
        self.noisePlot.setTitle(str(self.start_time) + " (Background)")
        
        # generate spectrogram
        if self.noise_data_item is not None:
            self.noiseSpecWidget.calc_spectrogram(self.noise_data_item.getData(), Fs=self.noise_fs)

    def clearWaveformPlots(self):
        self.waveform_data_item = None
        self.waveformPlot.clear()
        self.waveformPlot.setTitle(None)
        self.waveformPlot.clearPlotLabel()
        self.waveformPlot.setYRange(0, 1, padding=0)
        
        self.noise_data_item = None
        self.noisePlot.clear()
        self.noisePlot.setTitle(None)
        self.noisePlot.clearPlotLabel()
        self.noisePlot.setYRange(0,1,padding=0)
        


class IPSpectrogramWidget(IPPlotItem.IPPlotItem):

    def __init__(self, parent=None):
        super().__init__()

        self.buildUI()

    def buildUI(self):
        #self.spec_img = pg.ImageItem()
        self.spec_img = pg.ImageItem( image=np.eye(3), levels=(0,1) ) # create example image
        self.addItem(self.spec_img)

    def set_data(self, data, region):
        # data is a plot data item, range is the current range viewed
        self.data_item = data
        self.region = region

    def calc_spectrogram(self, data, nfft=1024, Fs=1.0):
        if data is None:
            return

        f, t, Sxx = signal.spectrogram(data[1], Fs, nfft=nfft)
        log_Sxx = np.log10(Sxx)
        print("f = {}".format(f))
        print("t = {}".format(t))
    
        #pg.setConfigOption('imageAxisOrder', 'col-major')

        self.spec_img.setImage(log_Sxx)

        # Colourmap
        # cmap = pg.colormap.get('CET-L9')
        cmap = cmap = pg.colormap.get('viridis')  # matplotlib style
        minv, maxv = np.nanmin(np.nanmin(log_Sxx[log_Sxx != -np.inf])), np.nanmax(np.nanmax(log_Sxx))
        bar = bar = pg.ColorBarItem(interactive=True, values=(minv, maxv), colorMap=cmap, label='magnitude [dB]')
        bar.setImageItem(self.spec_img, insert_in=self)
