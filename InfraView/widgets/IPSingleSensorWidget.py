from PyQt5.QtWidgets import (QWidget, QVBoxLayout)

from PyQt5.QtCore import pyqtSlot
from PyQt5 import QtGui

import pyqtgraph as pg
import numpy as np
from scipy import signal
from obspy.core import UTCDateTime

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

        self.signalSpecWidget = IPSpectrogramWidget(self)
        self.noiseSpecWidget = IPSpectrogramWidget(self)

        self.addItem(self.waveformPlot)
        self.nextRow()
        self.addItem(self.signalSpecWidget)
        self.nextRow()
        self.addItem(self.noisePlot)
        self.nextRow()
        self.addItem(self.noiseSpecWidget)

    def get_earliest_start_time(self):
        return self.parent.waveformWidget.get_earliest_start_time()

    @pyqtSlot(object)
    def signal_region_changed(self, lri):
        print('signal region changed {}'.format(lri.getRegion()))

    @pyqtSlot(object)
    def noise_region_changed(self, lri):
        print('noise region changed {}'.format(lri.getRegion()))

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
    def updateSignalRange(self, new_range):
        self.waveformPlot.setXRange(new_range[0], new_range[1], padding=0)
        # we want to set the title of the plot to reflect the current start time of the view
        self.start_time = self.get_earliest_start_time() + new_range[0]
        self.waveformPlot.setTitle(str(self.start_time) + " (Signal)")

        self.signalSpecWidget.set_xaxis(new_range, self.start_time)
        self.signalSpecWidget.auto_scale_yaxis()


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
    def updateNoiseRange(self, new_range):
        self.noisePlot.setXRange(new_range[0], new_range[1], padding=0)
        # we want to set the title of the plot to reflect the current start time of the view
        self.start_time = self.get_earliest_start_time() + new_range[0]
        self.noisePlot.setTitle(str(self.start_time) + " (Background)")

        self.noiseSpecWidget.set_xaxis(new_range, self.start_time)
        self.noiseSpecWidget.auto_scale_yaxis()
        
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

    transform = None
    full_range_y = []
    color_bar = None

    def __init__(self, parent, est=None):
        super().__init__(mode='spectrogram')
        self.parent = parent
        self.buildUI()

    def buildUI(self):
        self.spec_img = pg.ImageItem( image=np.eye(3), levels=(0,1) ) # create example image
        self.addItem(self.spec_img)

    def set_data(self, data, region):
        # data is a plot data item, range is the current range viewed
        self.data_item = data
        self.region = region

    def set_xaxis(self, range, start_time):
        self.setXRange(range[0], range[1], padding=0)
        
        # The tics will need to be adjusted for the new start time...
        self.getAxis('bottom').set_start_time(start_time)
    
    def set_yaxis(self, range):
        self.setYRange(range[0], range[1], padding=0)

    def auto_scale_yaxis(self):
        self.set_yaxis([self.full_range_y[0], self.full_range_y[1]])

    def calc_spectrogram(self, data, nfft=1024*5, Fs=1.0, noverlap=100, nperseg=128):
        if data is None:
            return

        f, t, Sxx = signal.spectrogram(data[1], Fs, nfft=nfft, noverlap=noverlap, nperseg=nperseg)
        log_Sxx = 10 * np.log10(Sxx)
        
        # the spectrogram comes in with the time axis as the y-axis... transpose it for sane plotting
        log_Sxx = np.transpose(log_Sxx)

        self.plot_spectrogram(f, t, log_Sxx)

    def plot_spectrogram(self, f, t, Sxx):

        self.transform = QtGui.QTransform()
        yscale = f[-1]/Sxx.shape[1]
        xscale = t[-1]/Sxx.shape[0]
        self.transform.scale(xscale, yscale)
        self.spec_img.setTransform(self.transform)

        # colormap
        color_maps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis']
        cmap = pg.colormap.get(color_maps[2])

        minv, maxv = np.nanmin(np.nanmin(Sxx[Sxx != -np.inf])), np.nanmax(np.nanmax(Sxx))

        if self.color_bar is None:
            self.color_bar = pg.ColorBarItem(interactive=True, values=(minv, maxv), colorMap=cmap, label='magnitude [dB]')
            self.color_bar.setImageItem(self.spec_img, insert_in=self)
        else:
            self.color_bar.setLevels(low=minv, high=maxv)
            
        self.spec_img.setColorMap(cmap)
        
        self.spec_img.setImage(Sxx)

        self.full_range_y = [f[0], f[-1]]
        self.set_yaxis(self.full_range_y)


    

        
