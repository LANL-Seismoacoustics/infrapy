
from PyQt5.QtWidgets import (QWidget, QRadioButton, QCheckBox, QComboBox, QDoubleSpinBox, QSpinBox, QFormLayout, 
                             QGroupBox, QHBoxLayout, QLabel, QPushButton, QToolBar, QToolButton,
                             QVBoxLayout, QDialog)

from PyQt5.QtCore import pyqtSlot, pyqtSignal, QObject, QThread
from PyQt5 import QtGui

import pyqtgraph as pg

import numpy as np

import pywt
from scipy.signal import spectrogram, stft

from sklearn.cluster import DBSCAN
from obspy.core import UTCDateTime

from InfraView.widgets import IPPlotItem
from InfraView.widgets import IPUtils
from InfraView.widgets import IPBaseWidgets

from infrapy.detection import spectral

import traceback

class IPSingleSensorWidget(QWidget):

    waveform_data_item = None

    spec_overlap = 0.8
    fs = 1.0

    mp_pool = None

    signal_start_sd_calc = pyqtSignal()
    
    def __init__(self, parent, pool=None):
        super().__init__(parent)
        self.appWidget = parent
        self.mp_pool = pool

        self.buildUI()

    def set_controlling_widget(self, cw):
        self.spectrogram_settings_widget = cw
        self.connect_signals_and_slots()
        self.update_values()

    def connect_signals_and_slots(self):
        self.spectrogram_settings_widget.colormap_cb.currentTextChanged.connect(self.signalSpecWidget.set_colormap)
        self.spectrogram_settings_widget.colorbar_rb.toggled.connect(self.signalSpecWidget.show_hide_colorbar)
        self.spectrogram_settings_widget.update_button.clicked.connect(self.updateSpectrograms)

    def buildUI(self):
        main_layout = QVBoxLayout()

        ##### TOOLBAR
        self.toolbar = QToolBar()
        
        self.tool_runDetector_button = QToolButton()
        self.tool_runDetector_button.setText("Run Detector")
        self.tool_runDetector_button.clicked.connect(self.run_spectral_detector)

        self.tool_clearDetections_button = QToolButton()
        self.tool_clearDetections_button.setText("Clear Detections")
        self.tool_clearDetections_button.clicked.connect(self.clear_detection_plot)

        self.toolbar.addWidget(self.tool_runDetector_button)
        self.toolbar.addWidget(self.tool_clearDetections_button)

        ##### WAVEFORM PLOTS
        self.waveformPlot = IPPlotItem.IPPlotItem(mode='waveform', est=None, lris=False)
        self.waveformPlot.setLabel('left', 'Amplitude')
        self.waveformPlot.hideButtons()
        self.waveformPlot.setPlotLabel("Signal Waveform")

        ##### SPECTROGRAM PLOTS
        self.signalSpecWidget = IPSpectrogramWidget(self)
        self.signalSpecWidget.setPlotLabel('Signal')

        ##### DETECTION PLOT
        self.detectionPlot = IPDetectionPlotItem(self, start_time=None)
        

        # link all the plot x-axes so that when they rescale together
        self.signalSpecWidget.setXLink(self.waveformPlot)
        self.detectionPlot.setXLink(self.waveformPlot)

        ##### LAYOUT
        self.glWidget = pg.GraphicsLayoutWidget()
        self.glWidget.addItem(self.waveformPlot)
        self.glWidget.nextRow()
        self.glWidget.addItem(self.signalSpecWidget)
        self.glWidget.nextRow()
        self.glWidget.addItem(self.detectionPlot)

        main_layout.setMenuBar(self.toolbar)
        main_layout.addWidget(self.glWidget)

        self.setLayout(main_layout)

        self.sdThread = QThread()

    @pyqtSlot(str)
    def update_theme(self, t):
        if t == 'light':
            self.glWidget.setBackground((255,255,255))
        elif t == 'dark':
            self.glWidget.setBackground(IPUtils.ip_dark_grey)

    @pyqtSlot()
    def update_values(self):

        new_fmin = self.spectrogram_settings_widget.fmin_spin.value()
        self.spec_overlap = 0.8

        self.nperseg = int(6. * self.fs / new_fmin)
        if self.nperseg > 512:
            self.nperseg = 512

        self.noverlap = int(0.8 * self.nperseg)
        self.nfft = self.nperseg


    def get_earliest_start_time(self):
        return self.appWidget.waveformWidget.get_earliest_start_time()
    
    def run_spectral_detector(self):
        self.detectionPlot.clear()

        f_s, t_s, Sxx_log = self.signalSpecWidget.get_logdata()
        # Pull in the detector settings
        if self.spectrogram_settings_widget.spec_type_cb.currentText() == 'CWT':
            # pull info from CWT widget
            t_skip = int(self.nperseg * (1.0 - self.spec_overlap))
            pval = self.spectrogram_settings_widget.cwt_pval_spin.value()
            freq_band = [self.spectrogram_settings_widget.cwt_fmin_spin.value(), 
                         self.spectrogram_settings_widget.cwt_fmax_spin.value()]
            clustering_freq_scaling = self.spectrogram_settings_widget.cwt_clust_freq_scale_spin.value()
            clustering_eps = self.spectrogram_settings_widget.cwt_clust_eps_spin.value()
            clustering_min_samples = self.spectrogram_settings_widget.cwt_clust_min_samples_spin.value()

            adaptive_window_length = self.spectrogram_settings_widget.cwt_adaptive_win_len_spin.value()

        else:
            # pull settings from spectrogram/stft widget
            t_skip = 1
            pval = self.spectrogram_settings_widget.pval_spin.value()
            freq_band = [self.spectrogram_settings_widget.fmin_spin.value(), 
                         self.spectrogram_settings_widget.fmax_spin.value()]
         
            clustering_freq_scaling = self.spectrogram_settings_widget.clust_freq_scale_spin.value()
            clustering_eps = self.spectrogram_settings_widget.clust_eps_spin.value()
            clustering_min_samples = self.spectrogram_settings_widget.clust_min_samples_spin.value()

            adaptive_window_length = self.spectrogram_settings_widget.adaptive_win_len_spin.value()

        adaptive_window_step = adaptive_window_length/2

        signal_t_range = self.signalSpecWidget.get_xrange()
        signal_window_mask = np.logical_and(signal_t_range[0] <= t_s, t_s <= signal_t_range[1])
        signal_Sxx_window = Sxx_log[:, signal_window_mask]

        

        self.signal_t_window = t_s[signal_window_mask]
        self.t_s = t_s
        self.f_s = f_s

        self.sd_worker_object = IPSpectralDetectorWorkerObject(self.f_s, self.signal_t_window, 
                                                                signal_Sxx_window, freq_band, 
                                                                pval, adaptive_window_length , 
                                                                adaptive_window_step, clustering_freq_scaling, 
                                                                clustering_eps, clustering_min_samples,
                                                                self.mp_pool, t_skip)
        
        self.signal_start_sd_calc.connect(self.sd_worker_object.run)
        self.sd_worker_object.moveToThread(self.sdThread)
        self.sd_worker_object.signal_runFinished.connect(self.detection_run_finished)
        self.sdThread.start()
        self.signal_start_sd_calc.emit()

    @pyqtSlot(np.ndarray, list, str)
    def detection_run_finished(self, spec_dets, clustering_pts, err_msg):
        if err_msg:
            # didn't find anything exit tread.
            self.sdThread.exit()
            IPUtils.errorPopup(err_msg)
        else:
            self.detectionPlot.plot_data(spec_dets, 
                                    self.signal_t_window[1]-self.signal_t_window[0], 
                                    self.waveformPlot.get_start_time(), 
                                    [self.t_s[0],self.t_s[-1]], 
                                    [self.f_s[0], self.f_s[-1]],
                                    clustering_pts)

        
    @pyqtSlot(object)
    def signal_region_changed(self, lri):
        #print('signal region changed {}'.format(lri.getRegion()))
        pass

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
        new_fs = 1.0/(plotLine.xData[1] - plotLine.xData[0])
        if new_fs != self.fs:
            self.fs = new_fs
            self.update_values()

        if initial:
            # only need to add the item if it wasn't already added
            self.waveformPlot.addItem(self.waveform_data_item)
        if plot_label is not None:
            self.waveformPlot.setPlotLabel(plot_label)
        self.waveformPlot.setXRange(region[0], region[1], padding=0)
        self.updateSignalSpectrogram()

    @pyqtSlot(tuple)
    def updateSignalRange(self, new_range):
        self.waveformPlot.setXRange(new_range[0], new_range[1], padding=0)
        # we want to set the title of the plot to reflect the current start time of the view
        self.start_time = self.get_earliest_start_time() + new_range[0]
        self.waveformPlot.setTitle(str(self.start_time) + " (Signal)")

        self.signalSpecWidget.set_xaxis(new_range)
        self.signalSpecWidget.auto_scale_yaxis()

    def updateSignalSpectrogram(self):
         # generate spectrogram
        if self.waveform_data_item is not None:
            # calculate the values used in the spectrograms.  These will be used in a few different places
            self.signalSpecWidget.calc_spectrogram(self.waveform_data_item.getOriginalDataset(), 
                                                  Fs=self.fs, 
                                                  nfft=self.nfft,
                                                  nperseg=self.nperseg,
                                                  noverlap=self.noverlap,
                                                  spec_type=self.spectrogram_settings_widget.spec_type_cb.currentText(),
                                                  morlet_o=self.spectrogram_settings_widget.omega0_spin.value())
            self.signalSpecWidget.set_start_time(self.get_earliest_start_time())

    
    @pyqtSlot()
    def updateSpectrograms(self):
        if self.spectrogram_settings_widget.spec_type_cb.currentText() != self.spectrogram_settings_widget.last_spec_type:
            self.spectrogram_settings_widget.last_spec_type = self.spectrogram_settings_widget.spec_type_cb.currentText()
            self.detectionPlot.spi.clear()
        self.updateSignalSpectrogram()
        
    def clearWaveformPlots(self):
        self.waveform_data_item = None
        self.waveformPlot.clear()
        self.waveformPlot.setTitle("")
        self.waveformPlot.clearPlotLabel()
        self.waveformPlot.setYRange(0, 1, padding=0)

        self.signalSpecWidget.clear_spectrogram()
        
        self.clear_detection_plot()

    def clear_detection_plot(self):
        self.detectionPlot.clear()
        self.detectionPlot.spi.clear()
        self.detectionPlot.setXLink(self.waveformPlot)

class IPSpectralDetectorWorkerObject(QObject):

    signal_runFinished = pyqtSignal(np.ndarray, list, str)

    def __init__(self, 
                 f_s, 
                 signal_t_window, 
                 signal_Sxx_window, 
                 freq_band, 
                 pval, 
                 adaptive_window_length , 
                 adaptive_window_step, 
                 clustering_freq_scaling, 
                 clustering_eps, 
                 clustering_min_samples, 
                 mp_pool, 
                 t_skip):
        
        super().__init__()
        self.f_s = f_s
        self.signal_t_window = signal_t_window
        self.signal_Sxx_window = signal_Sxx_window
        self.freq_band = freq_band
        self.pval = pval
        self.adaptive_window_length = adaptive_window_length
        self.adaptive_window_step = adaptive_window_step
        self.clustering_freq_scaling = clustering_freq_scaling
        self.clustering_eps = clustering_eps
        self.clustering_min_samples = clustering_min_samples
        self.mp_pool = mp_pool
        self.t_skip = t_skip

        self.thread_stopped = True

    @pyqtSlot()
    def run(self):
        self.thread_stopped = False

        try:
            # NOTE: the clustering window is hardcoded here to 600.  TODO read this from config file?
            clustering_window_len = 600
            spec_dets = np.empty(0)
            clustering_pts = []
            err_msg = ""
            spec_dets, clustering_pts, _ = spectral.run_sd(self.f_s, 
                                                           self.signal_t_window, 
                                                           self.signal_Sxx_window, 
                                                            self.freq_band, 
                                                            self.pval, 
                                                            self.adaptive_window_length , 
                                                            self.adaptive_window_step, 
                                                            self.clustering_freq_scaling, 
                                                            self.clustering_eps, 
                                                            self.clustering_min_samples, 
                                                            clustering_window_len, 
                                                            self.mp_pool, 
                                                            self.t_skip)
            

        except ValueError:
            err_msg = "Found no detections.  If you believe something is there, try expanding your analysis window to get better statistics for the background calculation."
        except Exception as e:
            err_msg = f"Error calculating detections:\n\t{str(e)}"

        self.signal_runFinished.emit(spec_dets, clustering_pts, err_msg)
        self.stop()            

    @pyqtSlot()
    def stop(self):
        self.thread_stopped = True


class IPSpectrogramWidget(IPPlotItem.IPPlotItem):

    sig_start_spec_calc = pyqtSignal()
    sig_start_stft_calc = pyqtSignal()

    def __init__(self, parent, est=None):
        super().__init__(mode='spectrogram')
        self.singleStationWidget = parent

        self.labi = None
        self.transform = None
        self.full_range_y = None
        self.spec_img = None     # image that holds the spectrogram
        self.histogram = None

        self.f = None
        self.t = None
        self.Sxx = None

        self.buildUI()

    def buildUI(self):
        self.setLabel(axis='left', text='Frequency (Hz)')
        #self.setLabel(axis='bottom', text='Time') # probably not needed, i think everyone knows what this axis is
        self.spec_img = pg.ImageItem( image=np.eye(3), levels=(0,1) ) # create example image
        self.addItem(self.spec_img)

        self.color_bar = pg.ColorBarItem()
        self.color_bar.setImageItem(self.spec_img, insert_in=self)
        self.color_bar.setVisible(False)

        # horizontal lines indicating filter min and max
        # initial_pos = self.singleStationWidget.spectrogram_settings_widget.fmin_spin.value()
        self.filter_min_line = pg.InfiniteLine(angle=0, pen='k', label="")
        self.filter_min_line.label.setPosition(0.9)
        self.filter_min_line.label.setColor((128,128,128))

        # initial_pos = self.singleStationWidget.spectrogram_settings_widget.fmax_spin.value()
        self.filter_max_line = pg.InfiniteLine(angle=0, pen='k', label="")
        self.filter_max_line.label.setPosition(0.9)
        self.filter_max_line.label.setColor((128,128,128))

        self.spec_img.getViewBox().addItem(self.filter_min_line)
        self.spec_img.getViewBox().addItem(self.filter_max_line)

        self.calc_spec_thread = QThread()

    @pyqtSlot(float)
    def fmin_changed(self, new_val):
        self.filter_min_line.setValue(new_val)
        self.filter_min_line.label.setText(f'fmin = {new_val}')

    @pyqtSlot(float)
    def fmax_changed(self, new_val):
        self.filter_max_line.setValue(new_val)
        self.filter_max_line.label.setText(f'fmax = {new_val}')

    def setPlotLabel(self, text):
        if self.labi is not None:
            self.vb.removeItem(self.labi)
        self.labi = pg.LabelItem(text=text)
        self.labi.setParentItem(self.vb)
        self.labi.anchor(itemPos=(0.0,0.0), parentPos=(0.0,0.0))

    def set_data(self, data, region):
        # data is a plot data item, range is the current range viewed
        self.data_item = data
        self.region = region

    def get_logdata(self):
        # return the f,t,and log10(Sxx) data for further use
        if self.Sxx is None:
            return None, None, None
        return self.f, self.t, 10 * np.log10(self.Sxx)
    
    def get_xrange(self):
        return self.viewRange()[0]
    
    def get_yrange(self):
        return self.viewRange()[1]
    
    def get_fmax(self):
        # needed by the detector, we need to return the maximum frequency of the spectrogram
        return self.f[-1]

    def set_xaxis(self, range=None):
        if range is None:
            self.setXRange(0, 1, padding=0)
        else:
            self.setXRange(range[0], range[1], padding=0)
    
    def set_yaxis(self, range=None):
        if range is None:
            self.setYRange(0, 1, padding=0)
        else:
            self.setYRange(range[0], range[1], padding=0)

    def auto_scale_yaxis(self):
        if self.full_range_y is not None:
            self.set_yaxis([self.full_range_y[0], self.full_range_y[1]])

    def clear_spectrogram(self):
        self.f = None
        self.t = None 
        self.Sxx = None
        self.spec_img = pg.ImageItem( image=np.eye(3), levels=(0,1) ) # create example image
        self.clear()
        self.addItem(self.spec_img)

    def calc_spectrogram(self, data, nfft, Fs, noverlap, nperseg, spec_type, morlet_o=None):
        if data is None:
            return
        
        self.calc_spec_worker_object = IPSpectrogramCalcWorker(data=data, fs=Fs, nfft=nfft, nperseg=nperseg, noverlap=noverlap, spec_type=spec_type, morlet_o=morlet_o)
        self.sig_start_spec_calc.connect(self.calc_spec_worker_object.run)
        self.calc_spec_worker_object.signal_runFinished.connect(self.run_finished)
        self.calc_spec_worker_object.moveToThread(self.calc_spec_thread)
        self.calc_spec_thread.start()
        self.sig_start_spec_calc.emit()

    @pyqtSlot(bool)
    def run_finished(self, success):
        if success:
            #keep Sxx in non-log form.  We will take the log of it later as needed
            self.f, self.t, self.Sxx = self.calc_spec_worker_object.get_results()
            if self.f is not None:
                self.plot_spectrogram(self.f, self.t, self.Sxx)
        else:
            IPUtils.errorPopup("Error while calculating the spectrogram")
    
    def set_start_time(self, st):
        self.getAxis(name='bottom').set_start_time(st)

    def get_start_time(self):
        return self.getAxis(name='bottom').get_start_time()
    
    @pyqtSlot(str)
    def set_colormap(self, color_map_str):
        cmap = pg.colormap.get(color_map_str, source='matplotlib')
        self.spec_img.setColorMap(cmap)
        self.color_bar.setColorMap(cmap)

    @pyqtSlot(bool)
    def show_hide_colorbar(self, checked):
        if self.color_bar is not None:
            self.color_bar.setVisible(checked)

    def plot_spectrogram(self, f, t, Sxx):

        # calc the db of the data (ref at 1.0 - if data is in Pascals, we'll want to do 20log10(P/Pref), but not sure how to know that right now)
        Sxx = 10 * np.log10(Sxx)
        Sxx = np.transpose(Sxx)
        Sxx[Sxx==-np.inf] = 0.000001
        #####SET UP THE TRANSFORM
        self.transform = QtGui.QTransform()
        yscale = f[-1]/Sxx.shape[1]
        xscale = t[-1]/Sxx.shape[0]
        self.transform.scale(xscale, yscale)
        self.spec_img.setTransform(self.transform)

        #####COLORMAP
        cmap = self.singleStationWidget.spectrogram_settings_widget.colormap_cb.currentText()
        self.set_colormap(cmap)

        #####SCALE WIDGET
        minv, maxv = np.nanmin(np.nanmin(Sxx[Sxx != -np.inf])), np.nanmax(np.nanmax(Sxx))

        #####BUILD THE COLORBAR
        self.color_bar = pg.ColorBarItem(interactive=True, label='magnitude [dB]', values=(-80,80))
        self.color_bar.setColorMap(cmap)
        self.color_bar.setImageItem(self.spec_img, insert_in=self)
        self.color_bar.setLevels(low=minv, high=maxv)
        colorbar_visible = self.singleStationWidget.spectrogram_settings_widget.colorbar_rb.isChecked()
        self.show_hide_colorbar(colorbar_visible)

        # ADD THE DATA TO MAKE THE IMAGE
        self.spec_img.setImage(Sxx)

        # AXIS LIMITS
        # self.setLimits(xMin=0, xMax=t[-1], yMin=0, yMax=f[-1])
        self.full_range_y = [f[0], f[-1]]
        self.set_yaxis(self.full_range_y)
        


class IPDetectionStatusDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.buildUI()

    def buildUI(self):
        
        self.threshold_label = QLabel("Threshold: ")
        self.detection_label = QLabel("Detections: ")
        main_layout = QVBoxLayout()
        main_layout.addWidget(self.threshold_label)
        main_layout.addWidget(self.detection_label)

        self.setLayout(main_layout)

    def exec_(self):
        super().exec_()

    def calculating_threshold(self):
        self.threshold_label.setText("Theshold: Calculating...")

    def finished_threshold(self):
        self.threshold_label.setText("Threshold: complete") 

    def calculating_detections(self):
        self.detection_label.setText("Detections: Calculating...")

    def finished_detections(self, value):
        self.detection_label.setText("Detections: " + str(value))


class IPScatterPlotTimeAxis(pg.AxisItem):
    # subclass the basic axis item, mainly to make custom time axis
    start_time = UTCDateTime(0)

    def __init__(self, start_time, *args, **kwargs):
        super().__init__(orientation='bottom', *args, **kwargs)
        # est is the "earliest_start_time"
        self.set_start_time(start_time)

    def tickStrings(self, values, scale, spacing):
        return [(self.start_time + value).strftime("%H:%M:%S") for value in values]

    def set_start_time(self, st):
        self.start_time = st

    def get_start_time(self):
        return self.start_time


class IPDetectionPlotItem(pg.PlotItem):
    data_item = None
    def __init__(self, parent, start_time=None):

        if start_time is None:
            start_time = UTCDateTime(0)
        super().__init__(axisItems={'bottom': IPPlotItem.IPSpectrogramTimeAxis()})

        self.showAxis('right')
        self.getAxis('right').setTicks('')
        self.showAxis('top')
        self.getAxis('top').setTicks('')

        # go ahead and create our scatterplotitem here for later use
        self.spi = pg.ScatterPlotItem(pxMode=True)
        self.addItem(self.spi)

        self.getAxis('left').setWidth(80)
        self.setLabel(axis='left', text='Frequency (Hz)')

    def testplot(self):
        # for testing
        spots3 = []
        for i in range(10):
            for j in range(10):
                spots3.append({'pos': (1e-3*i, 1e-3*j), 'size': 1e-3, 'pen': {'color': 'w', 'width': 2}, 'brush':pg.intColor(i*10+j, 100)})
        self.spi.addPoints(spots3)

    def set_yaxis(self, range=None):
        if range is None:
            self.setYRange(0, 1, padding=0)
        else:
            self.setYRange(range[0], range[1], padding=0)

    def clear(self):
        self.removeItem(self.spi)

    def plot_data(self, detections, dt, start_time, t_range, f_range, cluster_pts):
        #detection data comes in as a list of lists with each element being [t,f,value]
        # first clear any old points
        self.spi.clear()

        # initialize the start time of the axis
        self.getAxis('bottom').set_start_time(start_time)

        colors = IPUtils.blue_to_red
        spots = []
        for idx, det_points in enumerate(cluster_pts):
            try:
                col = colors[idx]
            except IndexError:
                IPUtils.errorPopup("More detections than color options. Perhaps run on a smaller timespan, or increase the linkage (EPS).")
                return
            
            for pts in det_points:
                spots.append({'pos': (pts[0], pts[1]), 'pen': {'color': col}, 'brush': col, 'symbol': 's', 'size':5.6*dt})

        self.spi.addPoints(spots)
        self.addItem(self.spi)

        # set axis limits
        self.setLimits(xMin=t_range[0], xMax=t_range[1], yMin=f_range[0], yMax=f_range[1])
        self.full_range_y = [f_range[0], f_range[1]]
        self.set_yaxis(self.full_range_y)


class IPSpectrogramCalcWorker(QObject):
    signal_runFinished = pyqtSignal(bool)

    f = None
    t = None
    Sxx = None

    def __init__(self, data, fs, nfft, nperseg, noverlap, spec_type, morlet_o=None):
        ''' options for type are "Spectrogram", "STFT", or "CWT" '''
        super().__init__()
        self.data = data
        self.fs = fs
        self.nfft = nfft
        self.nperseg = nperseg
        self.noverlap = noverlap
        self.spec_type = spec_type
        self.morlet_omega0 = morlet_o

    @pyqtSlot()
    def run(self):
        if self.spec_type.upper() == 'SPECTROGRAM':
            try:
                self.f, self.t, self.Sxx = spectrogram(self.data[1], 
                                        self.fs, 
                                        nfft=self.nfft, 
                                        noverlap=self.noverlap, 
                                        nperseg=self.nperseg)
            except:
                self.signal_runFinished.emit(False)
                return

        elif self.spec_type.upper() == 'STFT':
            try:
                self.f, self.t, self.Sxx = stft(self.data[1], 
                                            fs=self.fs, 
                                            nperseg=self.nperseg, 
                                            noverlap=self.noverlap, 
                                            nfft=self.nfft)
                
            except:
                self.signal_runFinished.emit(False)
                return 

        elif self.spec_type.upper() == 'CWT':
            if self.morlet_omega0 is None:
                IPUtils.errorPopup("Morelet Omega0 needed to calculate CWT")
                self.signal_runFinished.emit(False)
                return
            try:
                self.f, _, _ = spectrogram(self.data[1], 
                                           self.fs, 
                                           nperseg=self.nperseg, 
                                           noverlap=self.noverlap,
                                           nfft=self.nfft)
                self.t = self.data[0]
        
                '''
                self.f = self.f[1:]
                widths = self.morlet_omega0 / (2 * np.pi * self.f) * self.fs
                self.Sxx = cwt(self.data[1], morlet2, widths, w=self.morlet_omega0)
                '''

                wavelet_name = 'cmor1.0-' + str(self.morlet_omega0 / (2 * np.pi))
                scales = pywt.frequency2scale(wavelet_name, self.f[1:] / self.fs)      

                self.Sxx, self.f = pywt.cwt(self.data[1],
                                            scales,
                                            wavelet_name,
                                            sampling_period=1.0 / self.fs)

            except:
                self.signal_runFinished.emit(False)

        self.signal_runFinished.emit(True)


    def get_results(self):
        if self.Sxx is None:
            return None, None, None
        return self.f, self.t, np.abs(self.Sxx)

        