
from PyQt5.QtWidgets import (QWidget, QRadioButton, QCheckBox, QComboBox, QDoubleSpinBox, QSpinBox, QFormLayout, 
                             QGroupBox, QHBoxLayout, QLabel, QPushButton, QToolBar, QToolButton,
                             QVBoxLayout, QDialog)

from PyQt5.QtCore import pyqtSlot, pyqtSignal, QObject, QThread
from PyQt5 import QtGui

import pyqtgraph as pg

import numpy as np

from scipy.signal import spectrogram, stft, cwt, morlet2

from sklearn.cluster import DBSCAN
from obspy.core import UTCDateTime

from InfraView.widgets import IPPlotItem
from InfraView.widgets import IPUtils

from infrapy.detection import spectral

class IPSingleSensorWidget(QWidget):

    waveform_data_item = None
    noise_data_item = None

    spec_overlap = 0.8
    fs = 1.0

    mp_pool = None
    
    def __init__(self, parent, pool=None):
        super().__init__(parent)
        self.appWidget = parent
        self.mp_pool = pool

        self.buildUI()

        self.update_values()


    def buildUI(self):
        main_layout = QVBoxLayout()

        ##### TOOLBAR
        self.toolbar = QToolBar()
        self.toolbar.setStyleSheet("QToolBar{background-color: silver; }")
        #self.toolbar.setStyleSheet("QToolBar { border-bottom: 1px solid; } ")
        self.tool_runDetector_button = QToolButton()
        self.tool_runDetector_button.setText("Run Detector")
        self.tool_runDetector_button.clicked.connect(self.run_spectral_detector)

        self.tool_settings_button = QToolButton()
        self.tool_settings_button.setText("Settings...")
        self.tool_settings_button.clicked.connect(self.show_hide_spectrogram_settings)


        self.toolbar.addWidget(self.tool_runDetector_button)
        self.toolbar.addWidget(self.tool_settings_button)

        ##### SETTINGS WIDGET
        self.spectrogram_settings_widget = IPSpectrogramSettingsWidget(self)
        self.spectrogram_settings_widget.setVisible(False)

        ##### WAVEFORM PLOTS
        self.waveformPlot = IPPlotItem.IPPlotItem(mode='waveform', est=None, lris=False)
        self.waveformPlot.setLabel('left', 'Amplitude')
        self.waveformPlot.hideButtons()
        self.waveformPlot.setPlotLabel("Signal Waveform")
        self.waveformPlot.setVisible(self.spectrogram_settings_widget.show_signal_waveform_cb.isChecked())

        self.noisePlot = IPPlotItem.IPPlotItem(mode='waveform', est=None, lris=False)
        self.noisePlot.setLabel('left', 'Amplitude')
        self.noisePlot.hideButtons()
        self.noisePlot.setPlotLabel("Background Waveform")
        self.noisePlot.setVisible(self.spectrogram_settings_widget.show_noise_waveform_cb.isChecked())

        ##### SPECTROGRAM PLOTS
        self.signalSpecWidget = IPSpectrogramWidget(self)
        self.signalSpecWidget.setPlotLabel('Signal')
        #self.signalSpecWidget.sig_fmax_changed.connect(self.spectrogram_settings_widget.fmax_spin.setValue)
        #self.signalSpecWidget.sig_fmax_changed.connect(self.spectrogram_settings_widget.fmax_spin.setMaximum)
        self.signalSpecWidget.setVisible(self.spectrogram_settings_widget.show_signal_spectrogram_cb.isChecked())

        self.noiseSpecWidget = IPSpectrogramWidget(self)
        self.noiseSpecWidget.setPlotLabel("Background")
        self.noiseSpecWidget.setVisible(self.spectrogram_settings_widget.show_noise_spectrogram_cb.isChecked())

        ##### DETECTION PLOT
        self.detectionPlot = IPDetectionPlotItem(self, start_time=None)
        

        # link all the plot x-axes so that when they rescale together
        self.signalSpecWidget.setXLink(self.waveformPlot)
        self.noiseSpecWidget.setXLink(self.noisePlot)
        self.detectionPlot.setXLink(self.waveformPlot)

        ##### LAYOUT
        glWidget = pg.GraphicsLayoutWidget()
        glWidget.addItem(self.waveformPlot)
        glWidget.nextRow()
        glWidget.addItem(self.signalSpecWidget)
        glWidget.nextRow()
        glWidget.addItem(self.noisePlot)
        glWidget.nextRow()
        glWidget.addItem(self.noiseSpecWidget)
        glWidget.nextRow()
        glWidget.addItem(self.detectionPlot)

        main_layout.setMenuBar(self.toolbar)
        main_layout.addWidget(self.spectrogram_settings_widget)
        main_layout.addWidget(glWidget)

        self.setLayout(main_layout)

    @pyqtSlot(float)
    def update_values(self):
        # self.nperseg = 2048
        # self.noverlap = self.nperseg * self.spec_overlap
        # self.nfft = self.nperseg

        new_fmin = self.spectrogram_settings_widget.fmin_spin.value()
        self.spec_overlap = 0.8

        self.nperseg = int(5. * self.fs / new_fmin)
        if self.nperseg > 512:
            self.nperseg = 512

        self.noverlap = int(0.8 * self.nperseg)

        self.nfft = self.nperseg

    def show_hide_spectrogram_settings(self):
        self.spectrogram_settings_widget.setVisible(self.spectrogram_settings_widget.isHidden())

    def get_earliest_start_time(self):
        return self.appWidget.waveformWidget.get_earliest_start_time()
    
    def run_spectral_detector(self):
        # before we do anything, pull in spectrogram data and make sure there is something to process
        # pull in the spectrogram data
        f_n, t_n, Sxx_log_n = self.noiseSpecWidget.get_logdata()
        if f_n is None or t_n is None or Sxx_log_n is None:
            IPUtils.errorPopup("You must have data loaded to run the detector.")
            return      
        
        f_s, t_s, Sxx_log_s = self.signalSpecWidget.get_logdata()
        if f_s is None or t_s is None or Sxx_log_s is None:
            IPUtils.errorPopup("You must have signal data loaded to run the detector.")
            return

        # Pull in the detector settings
        pval = self.spectrogram_settings_widget.pval_spin.value()
        freq_band = [self.spectrogram_settings_widget.fmin_spin.value(), self.spectrogram_settings_widget.fmax_spin.value()] 
        clustering_freq_scaling = self.spectrogram_settings_widget.clust_freq_scale_spin.value()
        clustering_eps = self.spectrogram_settings_widget.clust_eps_spin.value()
        clustering_min_samples = self.spectrogram_settings_widget.clust_min_samples_spin.value()

        # do a few checks
        if freq_band[0] >= freq_band[1]:
            IPUtils.errorPopup("Frequency min must be less than frequency max")
            return

        noise_t_range = self.noiseSpecWidget.get_xrange()
        noise_window_mask = np.logical_and(noise_t_range[0] <= t_n, t_n <= noise_t_range[1])
        noise_t_window = t_n[noise_window_mask]
        noise_Sxx_window = Sxx_log_n[:, noise_window_mask]

        signal_t_range = self.signalSpecWidget.get_xrange()
        signal_window_mask = np.logical_and(signal_t_range[0] <= t_s, t_s <= signal_t_range[1])
        signal_t_window = t_s[signal_window_mask]
        signal_Sxx_window = Sxx_log_s[:, signal_window_mask]

        # there are way too many points, skip some
        if self.spectrogram_settings_widget.spec_type_cb.currentText() == 'CWT':
            t_skip = int(self.nperseg * (1.0 - self.spec_overlap))
            signal_Sxx_window = signal_Sxx_window[:,::t_skip]
            signal_t_window = signal_t_window[::t_skip]
        else:
            t_skip = 1
            

        freq_band_mask = np.logical_and(freq_band[0] < f_n, f_n < freq_band[1])
        
        if self.mp_pool is not None:
            args = [[noise_Sxx_window[:,::t_skip][fn], pval] for fn in range(len(f_n))]
            temp = self.mp_pool.map(spectral.calc_thresh_wrapper, args)
        else:
            temp = np.array([spectral.calc_thresh(noise_Sxx_window[:,::t_skip][fn], pval) for fn in range(len(f_n))])

        threshold  = np.array(temp)[:,0]
        
        spec_dets = []
        for idx, freq in enumerate(f_s):
            if freq_band[0] < freq and freq < freq_band[1]:
                # spec_dets = spec_dets + [[signal_t_window[tk], freq, signal_Sxx_window[idx][tk]] for tk in range(len(signal_t_window)) if signal_Sxx_window[idx][tk] >= threshold[idx]]
                for tk in range(len(signal_t_window)):
                    if signal_Sxx_window[idx][tk] >= threshold[idx]:
                        spec_dets = spec_dets + [[signal_t_window[tk], freq, signal_Sxx_window[idx][tk]]]

        spec_dets = np.unique(np.array(spec_dets), axis=0)

        # Cluster into detections
        spec_dets_logf = None
        spec_dets_logf = np.stack((spec_dets[:, 0], clustering_freq_scaling * np.log10(spec_dets[:, 1]))).T
        clustering = DBSCAN(eps=clustering_eps, min_samples=clustering_min_samples).fit(spec_dets_logf)
        self.detectionPlot.plot_data(spec_dets, 
                                    signal_t_window[1]-signal_t_window[0], 
                                    self.waveformPlot.get_start_time(), 
                                    [t_s[0],t_s[-1]], 
                                    [f_s[0], f_s[-1]],
                                    clustering)


    @pyqtSlot(object)
    def signal_region_changed(self, lri):
        #print('signal region changed {}'.format(lri.getRegion()))
        pass

    @pyqtSlot(object)
    def noise_region_changed(self, lri):
        #print('noise region changed {}'.format(lri.getRegion()))
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
        self.fs = 1.0/(plotLine.xData[1] - plotLine.xData[0])

        if initial:
            # only need to add the item if it wasn't already added
            self.noisePlot.addItem(self.noise_data_item)
        if plot_label is not None:
            self.noisePlot.setPlotLabel(plot_label)

        self.noisePlot.setXRange(region[0], region[1], padding=0)

        self.updateNoiseSpectrogram()


    @pyqtSlot(tuple)
    def updateNoiseRange(self, new_range):
        self.noisePlot.setXRange(new_range[0], new_range[1], padding=0)
        # we want to set the title of the plot to reflect the current start time of the view
        self.start_time = self.get_earliest_start_time() + new_range[0]
        self.noisePlot.setTitle(str(self.start_time) + " (Background)")

        self.noiseSpecWidget.set_xaxis(new_range)
        self.noiseSpecWidget.auto_scale_yaxis()

    def updateNoiseSpectrogram(self):
        if self.noise_data_item is not None:
            self.noiseSpecWidget.calc_spectrogram(self.noise_data_item.getOriginalDataset(), 
                                                  Fs=self.fs, 
                                                  nfft=self.nfft,
                                                  nperseg=self.nperseg,
                                                  noverlap=self.noverlap,
                                                  spec_type=self.spectrogram_settings_widget.spec_type_cb.currentText(),
                                                  morlet_o=self.spectrogram_settings_widget.omega0_spin.value())

            self.noiseSpecWidget.set_start_time(self.get_earliest_start_time())

    @pyqtSlot()
    def updateSpectrograms(self):
        if self.spectrogram_settings_widget.spec_type_cb.currentText() != self.spectrogram_settings_widget.last_spec_type:
            self.spectrogram_settings_widget.last_spec_type = self.spectrogram_settings_widget.spec_type_cb.currentText()
            self.detectionPlot.spi.clear()
        self.updateSignalSpectrogram()
        self.updateNoiseSpectrogram()

    @pyqtSlot()
    def updateVisibilities(self):
        self.waveformPlot.setVisible(self.spectrogram_settings_widget.show_signal_waveform_cb.isChecked())
        self.signalSpecWidget.setVisible(self.spectrogram_settings_widget.show_signal_spectrogram_cb.isChecked())
        self.noisePlot.setVisible(self.spectrogram_settings_widget.show_noise_waveform_cb.isChecked())
        self.noiseSpecWidget.setVisible(self.spectrogram_settings_widget.show_noise_spectrogram_cb.isChecked())
        
    def clearWaveformPlots(self):
        self.waveform_data_item = None
        self.waveformPlot.clear()
        self.waveformPlot.setTitle("")
        self.waveformPlot.clearPlotLabel()
        self.waveformPlot.setYRange(0, 1, padding=0)
        
        self.noise_data_item = None
        self.noisePlot.clear()
        self.noisePlot.setTitle("")
        self.noisePlot.clearPlotLabel()
        self.noisePlot.setYRange(0,1,padding=0)

        self.signalSpecWidget.clear_spectrogram()
        self.noiseSpecWidget.clear_spectrogram()
        self.detectionPlot.clear()

        self.detectionPlot.spi.clear()


class IPSpectrogramWidget(IPPlotItem.IPPlotItem):

    sig_start_spec_calc = pyqtSignal()
    sig_start_stft_calc = pyqtSignal()
    sig_fmax_changed = pyqtSignal(float)

    def __init__(self, parent, est=None):
        super().__init__(mode='spectrogram')
        self.singleStationWidget = parent

        self.labi = None
        self.transform = None
        self.full_range_y = None
        self.color_bar = None
        self.spec_img = None     # image that holds the spectrogram

        self.color_bar = None
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

        self.calc_spec_thread = QThread()

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
                self.sig_fmax_changed.emit(self.f[-1])
        else:
            IPUtils.errorPopup("Error while calculating the spectrogram")
    
    def set_start_time(self, st):
        self.getAxis(name='bottom').set_start_time(st)

    def get_start_time(self):
        return self.getAxis(name='bottom').get_start_time()

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
        color_map_str = self.singleStationWidget.spectrogram_settings_widget.colormap_cb.currentText()
        cmap = pg.colormap.get(color_map_str, source='matplotlib')
        self.spec_img.setColorMap(cmap)

        #####SCALE WIDGET
        scale_setting = self.singleStationWidget.spectrogram_settings_widget.get_scale_setting()
        minv, maxv = np.nanmin(np.nanmin(Sxx[Sxx != -np.inf])), np.nanmax(np.nanmax(Sxx))

        #####BUILD THE COLORBAR
        if scale_setting == 'cbar':
            if self.color_bar is None:
                self.color_bar = pg.ColorBarItem(interactive=True, label='magnitude [dB]', values=(-80,80))
            self.color_bar.setColorMap(cmap)
            self.color_bar.setImageItem(self.spec_img, insert_in=self)
            self.color_bar.setLevels(low=minv, high=maxv)
            self.color_bar.setVisible(True)

        #####HIDE EVERYTHING
        if scale_setting =='none':
            if self.color_bar is not None:
                self.color_bar.setVisible(False)

        # ADD THE DATA TO MAKE THE IMAGE
        self.spec_img.setImage(Sxx)

        # AXIS LIMITS
        # self.setLimits(xMin=0, xMax=t[-1], yMin=0, yMax=f[-1])
        self.full_range_y = [f[0], f[-1]]
        self.set_yaxis(self.full_range_y)
    


class IPSpectrogramSettingsWidget(QWidget):

    sig_spectrogram_changed = pyqtSignal()
  
    def __init__(self, parent):
        super().__init__(parent)
        self.singleStationWidget = parent
        self.last_spec_type = "Spectrogram"
        self.buildUI()

    def buildUI(self):

        spec_gb = QGroupBox("Spectrogram ")
        spec_layout = QHBoxLayout()

        self.update_button = QPushButton("Update")
        self.update_button.setMaximumWidth(100)
        self.update_button.setEnabled(False)
        self.update_button.clicked.connect(self.deactivate_update_button)
        self.update_button.clicked.connect(self.singleStationWidget.updateSpectrograms)
        self.update_button.clicked.connect(self.singleStationWidget.updateVisibilities)

        self.spec_type_cb = QComboBox(self)
        spec_type_label = QLabel("    Spectrogram type:")
        self.spec_type_cb.addItem("Spectrogram")
        self.spec_type_cb.addItem("STFT")
        self.spec_type_cb.addItem("CWT")
        self.spec_type_cb.currentTextChanged.connect(self.activate_update_button)

        self.omega0_label = QLabel("    omega0")
        self.omega0_spin = QSpinBox()
        self.omega0_spin.setMaximumWidth(100)
        self.omega0_spin.setMinimum(1)
        self.omega0_spin.setMaximum(100)
        self.omega0_spin.setValue(7)
        self.omega0_spin.setToolTip("Used by the CWT.  Number of periods in the wavelet.")
        self.omega0_spin.valueChanged.connect(self.activate_update_button)
        self.omega0_label.setVisible(False)
        self.omega0_spin.setVisible(False)

        colormap_label = QLabel("    Color Map: ")
        self.colormap_cb = QComboBox()
        available_maps = pg.colormap.listMaps(source='matplotlib')
        self.colormap_cb.addItems(available_maps)
        self.colormap_cb.setCurrentText('jet')
        self.colormap_cb.currentTextChanged.connect(self.activate_update_button)

        form1_layout = QFormLayout()
        form1_layout.addRow(spec_type_label, self.spec_type_cb)
        form1_layout.addRow(self.omega0_label, self.omega0_spin)
        form1_layout.addRow(colormap_label, self.colormap_cb)

        spec_layout.addLayout(form1_layout)
        spec_gb.setLayout(spec_layout)

        #######################
        color_scale_gb = QGroupBox("Scale")
        color_scale_layout = QVBoxLayout()
        self.none_rb = QRadioButton('None: ')
        self.none_rb.setChecked(True)
        self.none_rb.clicked.connect(self.activate_update_button)
        self.hist_rb = QRadioButton('Histogram: ')
        self.hist_rb.clicked.connect(self.activate_update_button)
        self.colorbar_rb = QRadioButton('Color Bar: ')
        self.colorbar_rb.clicked.connect(self.activate_update_button)
    
        #color_scale_layout.addWidget(self.hist_rb)
        color_scale_layout.addWidget(self.colorbar_rb)
        color_scale_layout.addWidget(self.none_rb)
        color_scale_gb.setLayout(color_scale_layout)

        #######################
        detector_gb = QGroupBox("Detector")

        pval_label = QLabel('pval: ')
        self.pval_spin = QDoubleSpinBox()
        self.pval_spin.valueChanged.connect(self.activate_update_button)
        self.pval_spin.setDecimals(4)
        self.pval_spin.setMinimum(0.0001)
        self.pval_spin.setMaximum(1.0)
        self.pval_spin.setValue(0.002)
        self.pval_spin.setSingleStep(0.001)
        self.pval_spin.setMinimumWidth(70)
        self.pval_spin.setMaximumWidth(70)

        fmin_label = QLabel("Freq min: ")
        self.fmin_spin = QDoubleSpinBox()
        self.fmin_spin.valueChanged.connect(self.activate_update_button)
        self.fmin_spin.setMinimum(0.001)
        self.fmin_spin.setMaximum(10000.0)
        self.fmin_spin.setValue(1.0)    
        self.fmin_spin.setMinimumWidth(70)
        self.fmin_spin.setMaximumWidth(70)

        fmax_label = QLabel("Freq max: ")
        self.fmax_spin = QDoubleSpinBox()
        self.fmax_spin.setMinimum(1.0)
        self.fmax_spin.setMaximum(10000.0)
        self.fmax_spin.setValue(10.0)    # this needs to be set when a spectrogram is created
        self.fmax_spin.setMinimumWidth(70)
        self.fmax_spin.setMaximumWidth(70)

        clust_freq_scale_label = QLabel("Cluster Freq. Scaling: ")
        self.clust_freq_scale_spin = QDoubleSpinBox()
        self.clust_freq_scale_spin.valueChanged.connect(self.activate_update_button)
        self.clust_freq_scale_spin.setMinimum(1.0)
        self.clust_freq_scale_spin.setMaximum(10000.0)
        self.clust_freq_scale_spin.setValue(35.0)
        self.clust_freq_scale_spin.setMinimumWidth(70)
        self.clust_freq_scale_spin.setMaximumWidth(70)

        clust_eps_label = QLabel("Clustering EPS: ")
        self.clust_eps_spin = QDoubleSpinBox()
        self.clust_eps_spin.valueChanged.connect(self.activate_update_button)
        self.clust_eps_spin.setMinimum(1.0)
        self.clust_eps_spin.setMaximum(10000.0)
        self.clust_eps_spin.setValue(10.0)
        self.clust_eps_spin.setMinimumWidth(70)
        self.clust_eps_spin.setMaximumWidth(70)

        clust_min_samples_label = QLabel("Clustering Min Samples: ")
        self.clust_min_samples_spin = QSpinBox()
        self.clust_min_samples_spin.valueChanged.connect(self.activate_update_button)
        self.clust_min_samples_spin.setMinimum(1)
        self.clust_min_samples_spin.setMaximum(10000)
        self.clust_min_samples_spin.setValue(40)
        self.clust_min_samples_spin.setMinimumWidth(70)
        self.clust_min_samples_spin.setMaximumWidth(70)

        form2_layout = QFormLayout()
        form2_layout.addRow(pval_label, self.pval_spin)
        form2_layout.addRow(fmin_label, self.fmin_spin)
        form2_layout.addRow(fmax_label, self.fmax_spin)

        form3_layout = QFormLayout()
        form3_layout.addRow(clust_freq_scale_label, self.clust_freq_scale_spin)
        form3_layout.addRow(clust_eps_label, self.clust_eps_spin)
        form3_layout.addRow(clust_min_samples_label, self.clust_min_samples_spin)
        
        det_layout = QHBoxLayout()
        det_layout.addLayout(form2_layout)
        det_layout.addLayout(form3_layout)

        detector_gb.setLayout(det_layout)

        #####SHOW HIDE
        showhide_gb = QGroupBox("Show/Hide")
        self.show_signal_waveform_cb = QCheckBox("Signal Waveform: ")
        self.show_signal_waveform_cb.setChecked(True)
        self.show_signal_waveform_cb.clicked.connect(self.activate_update_button)

        self.show_signal_spectrogram_cb = QCheckBox("Signal Spectrogram: ")
        self.show_signal_spectrogram_cb.setChecked(True)
        self.show_signal_spectrogram_cb.clicked.connect(self.activate_update_button)

        self.show_noise_waveform_cb = QCheckBox("Background Waveform: ")
        self.show_noise_waveform_cb.setChecked(False)
        self.show_noise_waveform_cb.clicked.connect(self.activate_update_button)

        self.show_noise_spectrogram_cb = QCheckBox("Background Spectrogram")
        self.show_noise_spectrogram_cb.setChecked(False)
        self.show_noise_spectrogram_cb.clicked.connect(self.activate_update_button)

        showhide_layout = QVBoxLayout()
        showhide_layout.addWidget(self.show_signal_waveform_cb)
        showhide_layout.addWidget(self.show_signal_spectrogram_cb)
        showhide_layout.addWidget(self.show_noise_waveform_cb)
        showhide_layout.addWidget(self.show_noise_spectrogram_cb)
        showhide_gb.setLayout(showhide_layout)

        self.hide_button = QPushButton("Hide")
        self.hide_button.setMaximumWidth(60)
        self.hide_button.clicked.connect(self.hide)

        h_layout = QHBoxLayout()
        h_layout.addWidget(detector_gb)        
        h_layout.addWidget(spec_gb)
        h_layout.addWidget(color_scale_gb)
        h_layout.addWidget(showhide_gb)
        h_layout.addWidget(self.update_button)
        h_layout.addStretch()
        h_layout.addWidget(self.hide_button)
        h_layout.setContentsMargins(0,0,0,0)
        self.setLayout(h_layout) 

    def get_scale_setting(self):
        if self.none_rb.isChecked():
            return 'none'
        elif self.hist_rb.isChecked():
            return 'hist'
        elif self.colorbar_rb.isChecked():
            return 'cbar'

    def activate_update_button(self):
        self.omega0_label.setVisible(self.spec_type_cb.currentText() == 'CWT')
        self.omega0_spin.setVisible(self.spec_type_cb.currentText() == 'CWT')
            
        self.update_button.setEnabled(True)

    def deactivate_update_button(self):
        self.update_button.setEnabled(False)

    @pyqtSlot()
    def hide(self):
        self.setVisible(False)


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

    def plot_data(self, detections, dt, start_time, t_range, f_range, db):
        #detection data comes in as a list of lists with each element being [t,f,value]
        # first clear any old points
        self.spi.clear()

        # initialize the start time of the axis
        self.getAxis('bottom').set_start_time(start_time)

        #we have to make the spots that will be drawn
        # spots = []
        # for detection in detections:
        #     spots.append({'pos': (detection[0], detection[1]), 'symbol': 's', 'size': 5*dt})

        # self.spi.addPoints(spots)

        #####PARSE THE CLUSTERING  following https://scikit-learn.org/stable/auto_examples/cluster/plot_dbscan.html#sphx-glr-auto-examples-cluster-plot-dbscan-py
        labels = db.labels_
        unique_labels = set(labels)
        core_samples_mask = np.zeros_like(labels, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True


        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)

        core_samples_mask = np.zeros_like(labels, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True

        colors = IPUtils.blue_to_red

        spots = []
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = pg.mkColor(0,0,0)
            class_member_mask = labels == k
           
            xy = detections[class_member_mask & core_samples_mask]
            for data in xy:
                spots.append({'pos': (data[0], data[1]), 'pen': {'color': col}, 'brush': col, 'symbol': 's', 'size':5.6*dt})

            xy = detections[class_member_mask & ~core_samples_mask]
            for data in xy:
                spots.append({'pos': (data[0], data[1]), 'pen': {'color': col}, 'brush': col, 'symbol': 's', 'size':5.6*dt})

        self.spi.addPoints(spots)

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
        if self.spec_type == 'Spectrogram':
            try:
                self.f, self.t, self.Sxx = spectrogram(self.data[1], 
                                        self.fs, 
                                        nfft=self.nfft, 
                                        noverlap=self.noverlap, 
                                        nperseg=self.nperseg)
            except:
                self.signal_runFinished.emit(False)
                return

        elif self.spec_type == 'STFT':
            try:
                self.f, self.t, self.Sxx = stft(self.data[1], 
                                            fs=self.fs, 
                                            nperseg=self.nperseg, 
                                            noverlap=self.noverlap, 
                                            nfft=self.nfft)
                
            except:
                self.signal_runFinished.emit(False)
                return 

        elif self.spec_type == 'CWT':
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
        
                self.f = self.f[1:]
                widths = self.morlet_omega0 / (2 * np.pi * self.f) * self.fs
                self.Sxx = cwt(self.data[1], morlet2, widths, w=self.morlet_omega0)

            except:
                self.signal_runFinished.emit(False)

        self.signal_runFinished.emit(True)


    def get_results(self):
        if self.Sxx is None:
            return None, None, None
        return self.f, self.t, np.abs(self.Sxx)

        