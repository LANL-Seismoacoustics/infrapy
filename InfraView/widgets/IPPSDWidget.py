from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QLabel, QSpinBox, QDoubleSpinBox, QAbstractSpinBox,
                             QGroupBox, QComboBox, QPushButton)
from PyQt5.QtCore import pyqtSignal
from PyQt5 import QtCore


import pyqtgraph as pg

import numpy as np
from scipy import signal

from InfraView.widgets.IPPlotItem import IPPlotItem


class IPPSDWidget(QWidget):

    windows = ['hann', 'hamming', 'boxcar', 'bartlett', 'blackman']

    noiseCurve = None
    signalCurve = None

    currentSignalData = None
    currentNoiseData = None

    blue_pen = pg.mkPen(color=(80, 159, 250), width=1)
    red_pen = pg.mkPen(color=(255, 71, 71), width=1)

    def __init__(self, parent):
        super().__init__()

        self.parent = parent
        self.buildUI()
        # self.updateFrequencyIndicators()
        self.show()

    def buildUI(self):
        self.setMinimumSize(100, 100)
        self.plotLayoutWidget = pg.GraphicsLayoutWidget()
        self.psdPlot = IPPlotItem(mode='PSD')

        self.psdPlot.enableAutoRange(self.psdPlot.xaxis(), enable=True)
        self.psdPlot.enableAutoRange(self.psdPlot.yaxis(), enable=True)
        
        self.psdPlot.getAxis('bottom').setRange(.1, 10)
        
        self.psdPlot.setLabel('bottom', 'f (Hz)')
        self.psdPlot.setLabel('left', 'Power Spectral Density (dB)')
        self.psdPlot.setTitle("...")
        
        self.psdPlot.setLogMode(x=True, y=False)

        initdata = np.array([1])

        self.noiseCurve = self.psdPlot.plot(x=initdata, y=initdata, pen=self.red_pen, name="Noise")
        self.signalCurve = self.psdPlot.plot(x=initdata, y=initdata, pen=self.blue_pen, name="Signal")

        self.plotLayoutWidget.addItem(self.psdPlot, 0, 0)

        label_fft_N = QLabel(self.tr('fft window (N): '))
        label_fft_N.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.fft_N_Spin = QSpinBox()
        self.fft_N_Spin.setMinimum(4)
        self.fft_N_Spin.setMaximum(2**20)
        self.fft_N_Spin.setValue(1024)

        label_fs = QLabel(self.tr('Sampling Freq.: '))
        label_fs.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.fs_Spin = QDoubleSpinBox()
        self.fs_Spin.setMaximum(1000000.0)
        self.fs_Spin.setMinimum(0.0)
        self.fs_Spin.setValue(20.0)
        self.fs_Spin.setReadOnly(True)
        self.fs_Spin.setSuffix(' Hz')
        self.fs_Spin.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.fs_Spin.setEnabled(False)

        label_fft_time = QLabel(self.tr('fft window: '))
        label_fft_time.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.fft_T_Spin = QDoubleSpinBox()
        self.fft_T_Spin.setMaximum(10000.)
        self.fft_T_Spin.setMinimum(0.1)
        self.fft_T_Spin.setValue(1.0)
        self.fft_T_Spin.setSuffix(' s')

        label_window = QLabel(self.tr('Window: '))
        label_window.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.window_cb = QComboBox()
        for window in self.windows:
            self.window_cb.addItem(window)

        label_f1 = QLabel('f1')
        label_f1.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.f1_Spin = QDoubleSpinBox()
        self.f1_Spin.setDecimals(5)
        self.f1_Spin.setMinimum(0.00001)
        self.f1_Spin.setMaximum(10000)
        self.f1_Spin.setSingleStep(0.1)
        self.f1_Spin.setSuffix(' Hz')
        self.f1_Spin.setValue(10**(self.psdPlot.getFreqRegion().getRegion()[0]))

        label_f2 = QLabel('f2')
        label_f2.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.f2_Spin = QDoubleSpinBox()
        self.f2_Spin.setDecimals(5)
        self.f2_Spin.setMinimum(0.00002)
        self.f2_Spin.setMaximum(10000)
        self.f2_Spin.setSingleStep(0.1)
        self.f2_Spin.setSuffix(' Hz')
        self.f2_Spin.setValue(10**(self.psdPlot.getFreqRegion().getRegion()[1]))

        set_filter_button = QPushButton('<- Set filter to ->')
        set_filter_button.clicked.connect(self.setFilterFromPSD)

        f_indicator_layout = QHBoxLayout()
        f_indicator_layout.addWidget(label_f1)
        f_indicator_layout.addWidget(self.f1_Spin)
        f_indicator_layout.addStretch()
        f_indicator_layout.addWidget(set_filter_button)
        f_indicator_layout.addStretch()
        f_indicator_layout.addWidget(label_f2)
        f_indicator_layout.addWidget(self.f2_Spin)

        parametersLayout = QGridLayout()
        parametersLayout.addWidget(label_fft_N, 0, 0)
        parametersLayout.addWidget(self.fft_N_Spin, 0, 1)
        parametersLayout.addWidget(label_fs, 0, 2)
        parametersLayout.addWidget(self.fs_Spin, 0, 3)
        parametersLayout.addWidget(label_fft_time, 1, 0)
        parametersLayout.addWidget(self.fft_T_Spin, 1, 1)
        parametersLayout.addWidget(label_window, 1, 2)
        parametersLayout.addWidget(self.window_cb, 1, 3)

        self.parametersGroup = QGroupBox()
        self.parametersGroup.setLayout(parametersLayout)

        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.plotLayoutWidget)
        mainLayout.addLayout(f_indicator_layout)
        mainLayout.addWidget(self.parametersGroup)

        self.setLayout(mainLayout)

        self.connectSignalsAndSlots()

    def connectSignalsAndSlots(self):
        self.fft_N_Spin.valueChanged.connect(self.updatePSDs)
        self.fft_N_Spin.valueChanged.connect(self.updateFFtT)
        self.fft_T_Spin.valueChanged.connect(self.updateFFtN)

        self.window_cb.currentIndexChanged.connect(self.updatePSDs)

        self.psdPlot.getFreqRegion().sigRegionChanged.connect(self.updateFrequencyIndicators)
        self.f1_Spin.editingFinished.connect(self.updateLinearFrequencyIndicators)
        self.f2_Spin.editingFinished.connect(self.updateLinearFrequencyIndicators)

        # update the beamformingwidget's settings
        # TODO RECONNECT
        # self.psdPlot.getFreqRegion().sigRegionChanged.connect(self.parent.beamformingWidget.bottomSettings.setFreqValues)

    def updateFFtT(self):
        self.fft_T_Spin.setValue(self.fft_N_Spin.value() / self.fs_Spin.value())

    def updateFFtN(self):
        self.fft_N_Spin.setValue(int(self.fft_T_Spin.value() * self.fs_Spin.value()))

    def updatePSDs(self):
        # Use this convenience method to update both PSDs when parameters in this widget are changed
        # Use the indivicual ones below for passing data when the waveform linearregionitems are resized
        # or when a different plot is activated

        self.updateSignalPSD()
        self.updateNoisePSD()

    def updateNoisePSD(self, data=None):
        # if new data is passed, use that, otherwise use what we have
        if data is not None:
            self.currentNoiseData = data.copy()

        if self.currentNoiseData is not None:
            f, pxx = self.calculate_psd(self.currentNoiseData)
            self.noiseCurve.setData(f, 10*np.log10(pxx), pen=self.red_pen)

    def updateSignalPSD(self, data=None):
        # if new data s passed, use that, otherwise use what we have
        if data is not None:
            self.currentSignalData = data.copy()

        if self.currentSignalData is not None:
            f, pxx = self.calculate_psd(self.currentSignalData)
            self.signalCurve.setData(f, 10*np.log10(pxx), pen=self.blue_pen)

    def calculate_psd(self, data):
        if data is not None:
            my_window = self.window_cb.currentText()

            my_nperseg = self.fft_N_Spin.value()
            if my_nperseg > len(data):
                my_nperseg = len(data)

            my_fs = self.fs_Spin.value()

            my_noverlap = int(my_nperseg / 2)

            f, pxx = signal.welch(data, my_fs, my_window, nperseg=my_nperseg, noverlap=my_noverlap)

            return f, pxx

    def set_title(self, title):
        self.psdPlot.setTitle(title)

    def set_fs(self, fs):
        self.fs_Spin.setValue(fs)
        self.updateFFtT()

    @QtCore.pyqtSlot(tuple)
    def updateFrequencyIndicators(self, region):
        r = region.getRegion()
        self.f1_Spin.setValue(10**r[0])
        self.f2_Spin.setValue(10**r[1])

    def updateLinearFrequencyIndicators(self):
        self.psdPlot.getFreqRegion().setRegion([np.log10(self.f1_Spin.value()), np.log10(self.f2_Spin.value())])

    def updateFreqRange(self, rgn):
        self.f1_Spin.setValue(rgn[0])
        self.f2_Spin.setValue(rgn[1])
        self.psdPlot.getFreqRegion().setRegion([np.log10(self.f1_Spin.value()), np.log10(self.f2_Spin.value())])

    def setFilterFromPSD(self):

        filterSettings = self.parent.filterSettingsWidget.get_filter_settings()
        filter_display_settings = self.parent.filterSettingsWidget.get_filter_display_settings()
        filterSettings['type'] = 'Band Pass'
        filterSettings['F_low'] = self.f2_Spin.value()
        filterSettings['F_high'] = self.f1_Spin.value()
        self.parent.filterSettingsWidget.set_filter_settings(filterSettings)

    def clearPlot(self):
        self.noiseCurve.setData([1], [1])
        self.signalCurve.setData([1], [1])
        self.set_title('...')
