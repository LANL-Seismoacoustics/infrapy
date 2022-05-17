from PyQt5.QtWidgets import (QWidget, QComboBox, QCheckBox, QAbstractSpinBox, QDoubleSpinBox, QSpinBox,
                             QHBoxLayout, QFormLayout, QFrame)
from PyQt5 import QtCore


class IPBeamformingSettingsWidget(QWidget):

    def __init__(self, parent):
        super().__init__()
        self.__parent = parent  # reference to the IPBeamformingWidget to which this belongs
        self.buildUI()

    def buildUI(self):

        self.windowLength_spin = QDoubleSpinBox()
        self.windowLength_spin.setMaximumWidth(60)
        self.windowLength_spin.setSuffix(' s')
        self.windowLength_spin.setMinimum(0.0)
        self.windowLength_spin.setMaximum(1000000)
        self.windowLength_spin.setValue(10.0)

        self.windowStep_spin = QDoubleSpinBox()
        self.windowStep_spin.setMaximumWidth(60)
        self.windowStep_spin.setSuffix(' s')
        self.windowStep_spin.setMinimum(0.0)
        self.windowStep_spin.setMaximum(1000000)
        self.windowStep_spin.setValue(2.5)

        self.fmin_spin = QDoubleSpinBox()
        self.fmin_spin.setMaximumWidth(60)
        self.fmin_spin.setSuffix(' Hz')
        self.fmin_spin.setMinimum(0.0)
        self.fmin_spin.setMaximum(10000)
        self.fmin_spin.setValue(0.5)    # if changed here, make sure you change in the IPPSDWidget
        self.fmin_spin.setReadOnly(True)
        self.fmin_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.fmax_spin = QDoubleSpinBox()
        self.fmax_spin.setMaximumWidth(60)
        self.fmax_spin.setSuffix(' Hz')
        self.fmax_spin.setMinimum(0.0)
        self.fmax_spin.setMaximum(10000)
        self.fmax_spin.setValue(5.0)    # if changed here, make sure you change in the IPPSDWidget
        self.fmax_spin.setReadOnly(True)
        self.fmax_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.numSigs_spin = QSpinBox()
        self.numSigs_spin.setMaximumWidth(40)
        self.numSigs_spin.setMinimum(1)
        self.numSigs_spin.setValue(1)
        self.numSigs_spin.setEnabled(False)

        self.method_cb = QComboBox()
        self.method_cb.addItem('bartlett')
        self.method_cb.addItem('gls')
        self.method_cb.addItem('bartlett_covar')
        self.method_cb.addItem('capon')
        self.method_cb.addItem('music')
        self.method_cb.currentTextChanged.connect(self.methodChanged)

        self.noiseStart_spin = QDoubleSpinBox()
        self.noiseStart_spin.setMinimum(0.0)
        self.noiseStart_spin.setMaximum(1000000)
        self.noiseStart_spin.setMaximumWidth(60)
        self.noiseStart_spin.setSuffix(' s')
        self.noiseStart_spin.setReadOnly(True)
        self.noiseStart_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.noiseDuration_spin = QDoubleSpinBox()
        self.noiseDuration_spin.setMaximumWidth(60)
        self.noiseDuration_spin.setMinimum(0.0)
        self.noiseDuration_spin.setMaximum(1000000)
        self.noiseDuration_spin.setSuffix(' s')
        self.noiseDuration_spin.setReadOnly(True)
        self.noiseDuration_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.sigStart_spin = QDoubleSpinBox()
        self.sigStart_spin.setMaximumWidth(60)
        self.sigStart_spin.setMinimum(0.0)
        self.sigStart_spin.setMaximum(1000000)
        self.sigStart_spin.setSuffix(' s')
        self.sigStart_spin.setReadOnly(True)
        self.sigStart_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.sigDuration_spin = QDoubleSpinBox()
        self.sigDuration_spin.setMaximumWidth(60)
        self.sigDuration_spin.setMinimum(0.0)
        self.sigDuration_spin.setMaximum(1000000)
        self.sigDuration_spin.setSuffix(' s')
        self.sigDuration_spin.setReadOnly(True)
        self.sigDuration_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.subwindow_cb = QCheckBox()
        self.subwindow_cb.setEnabled(False)
        self.subwindow_cb.stateChanged.connect(self.updateSubwindow)

        self.subWinLength_spin = QDoubleSpinBox()
        self.subWinLength_spin.setMaximumWidth(60)
        self.subWinLength_spin.setMinimum(0.0)
        self.subWinLength_spin.setMaximum(1000000)
        self.subWinLength_spin.setValue(self.windowLength_spin.value())
        self.subWinLength_spin.setSuffix(' s')
        self.subWinLength_spin.setEnabled(False)

        self.backaz_resol_spin = QDoubleSpinBox()
        self.backaz_resol_spin.setMinimum(0.1)
        self.backaz_resol_spin.setMaximum(20.)
        self.backaz_resol_spin.setValue(3.0)
        self.backaz_resol_spin.setSuffix(' deg')

        self.tracev_resol_spin = QDoubleSpinBox()
        self.tracev_resol_spin.setMinimum(0.1)
        self.tracev_resol_spin.setMaximum(500.)
        self.tracev_resol_spin.setValue(5.0)
        self.tracev_resol_spin.setSuffix(' m/s')

        self.backaz_start_spin = QDoubleSpinBox()
        self.backaz_start_spin.setMinimum(-180.0)
        self.backaz_start_spin.setMaximum(179.0)
        self.backaz_start_spin.setValue(-180.0)
        self.backaz_start_spin.setSuffix(' deg')
        self.backaz_start_spin.editingFinished.connect(self.checkBackAzRange)

        self.backaz_end_spin = QDoubleSpinBox()
        self.backaz_end_spin.setMinimum(-179.0)
        self.backaz_end_spin.setMaximum(180.0)
        self.backaz_end_spin.setValue(180.0)
        self.backaz_end_spin.setSuffix(' deg')
        self.backaz_end_spin.editingFinished.connect(self.checkBackAzRange)

        self.tracev_min_spin = QDoubleSpinBox()
        self.tracev_min_spin.setMinimum(1)
        self.tracev_min_spin.setMaximum(20000.)
        self.tracev_min_spin.setValue(300.0)
        self.tracev_min_spin.setSuffix(' m/s')
        self.tracev_min_spin.editingFinished.connect(self.checkTraceVRange)
        
        self.tracev_max_spin = QDoubleSpinBox()
        self.tracev_max_spin.setMinimum(2)
        self.tracev_max_spin.setMaximum(20000.)
        self.tracev_max_spin.setValue(750.0)
        self.tracev_max_spin.setSuffix(' m/s')
        self.tracev_max_spin.editingFinished.connect(self.checkTraceVRange)


        formlayout_col1 = QFormLayout()
        formlayout_col1.addRow("Method: ", self.method_cb)
        formlayout_col1.addRow("Num. of Signals: ", self.numSigs_spin)

        formlayout_col2 = QFormLayout()
        formlayout_col2.addRow("Window Length: ", self.windowLength_spin)
        formlayout_col2.addRow("Window Step: ", self.windowStep_spin)
        sub_win_layout = QHBoxLayout()
        sub_win_layout.addWidget(self.subWinLength_spin)
        sub_win_layout.addWidget(self.subwindow_cb)
        formlayout_col2.addRow("Subwindow Length: ", sub_win_layout)

        formlayout_col3 = QFormLayout()
        formlayout_col3.addRow("Freq Min: ", self.fmin_spin)
        formlayout_col3.addRow("Freq Max: ", self.fmax_spin)

        formlayout_col4 = QFormLayout()
        formlayout_col4.addRow("Noise Range Start: ", self.noiseStart_spin)
        formlayout_col4.addRow("Signal Range Start: ", self.sigStart_spin)

        formlayout_col5 = QFormLayout()
        formlayout_col5.addRow("Duration: ", self.noiseDuration_spin)
        formlayout_col5.addRow("Duration: ", self.sigDuration_spin)

        formlayout_col6 = QFormLayout()
        formlayout_col6.addRow("Back Azimuth Resolution: ", self.backaz_resol_spin)
        formlayout_col6.addRow("Back Azimuth Start Angle: ", self.backaz_start_spin)
        formlayout_col6.addRow("Back Azimuth End Angle: ", self.backaz_end_spin)

        formlayout_col7 = QFormLayout()
        formlayout_col7.addRow("Trace Vel. Resolution: ", self.tracev_resol_spin)
        formlayout_col7.addRow("Trace Vel Min: ", self.tracev_min_spin)
        formlayout_col7.addRow("Trace Vel Max: ", self.tracev_max_spin)

        horizLayout = QHBoxLayout()
        horizLayout.addLayout(formlayout_col1)
        horizLayout.addLayout(formlayout_col2)
        horizLayout.addLayout(formlayout_col3)
        horizLayout.addLayout(formlayout_col4)
        horizLayout.addLayout(formlayout_col5)
        horizLayout.addLayout(formlayout_col6)
        horizLayout.addLayout(formlayout_col7)
        horizLayout.addStretch()

        self.setLayout(horizLayout)

    def HLine(self):
        hl = QFrame()
        hl.setFrameShape(QFrame.HLine)
        hl.setFrameShadow(QFrame.Sunken)
        return hl

    def VLine(self):
        vl = QFrame()
        vl.setFrameShape(QFrame.VLine)
        vl.setFrameShadow(QFrame.Sunken)
        return vl

    def setFmin(self, min):
        self.fmin_spin.setValue(min)

    def setFmax(self, max):
        self.fmax_spin.setValue(max)

    def getNoiseRange(self):
        return (self.noiseStart_spin.value(), self.noiseStart_spin.value() + self.noiseDuration_spin.value())

    def getSignalRange(self):
        return (self.sigStart_spin.value(), self.sigStart_spin.value() + self.sigDuration_spin.value())

    def getFreqRange(self):
        return (self.fmin_spin.value(), self.fmax_spin.value())

    def getWinLength(self):
        return self.windowLength_spin.value()

    def getSubWinLength(self):
        if self.subWinLength_spin.isEnabled():
            return self.subWinLength_spin.value()
        else:
            return None

    def updateSubwindow(self, state):
        self.subWinLength_spin.setEnabled(state)

    def getNumSigs(self):
        if self.numSigs_spin.isEnabled():
            return self.numSigs_spin.value()
        else:
            return 1

    def getWinStep(self):
        return self.windowStep_spin.value()

    def getMethod(self):
        return self.method_cb.currentText()


    def getBackAzResolution(self):
        return self.backaz_resol_spin.value()

    def getBackAzRange(self):
        return (self.backaz_start_spin.value(), self.backaz_end_spin.value())

    @QtCore.pyqtSlot()
    def checkBackAzRange(self):
        start,stop = self.getBackAzRange()
        if start >= stop:
            self.backaz_start_spin.setStyleSheet("color: rgb(200,0,0); ")
            self.backaz_end_spin.setStyleSheet("color: rgb(200,0,0); ")
        else:
            self.backaz_start_spin.setStyleSheet("color: rgb(0, 0, 0);")
            self.backaz_end_spin.setStyleSheet("color: rgb(0, 0, 0);")

    def getTraceVelResolution(self):
        return self.tracev_resol_spin.value()

    def getTraceVRange(self):
        return (self.tracev_min_spin.value(), self.tracev_max_spin.value())

    @QtCore.pyqtSlot()
    def checkTraceVRange(self):
        min, max = self.getTraceVRange()
        if min >= max:
            self.tracev_min_spin.setStyleSheet("color: rgb(200,0,0);")
            self.tracev_max_spin.setStyleSheet("color: rgb(200,0,0);")
        else:
            self.tracev_min_spin.setStyleSheet("color: rgb(0,0,0);")
            self.tracev_max_spin.setStyleSheet("color: rgb(0,0,0);")

    @QtCore.pyqtSlot(tuple)
    def setNoiseValues(self, values):   # values is a tuple containing (start, stop)
        self.noiseStart_spin.setValue(values[0])
        self.noiseDuration_spin.setValue(values[1] - values[0])

    @QtCore.pyqtSlot(tuple)
    def setSignalValues(self, values):   # values is a tuple containing (start, stop)
        self.sigStart_spin.setValue(values[0])
        self.sigDuration_spin.setValue(values[1] - values[0])

    @QtCore.pyqtSlot(tuple)
    def setFreqValues(self, IPLinearRegionItem):
        values = IPLinearRegionItem.getRegion()
        self.fmin_spin.setValue(10**values[0])
        self.fmax_spin.setValue(10**values[1])

    @QtCore.pyqtSlot(str)
    def methodChanged(self, newMethod):
        if newMethod == 'music':
            self.numSigs_spin.setEnabled(True)
        else:
            self.numSigs_spin.setValue(1)
            self.numSigs_spin.setEnabled(False)

        if newMethod == 'music' or newMethod == 'capon' or newMethod == 'bartlett_covar':
            self.subwindow_cb.setEnabled(True)
            self.subWinLength_spin.setEnabled(self.subwindow_cb.isChecked())
        else:
            self.subwindow_cb.setChecked(False)
            self.subwindow_cb.setEnabled(False)
            self.subWinLength_spin.setEnabled(False)
