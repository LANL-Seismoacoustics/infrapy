from PyQt5.QtWidgets import (QWidget, QComboBox, QCheckBox, QLabel, QDoubleSpinBox, QSpinBox,
                             QHBoxLayout, QFormLayout, QFrame)
from PyQt5 import QtCore


class IPBeamformingSettingsWidget(QWidget):

    def __init__(self, parent):
        super().__init__()
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

        self.fmin_label = QLabel("0.5 Hz") # if changed here, make sure you change in the IPPSDWidget
        self.fmin_label.setMinimumWidth(90)
        self.fmax_label = QLabel("5.0 Hz")
        self.fmax_label.setMinimumWidth(90)   

        self.noiseStart_label = QLabel()
        self.noiseDuration_label = QLabel()

        self.sigStart_label = QLabel()
        self.sigDuration_label = QLabel()

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

        # removed col3

        formlayout_col4 = QFormLayout()
        formlayout_col4.addRow("Freq Min: ", self.fmin_label)
        formlayout_col4.addRow("Noise Range Start: ", self.noiseStart_label)
        formlayout_col4.addRow("Signal Range Start: ", self.sigStart_label)

        formlayout_col5 = QFormLayout()
        formlayout_col5.addRow("Freq Max: ", self.fmax_label)
        formlayout_col5.addRow("Duration: ", self.noiseDuration_label)
        formlayout_col5.addRow("Duration: ", self.sigDuration_label)

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
        self.fmin_label.setText("{:.5f} Hz".format(min))

    def setFmax(self, max):
        self.fmax_label.setText("{:.5f} Hz".format(max))

    def getNoiseRange(self):
        return (float(self.noiseStart_label.text()), float(self.noiseStart_label.text()) + float(self.noiseDuration_label.text()))

    def getSignalRange(self):
        return (float(self.sigStart_label.text()), float(self.sigStart_label.text()) + float(self.sigDuration_label.text()))

    def getFreqRange(self):
        # we need to remove the " Hz" from the label strings...
        min = float(self.fmin_label.text()[:-3])
        max = float(self.fmax_label.text()[:-3])
        return (min, max)

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
        self.noiseStart_label.setText("{:.2f}".format(values[0]))
        self.noiseDuration_label.setText("{:.2f}".format(values[1] - values[0]))

    @QtCore.pyqtSlot(tuple)
    def setSignalValues(self, values):   # values is a tuple containing (start, stop)
        self.sigStart_label.setText("{:.2f}".format(values[0]))
        self.sigDuration_label.setText("{:.2f}".format(values[1] - values[0]))

    @QtCore.pyqtSlot(tuple)
    def setFreqValues(self, IPLinearRegionItem):
        values = IPLinearRegionItem.getRegion()
        self.fmin_label.setText("{:.5f} Hz".format(10**values[0]))
        self.fmax_label.setText("{:.5f} Hz".format(10**values[1]))

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
