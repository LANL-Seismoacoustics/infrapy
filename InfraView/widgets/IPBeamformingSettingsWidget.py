from PyQt5.QtWidgets import (QWidget, QComboBox, QCheckBox, QLabel, QAbstractSpinBox, QDoubleSpinBox, QSpinBox,
                             QGridLayout,
                             QVBoxLayout, QHBoxLayout, QGroupBox,
                             QButtonGroup, QFrame)
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal, QSettings, Qt
import sys


class IPBeamformingSettingsWidget(QWidget):

    def __init__(self, parent):
        super().__init__()
        self.__parent = parent  # reference to the IPBeamformingWidget to which this belongs
        self.buildUI()

    def buildUI(self):

        label_windowLength = QLabel(self.tr('   Window Length: '))
        label_windowLength.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.windowLength_spin = QDoubleSpinBox()
        self.windowLength_spin.setMaximumWidth(60)
        self.windowLength_spin.setSuffix(' s')
        self.windowLength_spin.setMinimum(0.0)
        self.windowLength_spin.setMaximum(1000000)
        self.windowLength_spin.setValue(10.0)

        label_windowStep = QLabel(self.tr('   Window Step: '))
        label_windowStep.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.windowStep_spin = QDoubleSpinBox()
        self.windowStep_spin.setMaximumWidth(60)
        self.windowStep_spin.setSuffix(' s')
        self.windowStep_spin.setMinimum(0.0)
        self.windowStep_spin.setMaximum(1000000)
        self.windowStep_spin.setValue(2.5)

        label_fmin = QLabel(self.tr('   Freq. Min.: '))
        label_fmin.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.fmin_spin = QDoubleSpinBox()
        self.fmin_spin.setMaximumWidth(60)
        self.fmin_spin.setSuffix(' Hz')
        self.fmin_spin.setMinimum(0.0)
        self.fmin_spin.setValue(0.5)    # if changed here, make sure you change in the IPPSDWidget
        self.fmin_spin.setReadOnly(True)
        self.fmin_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        label_fmax = QLabel(self.tr('   Freq. Max.: '))
        label_fmax.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.fmax_spin = QDoubleSpinBox()
        self.fmax_spin.setMaximumWidth(60)
        self.fmax_spin.setSuffix(' Hz')
        self.fmax_spin.setMinimum(0.0)
        self.fmax_spin.setValue(5.0)    # if changed here, make sure you change in the IPPSDWidget
        self.fmax_spin.setReadOnly(True)
        self.fmax_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        label_numSigs = QLabel(self.tr('    Num. of Signals: '))
        label_numSigs.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.numSigs_spin = QSpinBox()
        self.numSigs_spin.setMaximumWidth(40)
        self.numSigs_spin.setMinimum(1)
        self.numSigs_spin.setValue(1)
        self.numSigs_spin.setEnabled(False)

        label_method = QLabel(self.tr('   Method: '))
        label_method.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.method_cb = QComboBox()
        self.method_cb.addItem('bartlett')
        self.method_cb.addItem('gls')
        self.method_cb.addItem('bartlett_covar')
        self.method_cb.addItem('capon')
        self.method_cb.addItem('music')
        self.method_cb.currentTextChanged.connect(self.methodChanged)

        label_noiseStart = QLabel(self.tr('     Noise Range Start: '))
        label_noiseStart.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.noiseStart_spin = QDoubleSpinBox()
        self.noiseStart_spin.setMinimum(0.0)
        self.noiseStart_spin.setMaximum(1000000)
        self.noiseStart_spin.setMaximumWidth(60)
        self.noiseStart_spin.setSuffix(' s')
        self.noiseStart_spin.setReadOnly(True)
        self.noiseStart_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        label_noiseDuration = QLabel(self.tr(' Duration: '))
        label_noiseDuration.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.noiseDuration_spin = QDoubleSpinBox()
        self.noiseDuration_spin.setMaximumWidth(60)
        self.noiseDuration_spin.setMinimum(0.0)
        self.noiseDuration_spin.setMaximum(1000000)
        self.noiseDuration_spin.setSuffix(' s')
        self.noiseDuration_spin.setReadOnly(True)
        self.noiseDuration_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        label_sigStart = QLabel(self.tr('     Signal Range Start: '))
        label_sigStart.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.sigStart_spin = QDoubleSpinBox()
        self.sigStart_spin.setMaximumWidth(60)
        self.sigStart_spin.setMinimum(0.0)
        self.sigStart_spin.setMaximum(1000000)
        self.sigStart_spin.setSuffix(' s')
        self.sigStart_spin.setReadOnly(True)
        self.sigStart_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        label_sigDuration = QLabel(self.tr(' Duration: '))
        label_sigDuration.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
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

        label_subwindow = QLabel(self.tr(' Subwindow Length: '))
        label_subwindow.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.subWinLength_spin = QDoubleSpinBox()
        self.subWinLength_spin.setMaximumWidth(60)
        self.subWinLength_spin.setMinimum(0.0)
        self.subWinLength_spin.setMaximum(1000000)
        self.subWinLength_spin.setValue(self.windowLength_spin.value())
        self.subWinLength_spin.setSuffix(' s')
        self.subWinLength_spin.setEnabled(False)

        label_backaz_resol = QLabel(self.tr(' Back Azimuth Resolution: '))
        self.backaz_resol_spin = QDoubleSpinBox()
        self.backaz_resol_spin.setMinimum(0.1)
        self.backaz_resol_spin.setMaximum(20.)
        self.backaz_resol_spin.setValue(3.0)
        self.backaz_resol_spin.setSuffix(' deg')

        label_tracev_resol = QLabel(self.tr(' Trace Velocity Resolution: '))
        self.tracev_resol_spin = QDoubleSpinBox()
        self.tracev_resol_spin.setMinimum(0.1)
        self.tracev_resol_spin.setMaximum(20.)
        self.tracev_resol_spin.setValue(5.0)
        self.tracev_resol_spin.setSuffix(' m/s')

        label_backaz_startF = QLabel(self.tr(' Back Azimuth Start F: '))
        self.backaz_startF_spin = QDoubleSpinBox()
        self.backaz_startF_spin.setMinimum(-180.0)
        self.backaz_startF_spin.setMaximum(179.0)
        self.backaz_startF_spin.setValue(-180.0)
        self.backaz_startF_spin.setSuffix(' Hz')

        label_backaz_endF = QLabel(self.tr(' Back Azimuth End F: '))
        self.backaz_endF_spin = QDoubleSpinBox()
        self.backaz_endF_spin.setMinimum(-179.0)
        self.backaz_endF_spin.setMaximum(180.0)
        self.backaz_endF_spin.setValue(180.0)
        self.backaz_endF_spin.setSuffix(' Hz')

        gridLayout = QGridLayout()
        gridLayout.addWidget(label_method, 0, 0)
        gridLayout.addWidget(self.method_cb, 0, 1)
        gridLayout.addWidget(label_numSigs, 1, 0)
        gridLayout.addWidget(self.numSigs_spin, 1, 1)

        gridLayout.addWidget(label_windowLength, 0, 2)
        gridLayout.addWidget(self.windowLength_spin, 0, 3)
        gridLayout.addWidget(label_windowStep, 1, 2)
        gridLayout.addWidget(self.windowStep_spin, 1, 3)
        gridLayout.addWidget(label_subwindow, 2, 2)
        gridLayout.addWidget(self.subWinLength_spin, 2, 3)
        gridLayout.addWidget(self.subwindow_cb, 2, 4)

        gridLayout.addWidget(label_fmin, 0, 5)
        gridLayout.addWidget(self.fmin_spin, 0, 6)
        gridLayout.addWidget(label_fmax, 1, 5)
        gridLayout.addWidget(self.fmax_spin, 1, 6)

        gridLayout.addWidget(label_noiseStart, 0, 7)
        gridLayout.addWidget(self.noiseStart_spin, 0, 8)
        gridLayout.addWidget(label_noiseDuration, 0, 9)
        gridLayout.addWidget(self.noiseDuration_spin, 0, 10)

        gridLayout.addWidget(label_sigStart, 1, 7)
        gridLayout.addWidget(self.sigStart_spin, 1, 8)
        gridLayout.addWidget(label_sigDuration, 1, 9)
        gridLayout.addWidget(self.sigDuration_spin, 1, 10)

        gridLayout.addWidget(label_backaz_resol, 0, 11)
        gridLayout.addWidget(self.backaz_resol_spin, 0, 12)
        gridLayout.addWidget(label_backaz_startF, 1, 11)
        gridLayout.addWidget(self.backaz_startF_spin, 1, 12)
        gridLayout.addWidget(label_backaz_endF, 2, 11)
        gridLayout.addWidget(self.backaz_endF_spin, 2, 12)
        gridLayout.addWidget(label_tracev_resol, 3, 11)
        gridLayout.addWidget(self.tracev_resol_spin, 3, 12)

        horizLayout = QHBoxLayout()
        horizLayout.addLayout(gridLayout)
        horizLayout.addStretch()

        vertLayout = QVBoxLayout()
        vertLayout.addLayout(horizLayout)
        vertLayout.addStretch()

        self.setLayout(vertLayout)

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

    def getTraceVelResolution(self):
        return self.tracev_resol_spin.value()

    def getBackAzResolution(self):
        return self.backaz_resol_spin.value()

    def getBackAzFreqRange(self):
        return (self.backaz_startF_spin.value(), self.backaz_endF_spin.value())

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
