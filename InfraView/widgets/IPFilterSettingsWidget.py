from PyQt5.QtWidgets import (QWidget, QCheckBox, QComboBox, QLabel, QDoubleSpinBox,
                             QSpinBox, QGridLayout, QGroupBox, QPushButton,
                             QVBoxLayout, QGroupBox)

from PyQt5.QtCore import pyqtSignal, QSettings, Qt


class IPFilterSettingsWidget(QWidget):

    filterChanged = pyqtSignal()
    sig_filter_changed = pyqtSignal(dict)
    sig_filter_display_changed = pyqtSignal(dict)

    _filter_display_settings = {}  # Holder of the current filter display settings

    _filter_display_settings_default = {'apply': False, 'showUnfiltered': False}

    _filter_settings = {}  # Holder of the current filter settings

    _filter_settings_default = {'type': 'Band Pass', 'F_low': 5.0, 'F_high': .5,
                                'order': 4, 'zphase': False}

    def __init__(self, parent=None):

        super().__init__(parent)

        self._filter_settings = self._filter_settings_default.copy()
        self._filter_display_settings = self._filter_display_settings_default.copy()
        self.__buildUI__()
        self.show()

    def __buildUI__(self):

        self.applyFilter_checkbox = QCheckBox('Apply Filter?')
        self.applyFilter_checkbox.setChecked(self._filter_display_settings['apply'])
        self.applyFilter_checkbox.stateChanged.connect(self.apply_filter)

        self.showUnfiltered = QCheckBox('Show Unfiltered?')
        self.showUnfiltered.setChecked(self._filter_display_settings['showUnfiltered'])
        self.showUnfiltered.stateChanged.connect(self.onActivated_showUnfiltered)

        self.cb_filter_type = QComboBox()
        self.cb_filter_type.addItem('Low Pass')
        self.cb_filter_type.addItem('High Pass')
        self.cb_filter_type.addItem('Band Pass')

        cb_idx = self.cb_filter_type.findText(self._filter_settings['type'])
        self.cb_filter_type.setCurrentIndex(cb_idx)
        self.cb_filter_type.currentIndexChanged[str].connect(self.onActivated_cb)

        self.label_lowpassFreq = QLabel(self.tr('Low Pass F: '))
        self.lowpassSpin = QDoubleSpinBox()
        self.lowpassSpin.setDecimals(3)
        self.lowpassSpin.setValue(self._filter_settings['F_low'])
        self.lowpassSpin.setMaximum(100000)
        self.lowpassSpin.setMinimum(.001)
        self.lowpassSpin.setSingleStep(.01)
        self.lowpassSpin.setSuffix(' Hz')
        self.lowpassSpin.valueChanged.connect(self.onActivated_lpSpin)

        self.label_highpassFreq = QLabel(self.tr('High Pass F: '))
        self.highpassSpin = QDoubleSpinBox()
        self.highpassSpin.setDecimals(3)
        self.highpassSpin.setValue(self._filter_settings['F_high'])
        self.highpassSpin.setMaximum(100000)
        self.highpassSpin.setMinimum(.001)
        self.highpassSpin.setSingleStep(.01)
        self.highpassSpin.setSuffix(' Hz')
        self.highpassSpin.valueChanged.connect(self.onActivated_hpSpin)

        self.label_order = QLabel(self.tr('Order: '))
        self.orderSpin = QSpinBox()
        self.orderSpin.setMinimum(1)
        self.orderSpin.setValue(self._filter_settings['order'])

        self.label_zeroPhase = QLabel(self.tr('Zero Phase'))
        self.zeroPhase = QCheckBox()
        self.zeroPhase.setChecked(self._filter_settings['zphase'])

        self.update_Button = QPushButton('Update')
        self.update_Button.setMaximumWidth(200)
        self.update_Button.clicked.connect(self.update_clicked)

        layout = QGridLayout(self)
        # layout.setMargin(20)
        layout.addWidget(QLabel(self.tr('Type: ')), 0, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.cb_filter_type, 0, 1)
        layout.addWidget(self.label_highpassFreq, 1, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.highpassSpin, 1, 1)
        layout.addWidget(self.label_lowpassFreq, 2, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.lowpassSpin, 2, 1)
        layout.addWidget(self.label_order, 3, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.orderSpin, 3, 1)
        layout.addWidget(self.label_zeroPhase, 4, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.zeroPhase, 4, 1)
        layout.addWidget(self.update_Button, 5, 1)

        self.settingsBox = QGroupBox("Settings: ")
        self.settingsBox.setLayout(layout)

        # layout.addWidget(self.showSpect_Button, 6, 1)  # Slow as hell

        qvbox = QVBoxLayout()
        qvbox.addWidget(self.applyFilter_checkbox)
        qvbox.addWidget(self.showUnfiltered)
        qvbox.addWidget(self.settingsBox)
        qvbox.addStretch()

        # default setting is to disable all inputs except for the applyFilter checkbox
        self.disableAll()

        self.setLayout(qvbox)

    def onActivated_cb(self, text):
        self._filter_settings['type'] = text

        if text == 'Low Pass' or text == 'Low Pass Cheby2' or text == 'Low Pass Fir':
            # disable the highpass freq. spin
            self.label_highpassFreq.setEnabled(False)
            self.highpassSpin.setEnabled(False)
            # make sure the lowpass freq. spin is enabled
            self.label_lowpassFreq.setEnabled(True)
            self.lowpassSpin.setEnabled(True)
        elif text == 'High Pass':
            # enable the high pass freq. spin
            self.label_highpassFreq.setEnabled(True)
            self.highpassSpin.setEnabled(True)
            # make sure the lowpass freq. spin is disabled
            self.label_lowpassFreq.setEnabled(False)
            self.lowpassSpin.setEnabled(False)
        elif text == 'Band Pass':
            # enable bothsts, sts_filtered, normalize=False
            self.label_highpassFreq.setEnabled(True)
            self.highpassSpin.setEnabled(True)
            self.label_lowpassFreq.setEnabled(True)
            self.lowpassSpin.setEnabled(True)
        else:
            # something weird is happening, just bail
            return

        # since something changed, we need to make sure the update button is enabled    
        self.update_Button.setEnabled(True)

    def onActivated_showUnfiltered(self):
        self._filter_display_settings['showUnfiltered'] = self.showUnfiltered.isChecked()
        self.sig_filter_display_changed.emit(self._filter_display_settings)

    def update_clicked(self):
        self.sig_filter_changed.emit(self._filter_settings)

    def onActivated_zeroPhase(self, int):
        self._filter_settings['zphase'] = self.zeroPhase.isChecked()

    def onActivated_lpSpin(self, float):
        self._filter_settings['F_low'] = self.lowpassSpin.value()

    def onActivated_hpSpin(self, float):
        self._filter_settings['F_high'] = self.highpassSpin.value()

    def enableAll(self):
        # This enables all of the appropriate inputs
        self.cb_filter_type.setEnabled(True)
        self.label_order.setEnabled(True)
        self.orderSpin.setEnabled(True)
        self.onActivated_cb(self._filter_settings['type'])
        self.label_zeroPhase.setEnabled(True)
        self.zeroPhase.setEnabled(True)
        self.showUnfiltered.setEnabled(True)
        self.update_Button.setEnabled(True)

    def disableAll(self):
        # This disables all the inputs EXCEPT for the applyFilter checkbox
        self.cb_filter_type.setEnabled(False)
        self.lowpassSpin.setEnabled(False)
        self.highpassSpin.setEnabled(False)
        self.label_highpassFreq.setEnabled(False)
        self.label_lowpassFreq.setEnabled(False)
        self.label_order.setEnabled(False)
        self.orderSpin.setEnabled(False)
        self.label_zeroPhase.setEnabled(False)
        self.zeroPhase.setEnabled(False)
        self.showUnfiltered.setEnabled(False)
        self.update_Button.setEnabled(False)

    def apply_filter(self, state):

        if state == 2:
            self._filter_display_settings['apply'] = True
            self.enableAll()
        else:
            self._filter_display_settings['apply'] = False

        self.sig_filter_changed.emit(self._filter_settings)
        self.sig_filter_display_changed.emit(self._filter_display_settings)

    def reset_filter_settings(self):
        self._filter_settings = self._filter_settings_default.copy()
        self._filter_display_settings = self._filter_display_settings_default.copy()

        self.updateWidget()
        self.disableAll()

    def get_filter_settings(self):
        return self._filter_settings

    def get_filter_display_settings(self):
        return self._filter_display_settings

    def set_filter_settings(self, settings):
        self._filter_settings = settings
        self.sig_filter_changed.emit(settings)
        self.update_widget()

    def update_widget(self):
        # when filter settings are changed programatically, update the widget to show current settings

        self.applyFilter_checkbox.setChecked(self._filter_display_settings['apply'])
        self.showUnfiltered.setChecked(self._filter_display_settings['showUnfiltered'])

        cb_idx = self.cb_filter_type.findText(self._filter_settings['type'])
        self.cb_filter_type.setCurrentIndex(cb_idx)

        self.lowpassSpin.setValue(self._filter_settings['F_low'])
        self.highpassSpin.setValue(self._filter_settings['F_high'])
        self.orderSpin.setValue(self._filter_settings['order'])
        self.zeroPhase.setChecked(self._filter_settings['zphase'])

    def save_current_filter(self):
        newSettings = QSettings('LANL', 'IPView')
        newFilterName = "Billy"
        newFilterKey = newFilterName + "/FilterType"
        newSettings.setValue("newFilterKey", "Butterworth")

    def refresh_filter_entries(self):
        self.cb_filter_type.setCurrentText(self._filter_settings['type'])
        self.lowpassSpin.setValue(self._filter_settings['F_low'])
        self.highpassSpin.setValue(self._filter_settings['F_high'])
        self.orderSpin.setValue(self._filter_settings['order'])
        self.zeroPhase.setChecked(self._filter_settings['zphase'])
