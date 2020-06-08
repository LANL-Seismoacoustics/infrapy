from PyQt5.QtWidgets import (QWidget, QComboBox, QCheckBox, QLabel, QAbstractSpinBox, QDoubleSpinBox, QSpinBox,
                             QFormLayout, QGridLayout,
                             QVBoxLayout, QHBoxLayout, QGroupBox,
                             QButtonGroup, QFrame)
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal, pyqtSlot, QSettings, Qt
import sys


class IPDetectorSettingsWidget(QWidget):

    auto_threshold_level = None

    def __init__(self, parent):
        super().__init__()
        self.__parent = parent  # reference to the IPBeamformingWidget to which this belongs
        self.buildUI()

    def buildUI(self):
        self.det_type_group = QButtonGroup(self)
        self.det_type_group.setExclusive(True)
        self.det_type_group.buttonClicked.connect(self.update_widget)

        self.auto_checkbox = QCheckBox()
        self.auto_checkbox.setChecked(True)
        self.auto_checkbox.setToolTip('If auto is selected, the program will estimate the background f-stat from the noise region selected on in the waveform window.')

        self.manual_checkbox = QCheckBox()
        self.manual_value = QDoubleSpinBox()
        self.manual_value.setValue(1.0)
        self.manual_value.setEnabled(False)

        main_layout = QGridLayout()

        form_layout = QFormLayout()

        self.det_type_group.addButton(self.auto_checkbox)
        self.det_type_group.addButton(self.manual_checkbox)

        form_layout.addRow("Automatically calculate threshold: ", self.auto_checkbox)
        form_layout.addRow("Manually calculate threshold: ", self.manual_checkbox)
        form_layout.addRow("Manual threshold level: ", self.manual_value)

        main_layout.addLayout(form_layout,0,0)

        self.setLayout(main_layout)

    @pyqtSlot()
    def update_widget(self):
        self.manual_value.setEnabled(self.manual_checkbox.isChecked())

    def is_auto_threshold(self):
        return self.auto_checkbox.isChecked()

    @pyqtSlot(float)
    def set_auto_threshold_level(self, level):
        self.auto_threshold_level = level

    def get_auto_threshold_level(self):
        return self.auto_threshold_level

    def get_manual_threshold_level(self):
        return self.manual_value.value()


        