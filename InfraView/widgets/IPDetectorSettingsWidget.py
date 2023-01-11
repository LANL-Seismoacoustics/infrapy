from PyQt5.QtWidgets import (QWidget, QCheckBox, QDoubleSpinBox, QSpinBox,
                             QFormLayout, QHBoxLayout, QButtonGroup)

from PyQt5.QtCore import pyqtSlot


class IPDetectorSettingsWidget(QWidget):

    auto_threshold_level = None

    def __init__(self, parent):
        super().__init__(parent)
        self.buildUI()

    def buildUI(self):
        self.det_type_group = QButtonGroup(self)
        self.det_type_group.setExclusive(True)
        self.det_type_group.buttonClicked.connect(self.update_widget)

        self.auto_checkbox = QCheckBox()
        self.auto_checkbox.setChecked(True)
        self.auto_checkbox.setToolTip('If auto is selected, the program will estimate the background f-stat from the noise region selected on in the waveform window.')

        self.back_az_limit = QDoubleSpinBox()
        self.back_az_limit.setSuffix(' deg')
        self.back_az_limit.setValue(10.0)
        self.back_az_limit.setMinimum(1.0)
        self.back_az_limit.setMaximum(360.0)
        self.back_az_limit.setToolTip("This is the range in degrees of variation in back azimuths used to determine if a value is part of a detection.")

        self.min_peak_width = QSpinBox()
        self.min_peak_width.setSuffix(' points')
        self.min_peak_width.setValue(5)
        self.min_peak_width.setMinimum(1)
        self.min_peak_width.setToolTip("This determines the number of points in the f-statistic plot required for a valid detection.")

        self.manual_checkbox = QCheckBox()
        self.manual_value = QDoubleSpinBox()
        self.manual_value.setValue(1.0)
        self.manual_value.setEnabled(False)

        main_layout = QHBoxLayout()

        form_layout_col1 = QFormLayout()
        form_layout_col2 = QFormLayout()

        self.det_type_group.addButton(self.auto_checkbox)
        self.det_type_group.addButton(self.manual_checkbox)

        form_layout_col1.addRow("Automatically calculate threshold: ", self.auto_checkbox)
        form_layout_col1.addRow("Manually calculate threshold: ", self.manual_checkbox)
        form_layout_col1.addRow("Manual threshold level: ", self.manual_value)
        
        form_layout_col2.addRow("Back azimuth limit: ", self.back_az_limit)
        form_layout_col2.addRow("Minimum peak width: ", self.min_peak_width)
        
        main_layout.addLayout(form_layout_col1)
        main_layout.addLayout(form_layout_col2)
        main_layout.addStretch()
        main_layout.insertSpacing(2, 20)

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


        