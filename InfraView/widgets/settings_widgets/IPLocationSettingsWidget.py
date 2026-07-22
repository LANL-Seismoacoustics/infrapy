from PyQt5.QtWidgets import (QCheckBox, QHBoxLayout, QVBoxLayout, QFormLayout, QColorDialog, QSpinBox,
                             QGroupBox, QComboBox, QLabel, QPushButton, QFileDialog, QDoubleSpinBox)

from PyQt5.QtCore import pyqtSignal, pyqtSlot, QSettings
from PyQt5.QtGui import QColor

import numpy as np

from InfraView.widgets import IPBaseWidgets

class IPLocationSettingsWidget(IPBaseWidgets.IPSettingsWidget):

    signal_colors_changed = pyqtSignal()
    signal_background_changed = pyqtSignal()
    signal_offline_directory_changed = pyqtSignal()
    signal_map_settings_changed = pyqtSignal()
    
    ocean_color = QColor(0, 107, 166)
    land_color = QColor(222, 222, 222)

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent

        self.buildUI()

    def buildUI(self):
        ###   feature settings   ###
        features_gb = IPBaseWidgets.IPSettingsGroupBox("Features")
        self.borders_checkbox = QCheckBox('Countries  ')
        self.states_checkbox = QCheckBox('States and Provinces  ')
        self.lakes_checkbox = QCheckBox('Lakes  ')
        self.rivers_checkbox = QCheckBox('Rivers  ')
        self.coast_checkbox = QCheckBox('Coastline  ')

        features_layout = QVBoxLayout()
        features_layout.addWidget(self.borders_checkbox)
        features_layout.addWidget(self.states_checkbox)
        features_layout.addWidget(self.lakes_checkbox)
        features_layout.addWidget(self.rivers_checkbox)
        features_layout.addWidget(self.coast_checkbox)
        features_layout.addStretch()
        
        features_gb.setLayout(features_layout)

        ###   color settings   ###
        ocean_color_label = QLabel("Oceans: ")
        self.ocean_color_button = IPBaseWidgets.IPColorButton(self.ocean_color)
        land_color_label = QLabel("Land: ")
        self.land_color_button = IPBaseWidgets.IPColorButton(self.land_color)

        colors_layout1 = QFormLayout()
        colors_layout1.addRow(ocean_color_label, self.ocean_color_button)
        colors_layout2 = QFormLayout()
        colors_layout2.addRow(land_color_label, self.land_color_button)
        colors_layout = QHBoxLayout()
        colors_layout.addLayout(colors_layout1)
        colors_layout.addLayout(colors_layout2)

        ###   grid settings
        self.show_grid_checkbox = QCheckBox("Show Grid Lines ")

        ###   resolution settings   ###
        label_resolution = QLabel(self.tr('Resolution'))
        self.resolution_cb = QComboBox()
        self.resolution_cb.addItem('50m')
        self.resolution_cb.addItem('110m')
        self.resolution_cb.setCurrentIndex(1)

        resolution_layout = QHBoxLayout()
        resolution_layout.addWidget(label_resolution)
        resolution_layout.addWidget(self.resolution_cb)

        ### background image ##
        self.backgroud_image_checkbox = QCheckBox('Use background image  ')

        ###   offline maps settings   ###
        self.offline_checkbox = QCheckBox('Use offline maps  ')
        self.offline_directory_label = QLabel("Use offline maps")
        # read in the offline_director from settings if there is one
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('LocationWidget')
        odd = settings.value('offline_maps_dir', '')
        odd_isChecked_str = settings.value('use_offline_cb', 'False')
        if type(odd_isChecked_str) is str:
            odd_isChecked = odd_isChecked_str.lower() == 'true'
        else:
            odd_isChecked = odd_isChecked_str
        settings.endGroup()

        self.offline_directory_label.setText(odd)
        # for now, if there is a directory in the offline_directory_label, assume they want to use that, and activate checkbox
        self.offline_checkbox.setChecked(odd_isChecked)
        self.offline_directory_label.setEnabled(odd_isChecked)

        self.offline_directory_select_button = QPushButton("Select Folder...")
        self.offline_directory_select_button.setEnabled(odd_isChecked)

        self.offline_file_dialog = QFileDialog()
        self.offline_file_dialog.setFileMode(QFileDialog.Directory)

        offdir_layout = QHBoxLayout()
        offdir_layout.addWidget(self.offline_checkbox)
        offdir_layout.addWidget(self.offline_directory_select_button)
        offline_layout = QVBoxLayout()
        offline_layout.addLayout(offdir_layout)
        offline_layout.addWidget(self.offline_directory_label)


        options_gb = QGroupBox("Options")
        options_layout = QVBoxLayout()
        options_layout.addLayout(resolution_layout)
        options_layout.addWidget(self.backgroud_image_checkbox)
        options_layout.addLayout(offline_layout)
        options_layout.addLayout(colors_layout)
        options_layout.addStretch()
        options_gb.setLayout(options_layout)

        ### extent ###
        self.extent_settings = IPExtentSettingsWidget(title="Extent")

        ### bisl ###
        self.bisl_settings = IPBISLSettingsWidget(parent=self, title="BISL")

        ###   layouts   ###
        boxes_layout = QHBoxLayout()
        boxes_layout.addWidget(features_gb)
        boxes_layout.addWidget(options_gb)
        boxes_layout.addWidget(self.extent_settings)
        boxes_layout.addWidget(self.bisl_settings)

        main_layout = QHBoxLayout()
        main_layout.addLayout(boxes_layout)
        main_layout.addStretch()
        main_layout.setContentsMargins(0,0,0,0)
        self.setLayout(main_layout)

        self.connect_signals_and_slots()

    def to_dict(self):
        s_dict = {}
        s_dict['gui_borders'] = self.borders_checkbox.isChecked()
        s_dict['gui_states'] = self.states_checkbox.isChecked()
        s_dict['gui_lakes'] = self.lakes_checkbox.isChecked()
        s_dict['gui_rivers'] = self.rivers_checkbox.isChecked()
        s_dict['gui_coastline'] = self.coast_checkbox.isChecked()
        s_dict['gui_show_grid'] = self.show_grid_checkbox.isChecked()

        s_dict['gui_ocean_color'] = self.ocean_color_button.color_str()
        s_dict['gui_land_color'] = self.land_color_button.color_str()
        s_dict['gui_resolution'] = self.resolution_cb.currentText()
        s_dict['gui_use_background'] = self.backgroud_image_checkbox.isChecked()
        s_dict['gui_offline'] = self.offline_checkbox.isChecked()
        s_dict['gui_offline_dir'] = self.offline_directory_label.text()



        s_dict['extent_dict'] = self.extent_settings.to_dict()

        s_dict['bisl_dict'] = self.bisl_settings.to_dict()

        return s_dict

    def from_dict(self, s_dict):
        bool_map = {"True": True, "False": False}
        sd = s_dict['location_widget']
        self.bisl_settings.from_dict(sd)

        self.borders_checkbox.setChecked(bool_map[sd['gui_borders']])
        self.states_checkbox.setChecked(bool_map[sd['gui_states']])
        self.lakes_checkbox.setChecked(bool_map[sd['gui_lakes']])
        self.rivers_checkbox.setChecked(bool_map[sd['gui_rivers']])
        self.coast_checkbox.setChecked(bool_map[sd['gui_coastline']])

        self.ocean_color_button.set_color_from_str(sd['gui_ocean_color'])
        self.land_color_button.set_color_from_str(sd['gui_land_color'])
        self.resolution_cb.setCurrentText(sd['gui_resolution'])
        self.backgroud_image_checkbox.setChecked(bool_map[sd['gui_use_background']])
        self.offline_checkbox.setChecked(bool_map[sd['gui_offline']])
        self.offline_directory_label.setText(sd['gui_offline_dir'])

        self.extent_settings.from_dict(sd)

        self.signal_colors_changed.emit()

    def connect_signals_and_slots(self):
        self.offline_checkbox.clicked.connect(self.offline_directory_select_button.setEnabled)
        self.offline_checkbox.clicked.connect(self.offline_directory_label.setEnabled)
        self.offline_checkbox.clicked.connect(self.update_settings)

        self.ocean_color_button.clicked.connect(self.update_ocean_color)
        self.land_color_button.clicked.connect(self.update_land_color)

        self.backgroud_image_checkbox.clicked.connect(self.toggle_background_image)
        self.show_grid_checkbox.clicked.connect(self.update_grid_lines)
        self.offline_directory_select_button.clicked.connect(self.select_offline_maps_directory)

    def toggle_background_image(self):
        self.land_color_button.setDisabled(self.backgroud_image_checkbox.isChecked())
        self.ocean_color_button.setDisabled(self.backgroud_image_checkbox.isChecked())
        self.signal_background_changed.emit()

    def update_grid_lines(self):
        self.signal_map_settings_changed.emit()

    def update_ocean_color(self):
        new_color = QColorDialog.getColor(self.ocean_color_button.color())
        if new_color.isValid():
            self.ocean_color_button.set_color(new_color)
            self.signal_colors_changed.emit()

    def update_land_color(self):
        new_color = QColorDialog.getColor(self.land_color_button.color())
        if new_color.isValid(): 
            self.land_color_button.set_color(new_color)
            self.signal_colors_changed.emit()

    def select_offline_maps_directory(self):
        curr_dir = self.offline_directory_label.text()
        
        new_dir = QFileDialog.getExistingDirectory()
        
        self.offline_directory_label.setText(new_dir) 
        self.signal_offline_directory_changed.emit()
        
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('LocationWidget')
        settings.setValue('offline_maps_dir', new_dir)
        settings.endGroup()

    def update_settings(self):
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('LocationWidget')
        settings.setValue('use_offline_cb', self.offline_checkbox.isChecked())
        settings.endGroup()


class IPBISLSettingsWidget(IPBaseWidgets.IPSettingsGroupBox):

    def __init__(self, title="", parent=None):
        super().__init__(title=title, parent=parent)

        self.earth_radius = 6378.1   # km

        self.parent = parent
        self.setTitle(self.tr(title))
        self.buildUI()

    def buildUI(self):

        # self.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)

        self.bm_width_edit = QDoubleSpinBox()
        self.bm_width_edit.setMinimum(2.5)
        self.bm_width_edit.setMaximum(45.0)
        self.bm_width_edit.setValue(10)
        self.bm_width_edit.setSuffix(' deg')
        self.bm_width_edit.valueChanged.connect(self.enable_update_dm_button)

        self.rng_max_edit = QSpinBox()
        self.rng_max_edit.setMinimum(100)
        self.rng_max_edit.setSingleStep(100)
        self.rng_max_edit.setMaximum(int(np.pi * self.earth_radius))
        self.rng_max_edit.setValue(3000)
        self.rng_max_edit.setSuffix(' km')
        self.rng_max_edit.valueChanged.connect(self.enable_update_dm_button)

        self.resolution_edit = QDoubleSpinBox()
        self.resolution_edit.setMinimum(.01)
        self.resolution_edit.setMaximum(10)
        self.resolution_edit.setValue(.05)
        self.resolution_edit.valueChanged.connect(self.enable_update_dm_button)

        self.tm_resolution_edit = QSpinBox()
        self.tm_resolution_edit.setMinimum(1)
        self.tm_resolution_edit.setMaximum(600)
        self.tm_resolution_edit.setValue(60)
        self.tm_resolution_edit.valueChanged.connect(self.enable_update_dm_button)

        self.confidence_edit = QSpinBox()
        self.confidence_edit.setMinimum(1)
        self.confidence_edit.setMaximum(99)
        self.confidence_edit.setValue(95)
        self.confidence_edit.setSuffix(' %')

        self.update_dm_button = QPushButton('Update Dist. Matrix')
        self.update_dm_button.setEnabled(False)

        layout1 = QFormLayout()
        layout1.addRow(self.tr('Beam Width: '), self.bm_width_edit)
        layout1.addRow(self.tr('Range Max.: '), self.rng_max_edit)
        layout1.addRow(self.tr('Lat/Lon Resolution'), self.resolution_edit)
        layout2 = QFormLayout()

        layout2.addRow(self.tr('Time Resolution'), self.tm_resolution_edit)
        layout2.addRow(self.tr('Confidence'), self.confidence_edit)
        layout2.addRow("", self.update_dm_button)

        layout = QHBoxLayout()
        layout.addLayout(layout1)
        layout.addLayout(layout2)
        self.setLayout(layout)

    def to_dict(self):
        s_dict = {}
        s_dict['bisl_bm_width'] = self.bm_width_edit.value()
        s_dict['bisl_rng_max'] = self.rng_max_edit.value()
        s_dict['bisl_resolution'] = self.resolution_edit.value()
        s_dict['bisl_t_resolution'] = self.tm_resolution_edit.value()
        s_dict['bisl_confidence'] = self.confidence_edit.value()
        return s_dict

    def from_dict(self, s_dict):
        sd = s_dict['bisl_dict']
        self.bm_width_edit.setValue(float(sd['bisl_bm_width']))
        self.rng_max_edit.setValue(int(float(sd['bisl_rng_max'])))
    
    @pyqtSlot(float)
    @pyqtSlot(int)
    def enable_update_dm_button(self, _):
        self.parent.bisl_settings.update_dm_button.setEnabled(True)


class IPExtentSettingsWidget(IPBaseWidgets.IPSettingsGroupBox):

    sig_extent_changed = pyqtSignal(list)
    sig_set_to_global = pyqtSignal()
    sig_autoscale = pyqtSignal()

    def __init__(self, title="", parent=None):
        super().__init__(title=title, parent=parent)

        self.setTitle(self.tr(title))
        self.buildUI()

    def buildUI(self):

        ll_label = QLabel(" Lower left:")
        self.ll_lat_spin = QDoubleSpinBox()
        self.ll_lat_spin.setMaximumWidth(100)
        self.ll_lat_spin.setRange(-90.0, 90.0)
        self.ll_lat_spin.setValue(-90.0)
        self.ll_lat_spin.setPrefix("Lat: ")
        self.ll_lat_spin.valueChanged.connect(self.activate_update_button)

        self.ll_lon_spin = QDoubleSpinBox()
        self.ll_lon_spin.setMaximumWidth(100)
        self.ll_lon_spin.setRange(-179.99, 180.0)
        self.ll_lon_spin.setValue(-179.99)
        self.ll_lon_spin.setPrefix("Lon: ")
        self.ll_lon_spin.valueChanged.connect(self.activate_update_button)

        ll_layout = QHBoxLayout()
        ll_layout.addWidget(ll_label)
        ll_layout.addWidget(self.ll_lon_spin)
        ll_layout.addWidget(self.ll_lat_spin)


        ur_label = QLabel(" Upper right:")
        self.ur_lat_spin = QDoubleSpinBox()
        self.ur_lat_spin.setMaximumWidth(100)
        self.ur_lat_spin.setRange(-90., 90.0)
        self.ur_lat_spin.setValue(90.0)
        self.ur_lat_spin.setPrefix("Lat: ")
        self.ur_lat_spin.valueChanged.connect(self.activate_update_button)

        self.ur_lon_spin = QDoubleSpinBox()
        self.ur_lon_spin.setMaximumWidth(100)
        self.ur_lon_spin.setRange(-179.99, 180.0)
        self.ur_lon_spin.setValue(180.0)
        self.ur_lon_spin.setPrefix("Lon: ")
        self.ur_lon_spin.valueChanged.connect(self.activate_update_button)

        ur_layout = QHBoxLayout()
        ur_layout.addWidget(ur_label)
        ur_layout.addWidget(self.ur_lon_spin)
        ur_layout.addWidget(self.ur_lat_spin)

        self.update_plot_button = QPushButton("Update")
        self.update_plot_button.setMaximumWidth(100)
        self.update_plot_button.setEnabled(False)
        self.update_plot_button.clicked.connect(self.deactivate_update_button)
        self.update_plot_button.clicked.connect(self.update_map_extent)

        self.set_to_global_button = QPushButton("Global")
        self.set_to_global_button.setMaximumWidth(100)
        self.set_to_global_button.clicked.connect(self.set_to_global)

        self.autoscale_button = QPushButton("Autoscale")
        self.autoscale_button.setMaximumWidth(100)
        self.autoscale_button.clicked.connect(self.autoscale_map)

        self.hide_button = QPushButton("Hide")
        self.hide_button.setMaximumWidth(60)

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.set_to_global_button)
        button_layout.addWidget(self.autoscale_button)
        button_layout.addWidget(self.update_plot_button)

        coord_layout = QVBoxLayout()
        coord_layout.addLayout(ll_layout)
        coord_layout.addLayout(ur_layout)
        coord_layout.addLayout(button_layout)
        coord_layout.addStretch()
        self.setLayout(coord_layout)

    def to_dict(self):
        s_dict = {}
        s_dict['ll_lat'] = self.ll_lat_spin.value()
        s_dict['ll_lon'] = self.ll_lon_spin.value()
        s_dict['ur_lat'] = self.ur_lat_spin.value()
        s_dict['ur_lon'] = self.ur_lon_spin.value()
        return s_dict
    
    def from_dict(self, s_dict):
        sd = s_dict['extent_dict']
        self.ll_lat_spin.setValue(float(sd['ll_lat']))
        self.ll_lon_spin.setValue(float(sd['ll_lon']))
        self.ur_lat_spin.setValue(float(sd['ur_lat']))
        self.ur_lon_spin.setValue(float(sd['ur_lon']))

        self.update_map_extent()
 
    def set_extent_spin_values(self, extent):
        # ll_lon: lower left longitude
        # ur_lat: upper right latitude
        # etc

        self.ll_lon_spin.setValue(extent[0])
        self.ll_lat_spin.setValue(extent[2])
        self.ur_lon_spin.setValue(extent[1])
        self.ur_lat_spin.setValue(extent[3])

    def set_to_global(self):
        self.sig_set_to_global.emit()

    def autoscale_map(self):
        self.sig_autoscale.emit()

    def activate_update_button(self):
        self.update_plot_button.setEnabled(True)

    def deactivate_update_button(self):
        self.update_plot_button.setEnabled(False)

    def update_map_extent(self):
        extent = [self.ll_lon_spin.value(), self.ur_lon_spin.value(), self.ll_lat_spin.value(), self.ur_lat_spin.value()]
        self.sig_extent_changed.emit(extent)