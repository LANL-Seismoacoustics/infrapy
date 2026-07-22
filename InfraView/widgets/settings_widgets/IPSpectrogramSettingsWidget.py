from PyQt5.QtWidgets import QPushButton, QComboBox, QLabel, QRadioButton, QFormLayout, QGroupBox, QHBoxLayout, QSpinBox, QDoubleSpinBox
from PyQt5.QtCore import pyqtSlot, pyqtSignal

import pyqtgraph as pg

from InfraView.widgets import IPBaseWidgets

class IPSpectrogramSettingsWidget(IPBaseWidgets.IPSettingsWidget):

    sig_spectrogram_changed = pyqtSignal()
  
    def __init__(self, parent):
        super().__init__(parent)
        self.last_spec_type = "Spectrogram"
        self.buildUI()

    def buildUI(self):
        ############ Spectrogram settings ##############
        spec_gb = QGroupBox("Spectrogram ")
        spec_layout = QHBoxLayout()

        self.update_button = QPushButton("Update")
        self.update_button.setMaximumWidth(100)
        self.update_button.setEnabled(False)

        self.spec_type_cb = QComboBox(self)
        spec_type_label = QLabel("Spectrogram Type:")
        self.spec_type_cb.addItem("SPECTROGRAM")
        self.spec_type_cb.addItem("STFT")
        self.spec_type_cb.addItem("CWT")

        #TODO Are we going to hardcode this value or what?
        # Currently it's hidden
        self.omega0_label = QLabel("    omega0")
        self.omega0_spin = QSpinBox()
        self.omega0_spin.setMaximumWidth(100)
        self.omega0_spin.setMinimum(1)
        self.omega0_spin.setMaximum(100)
        self.omega0_spin.setValue(7)
        self.omega0_spin.setToolTip("Used by the CWT.  Number of periods in the wavelet.")
        self.omega0_label.setVisible(False)
        self.omega0_spin.setVisible(False)

        colormap_label = QLabel("Color Map: ")
        self.colormap_cb = QComboBox()
        available_maps = pg.colormap.listMaps(source='matplotlib')
        self.colormap_cb.addItems(available_maps)
        self.colormap_cb.setCurrentText('jet')

        self.colorbar_rb = QRadioButton('')
        

        form1_layout = QFormLayout()
        form1_layout.addRow(self.spec_type_cb, self.update_button)
        form1_layout.addRow(self.omega0_label, self.omega0_spin)
        form1_layout.addRow(colormap_label, self.colormap_cb)
        form1_layout.addRow("Color Bar: ", self.colorbar_rb)
        
        spec_layout.addLayout(form1_layout)
        spec_gb.setLayout(spec_layout)

        ############ Detector settings ##############
        self.detector_gb = QGroupBox("Detector")

        pval_label = QLabel('pval: ')
        self.pval_spin = QDoubleSpinBox()
        self.pval_spin.setDecimals(4)
        self.pval_spin.setMinimum(0.001)
        self.pval_spin.setMaximum(1.0)
        self.pval_spin.setValue(0.01)
        self.pval_spin.setSingleStep(0.001)
        self.pval_spin.setMinimumWidth(70)
        self.pval_spin.setMaximumWidth(70)

        fmin_label = QLabel("Freq Min: ")
        self.fmin_spin = QDoubleSpinBox()
        self.fmin_spin.setMinimum(0.001)
        self.fmin_spin.setMaximum(10000.0)
        self.fmin_spin.setValue(10000)  # Set to strange value so will be guaranteed to trigger a valueChange when loaded.
        self.fmin_spin.setMinimumWidth(70)
        self.fmin_spin.setMaximumWidth(70)

        fmax_label = QLabel("Freq Max: ")
        self.fmax_spin = QDoubleSpinBox()
        self.fmax_spin.setMinimum(0.002)
        self.fmax_spin.setMaximum(10000.0)
        self.fmax_spin.setValue(0.002)    # this needs to be set when a spectrogram is created
        self.fmax_spin.setMinimumWidth(70)
        self.fmax_spin.setMaximumWidth(70)

        clust_freq_scale_label = QLabel("Cluster Freq. Scaling: ")
        self.clust_freq_scale_spin = QDoubleSpinBox()
        self.clust_freq_scale_spin.setMinimum(1.0)
        self.clust_freq_scale_spin.setMaximum(10000.0)
        self.clust_freq_scale_spin.setValue(35.0)
        self.clust_freq_scale_spin.setMinimumWidth(70)
        self.clust_freq_scale_spin.setMaximumWidth(70)

        clust_eps_label = QLabel("Clustering EPS: ")
        self.clust_eps_spin = QSpinBox()
        self.clust_eps_spin.setMinimum(1)
        self.clust_eps_spin.setMaximum(10000)
        self.clust_eps_spin.setValue(10)
        self.clust_eps_spin.setMinimumWidth(70)
        self.clust_eps_spin.setMaximumWidth(70)

        clust_min_samples_label = QLabel("Clustering Min Samples: ")
        self.clust_min_samples_spin = QSpinBox()
        self.clust_min_samples_spin.setMinimum(1)
        self.clust_min_samples_spin.setMaximum(10000)
        self.clust_min_samples_spin.setValue(40)
        self.clust_min_samples_spin.setMinimumWidth(70)
        self.clust_min_samples_spin.setMaximumWidth(70)

        adaptive_win_len_label = QLabel("Adaptive Window Length (s)")
        self.adaptive_win_len_spin = QSpinBox()
        self.adaptive_win_len_spin.setMinimum(1)
        self.adaptive_win_len_spin.setMaximum(10000)
        self.adaptive_win_len_spin.setValue(600)
        self.adaptive_win_len_spin.setMinimumWidth(70)
        self.adaptive_win_len_spin.setMaximumWidth(70)

        form2_layout = QFormLayout()
        form2_layout.addRow(pval_label, self.pval_spin)
        form2_layout.addRow(fmin_label, self.fmin_spin)
        form2_layout.addRow(fmax_label, self.fmax_spin)

        form3_layout = QFormLayout()
        form3_layout.addRow(clust_freq_scale_label, self.clust_freq_scale_spin)
        form3_layout.addRow(clust_eps_label, self.clust_eps_spin)
        form3_layout.addRow(clust_min_samples_label, self.clust_min_samples_spin)
        form3_layout.addRow(adaptive_win_len_label, self.adaptive_win_len_spin)
        
        det_layout = QHBoxLayout()
        det_layout.addLayout(form2_layout)
        det_layout.addLayout(form3_layout)

        self.detector_gb.setLayout(det_layout)

        #######################
        self.cwt_detector_gb = QGroupBox("CWT Detector")
        self.cwt_detector_gb.setVisible(False)
        cwt_pval_label = QLabel('pval: ')
        self.cwt_pval_spin = QDoubleSpinBox()
        self.cwt_pval_spin.setDecimals(4)
        self.cwt_pval_spin.setMinimum(0.001)
        self.cwt_pval_spin.setMaximum(1.0)
        self.cwt_pval_spin.setValue(0.01)
        self.cwt_pval_spin.setSingleStep(0.001)
        self.cwt_pval_spin.setMinimumWidth(70)
        self.cwt_pval_spin.setMaximumWidth(70)

        cwt_fmin_label = QLabel("Freq min: ")
        self.cwt_fmin_spin = QDoubleSpinBox()
        self.cwt_fmin_spin.setMinimum(0.001)
        self.cwt_fmin_spin.setMaximum(10000.0)
        self.cwt_fmin_spin.setValue(10000.0)    # Set to strange value so will be guaranteed to trigger a valueChange when loaded.
        self.cwt_fmin_spin.setMinimumWidth(70)
        self.cwt_fmin_spin.setMaximumWidth(70)

        cwt_fmax_label = QLabel("Freq max: ")
        self.cwt_fmax_spin = QDoubleSpinBox()
        self.cwt_fmax_spin.setMinimum(1.0)
        self.cwt_fmax_spin.setMaximum(10000.0)
        self.cwt_fmax_spin.setValue(10.0)    # this needs to be set when a spectrogram is created
        self.cwt_fmax_spin.setMinimumWidth(70)
        self.cwt_fmax_spin.setMaximumWidth(70)

        cwt_clust_freq_scale_label = QLabel("Cluster Freq. Scaling: ")
        self.cwt_clust_freq_scale_spin = QDoubleSpinBox()
        self.cwt_clust_freq_scale_spin.setMinimum(1.0)
        self.cwt_clust_freq_scale_spin.setMaximum(10000.0)
        self.cwt_clust_freq_scale_spin.setValue(35.0)
        self.cwt_clust_freq_scale_spin.setMinimumWidth(70)
        self.cwt_clust_freq_scale_spin.setMaximumWidth(70)

        cwt_clust_eps_label = QLabel("Clustering EPS: ")
        self.cwt_clust_eps_spin = QSpinBox()
        self.cwt_clust_eps_spin.setMinimum(1)
        self.cwt_clust_eps_spin.setMaximum(10000)
        self.cwt_clust_eps_spin.setValue(5)
        self.cwt_clust_eps_spin.setMinimumWidth(70)
        self.cwt_clust_eps_spin.setMaximumWidth(70)

        cwt_clust_min_samples_label = QLabel("Clustering Min Samples: ")
        self.cwt_clust_min_samples_spin = QSpinBox()
        self.cwt_clust_min_samples_spin.setMinimum(1)
        self.cwt_clust_min_samples_spin.setMaximum(10000)
        self.cwt_clust_min_samples_spin.setValue(500)
        self.cwt_clust_min_samples_spin.setMinimumWidth(70)
        self.cwt_clust_min_samples_spin.setMaximumWidth(70)

        cwt_adaptive_win_len_label = QLabel("Adaptive window length (s)")
        self.cwt_adaptive_win_len_spin = QSpinBox()
        self.cwt_adaptive_win_len_spin.setMinimum(1)
        self.cwt_adaptive_win_len_spin.setMaximum(10000)
        self.cwt_adaptive_win_len_spin.setValue(600)
        self.cwt_adaptive_win_len_spin.setMinimumWidth(70)
        self.cwt_adaptive_win_len_spin.setMaximumWidth(70)

        form4_layout = QFormLayout()
        form4_layout.addRow(cwt_pval_label, self.cwt_pval_spin)
        form4_layout.addRow(cwt_fmin_label, self.cwt_fmin_spin)
        form4_layout.addRow(cwt_fmax_label, self.cwt_fmax_spin)

        form5_layout = QFormLayout()
        form5_layout.addRow(cwt_clust_freq_scale_label, self.cwt_clust_freq_scale_spin)
        form5_layout.addRow(cwt_clust_eps_label, self.cwt_clust_eps_spin)
        form5_layout.addRow(cwt_clust_min_samples_label, self.cwt_clust_min_samples_spin)
        form5_layout.addRow(cwt_adaptive_win_len_label, self.cwt_adaptive_win_len_spin)

        cwt_det_layout = QHBoxLayout()
        cwt_det_layout.addLayout(form4_layout)
        cwt_det_layout.addLayout(form5_layout)

        self.cwt_detector_gb.setLayout(cwt_det_layout)

        self.hide_button = QPushButton("Hide")
        self.hide_button.setMaximumWidth(60)
        self.hide_button.clicked.connect(self.hide)

        h_layout = QHBoxLayout()
        h_layout.addWidget(spec_gb)
        h_layout.addWidget(self.detector_gb)       
        h_layout.addWidget(self.cwt_detector_gb) 
        h_layout.addStretch()
        h_layout.setContentsMargins(0,0,0,0)
        self.setLayout(h_layout) 

    def to_dict(self):
        s_dict = {}
        s_dict['spec_type'] = self.spec_type_cb.currentText().lower()
        s_dict['gui_colormap'] = self.colormap_cb.currentText()
        s_dict['gui_colorbar'] = self.colorbar_rb.isChecked()

        s_dict['omega0'] = self.omega0_spin.value()
        s_dict['pval'] = self.pval_spin.value()
        s_dict['fmin'] = self.fmin_spin.value()
        s_dict['fmax'] = self.fmax_spin.value()
        s_dict['cluster_freq_scale'] = self.clust_freq_scale_spin.value()
        s_dict['cluster_eps'] = self.clust_eps_spin.value()
        s_dict['cluster_min_samples'] = self.clust_min_samples_spin.value()
        s_dict['adapt_win_len'] = self.adaptive_win_len_spin.value()

        s_dict['cwt_pval'] = self.cwt_pval_spin.value()
        s_dict['cwt_fmin'] = self.cwt_fmin_spin.value()
        s_dict['cwt_fmax'] = self.cwt_fmax_spin.value()
        s_dict['cwt_cluster_freq_scale'] = self.cwt_clust_freq_scale_spin.value()
        s_dict['cwt_cluster_eps'] = self.cwt_clust_eps_spin.value()
        s_dict['cwt_cluster_min_samples'] = self.cwt_clust_min_samples_spin.value()
        s_dict['cwt_adapt_win_len'] = self.cwt_adaptive_win_len_spin.value()

        return s_dict

    def from_dict(self, s_dict):
        bool_map = {"True": True, "False": False}

        sd = s_dict['spectral_widget']

        self.spec_type_cb.setCurrentText(sd['spec_type'].upper())
        self.colormap_cb.setCurrentText(sd['gui_colormap'])
        self.colorbar_rb.setChecked(bool_map[sd['gui_colorbar']])
        
        self.omega0_spin.setValue(int(float(sd['omega0'])))

        self.pval_spin.setValue(float(sd['pval']))
        self.fmin_spin.setValue(float(sd['fmin']))
        self.fmax_spin.setValue(float(sd['fmax']))
        self.clust_freq_scale_spin.setValue(float(sd['cluster_freq_scale']))
        self.clust_eps_spin.setValue(int(float(sd['cluster_eps'])))
        self.clust_min_samples_spin.setValue(int(sd['cluster_min_samples']))
        self.adaptive_win_len_spin.setValue(int(float(sd['adapt_win_len'])))

        self.cwt_pval_spin.setValue(float(sd['cwt_pval']))
        self.cwt_fmin_spin.setValue(float(sd['cwt_fmin']))
        self.cwt_fmax_spin.setValue(float(sd['cwt_fmax']))
        self.cwt_clust_freq_scale_spin.setValue(float(sd['cwt_cluster_freq_scale']))
        self.cwt_clust_eps_spin.setValue(int(float(sd['cwt_cluster_eps'])))
        self.cwt_clust_min_samples_spin.setValue(int(sd['cwt_cluster_min_samples']))
        self.cwt_adaptive_win_len_spin.setValue(int(float(sd['cwt_adapt_win_len'])))

        self.update_button.clicked.emit()

    def set_controlled_widget(self, widget):
        super().set_controlled_widget(widget)
        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        # This method is typically called exactly once.  set_controlled_widget MUST be called (successfully) before you call this.
        
        self.update_button.clicked.connect(self.deactivate_update_button)

        self.spec_type_cb.currentTextChanged.connect(self.activate_update_button)
        self.spec_type_cb.currentTextChanged.connect(self.detector_selector)

        #TODO these should update the Run Detector button, not the update button --- WAIT, Actually this might not be correct...
        self.pval_spin.valueChanged.connect(self.activate_update_button)

        self.fmin_spin.valueChanged.connect(self.activate_update_button)
        self.fmin_spin.valueChanged.connect(self.controlled_widget.signalSpecWidget.fmin_changed)

        self.fmax_spin.valueChanged.connect(self.activate_update_button)
        self.fmax_spin.valueChanged.connect(self.controlled_widget.signalSpecWidget.fmax_changed)

        self.clust_freq_scale_spin.valueChanged.connect(self.activate_update_button)
        self.clust_eps_spin.valueChanged.connect(self.activate_update_button)
        self.clust_min_samples_spin.valueChanged.connect(self.activate_update_button)

        #TODO CWT should also update the Run Detector button --- WAIT, Actually this might not be correct...
        self.cwt_pval_spin.valueChanged.connect(self.activate_update_button)

        self.cwt_fmin_spin.valueChanged.connect(self.activate_update_button)
        self.cwt_fmin_spin.valueChanged.connect(self.controlled_widget.signalSpecWidget.fmin_changed)

        self.cwt_fmax_spin.valueChanged.connect(self.activate_update_button)
        self.cwt_fmax_spin.valueChanged.connect(self.controlled_widget.signalSpecWidget.fmax_changed)
        
        self.cwt_clust_freq_scale_spin.valueChanged.connect(self.activate_update_button)
        self.cwt_clust_eps_spin.valueChanged.connect(self.activate_update_button)
        self.cwt_clust_min_samples_spin.valueChanged.connect(self.activate_update_button)

    @pyqtSlot(str)
    def detector_selector(self, det_type):
        if det_type == 'CWT':
            self.cwt_detector_gb.setVisible(True)
            self.detector_gb.setVisible(False)
        else:
            self.detector_gb.setVisible(True)
            self.cwt_detector_gb.setVisible(False)

    def get_scale_setting(self):

        if self.colorbar_rb.isChecked():
            return 'cbar'
        else:
            return 'none'
        

    def activate_update_button(self):
        # self.omega0_label.setVisible(self.spec_type_cb.currentText() == 'CWT')
        # self.omega0_spin.setVisible(self.spec_type_cb.currentText() == 'CWT')
            
        self.update_button.setEnabled(True)

    def deactivate_update_button(self):
        self.update_button.setEnabled(False)

    @pyqtSlot()
    def hide(self):
        self.setVisible(False)