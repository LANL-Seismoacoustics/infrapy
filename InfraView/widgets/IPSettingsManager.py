from PyQt5.QtWidgets import QStackedWidget, QWidget, QPushButton, QHBoxLayout, QVBoxLayout, QFrame, QFileDialog
from PyQt5.QtCore import Qt, pyqtSlot

import traceback, yaml, errno, os
import configparser

from pathlib import Path

from InfraView.widgets.settings_widgets import IPBeamformingSettingsWidget
from InfraView.widgets.settings_widgets import IPSpectrogramSettingsWidget
from InfraView.widgets.settings_widgets import IPLocationSettingsWidget
from InfraView.widgets.settings_widgets import IPWaveformSettingsWidget
from InfraView.widgets.settings_widgets import IPDatabaseSettingsWidget


class IPSettingsManager(QFrame):
    def __init__(self, parent):
        super().__init__(parent)

        self.setObjectName("settingsManager")
        self.setWindowFlags(Qt.FramelessWindowHint)
        self.setStyleSheet("#settingsManager {border: 1px solid #777;} ")

        self.default_config_path = Path(__file__).parent.parent.parent /"infrapy" / "resources" / "default.config"
        if self.default_config_path.exists():
            self.current_config_path = self.default_config_path
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), str(self.default_config_path))

        self.settings_widget_dict = {}

        self.initialize_settings_widgets()

        self.settings_stack = QStackedWidget(self)
        self.insert_settings_widgets()

        hide_layout = QVBoxLayout()
        # self.hide_button = QPushButton('Hide')
        # self.hide_button.clicked.connect(self.hide)
        self.load_button = QPushButton('Load')
        self.load_button.clicked.connect(self.load_settings)
        self.save_button = QPushButton('Save')
        self.save_button.clicked.connect(self.save_settings)
        self.save_as_button = QPushButton('Save As...')
        self.save_as_button.clicked.connect(self.save_settings_as)
        self.load_defaults_button = QPushButton('Defaults')
        self.load_defaults_button.clicked.connect(self.load_default_settings)
        self.load_defaults_button.setToolTip("Load default settings")
        # hide_layout.addWidget(self.hide_button)
        hide_layout.addWidget(self.load_button)
        hide_layout.addWidget(self.load_defaults_button)
        hide_layout.addWidget(self.save_button)
        hide_layout.addWidget(self.save_as_button)
        hide_layout.addStretch()

        layout = QHBoxLayout()
        layout.addWidget(self.settings_stack)
        layout.addStretch()
        layout.addLayout(hide_layout)

        self.setLayout(layout)

        self.setVisible(False)

    def hide(self):
        self.setVisible(False)

    def initialize_settings_widgets(self):
        
        self.settings_widget_dict['database'] = IPDatabaseSettingsWidget.IPDatabaseSettingsWidget(self)
        self.settings_widget_dict['waveforms'] = IPWaveformSettingsWidget.IPWaveformSettingsWidget(self)
        self.settings_widget_dict['spectral'] = IPSpectrogramSettingsWidget.IPSpectrogramSettingsWidget(self)
        self.settings_widget_dict['location'] = IPLocationSettingsWidget.IPLocationSettingsWidget(parent=self)
        self.settings_widget_dict['beamforming'] = IPBeamformingSettingsWidget.IPBeamformingSettingsWidget(self)


    def connect_widgets_and_settings(self, widget_dict):
        # settings widgets need to have a reference to the widgets they control, and vice-versa
        for key, value in widget_dict.items():
            try:
                self.settings_widget_dict[key].set_controlled_widget(value)
            except  KeyError:
                pass    # probably the App window, which isn't a settings widget

        for key, value in self.settings_widget_dict.items():
            try: 
                widget_dict[key].set_controlling_widget(value)
            except Exception as e:
                import traceback
                traceback.print_exc()

                print("{} doesn't have set_controlling_widget method yet".format(type(widget_dict[key])))

    def insert_settings_widgets(self):
        for _, value in self.settings_widget_dict.items():
            self.settings_stack.addWidget(value)

        self.settings_stack.setCurrentWidget(self.settings_widget_dict['waveforms'])

    @pyqtSlot(str)
    def widget_changed(self, widget_name):
        # someone clicked a action to change the active, so we need to change the settings widget to match
        try:
            self.settings_stack.setCurrentWidget(self.settings_widget_dict[widget_name.lower()])
        except KeyError:
            print("{} settings not found".format(widget_name))

    def toggle_visibility(self):
        self.setVisible(self.isHidden())


    def set_current_config(self, new_path):
        self.current_config_path = new_path
        self.save_button.setToolTip(str(new_path))

    @pyqtSlot()
    def save_settings(self):
        # filepath is a pathlib.Path representing the absolute path + filname to the file

        # collect settings dictionaries from all of the widgets, and save them to an infrapy file
        settings_dict = {}
        for key, value in self.settings_widget_dict.items():
            try:
                settings_dict[key+'_widget'] = value.to_dict()
            except AttributeError as e:
                #print(traceback.format_exc())
                print("{} doesn't have a to_dict method yet".format(value))
        cli_dict = self.map_infraview_settings_to_cli(settings_dict)

        # if self.current_cnfg is the self.default_cnfg, prompt for a new filename. (ie don't allow anyone to write over default config)
        # otherwise, save to the current config path
        if self.current_config_path == self.default_config_path:
            filename = QFileDialog.getSaveFileName(self, "Config File", "", "Config Files (*.ini *.config)" )
            self.set_current_config(Path(filename[0]))

        self.dict_to_ini(cli_dict, self.current_config_path)

    @pyqtSlot()
    def save_settings_as(self):

        settings_dict = {}
        for key, value in self.settings_widget_dict.items():
            try:
                settings_dict[key+'_widget'] = value.to_dict()
            except AttributeError as e:
                #print(traceback.format_exc())
                print("{} doesn't have a to_dict method yet".format(value))
        cli_dict = self.map_infraview_settings_to_cli(settings_dict)

        # force the filedialog
        filename = QFileDialog.getSaveFileName(self, "Config File", "", "Config Files (*.ini *.config)")

        self.set_current_config(Path(filename[0]))

        self.dict_to_ini(cli_dict, self.current_config_path)

    def load_default_settings(self):
        filename = str(self.default_config_path)
        settings_dict = self.ini_to_dict(filename)
        s_dict = self.map_cli_to_infraview_settings(settings_dict)
        self.set_current_config(self.default_config_path)

        for key, value in self.settings_widget_dict.items():
            try:
                settings_dict[key+'_widget'] = value.from_dict(s_dict)
            except AttributeError as e:
                print(traceback.format_exc())
                print("{} doesn't have a from_dict method yet".format(value))

        

    def load_settings(self):
        filename = QFileDialog.getOpenFileName(self, "Config File", "", "Config Files (*.ini *.config)" )[0]
        if not filename:
            return

        settings_dict = self.ini_to_dict(filename)
        s_dict = self.map_cli_to_infraview_settings(settings_dict)
        self.set_current_config(Path(filename))
                                
        for key, value in self.settings_widget_dict.items():
            try:
                settings_dict[key+'_widget'] = value.from_dict(s_dict)
            except AttributeError as e:
                print(traceback.format_exc())
                print("{} doesn't have a from_dict method yet".format(value))


    def ini_to_dict(self, ini_filename):
        # ini_filename is the absolute path to the ini file to read
        config = configparser.ConfigParser()
        config.read(ini_filename)
        
        ini_dict = {}
        for section in config.sections():
            ini_dict[section] = {}
            for option in config.options(section):
                ini_dict[section][option] = config.get(section, option)
        return ini_dict

    def dict_to_ini(self, dict, ini_filepath):
        config = configparser.ConfigParser()
        sections = dict.keys()

        for section in sections:
            config.add_section(section)

        for section in sections:
            sub_dict = dict[section]
            fields = sub_dict.keys()
            for field in fields:
                value = sub_dict[field]
                config.set(section, field, str(value))

        with open(str(ini_filepath), 'w') as f:
            config.write(f)

        self.current_config_path = ini_filepath

    def map_cli_to_infraview_settings(self, cli_dict):
        # Used whem reading in a configuration file.
        # The command line config files are read into a dictionary, this function
        # converts that dictionary to the one that the settings widget uses.
        settings_dict = {'waveforms_widget': {'filter_dict': {}, 'psd_dict': {}}, 
                         'beamforming_widget': {'detector_settings':{}},
                         'spectral_widget': {},
                         'location_widget': {'bisl_dict': {}, 'extent_dict': {}}}

        #FK
        settings_dict["waveforms_widget"]["filter_dict"]["highpass"] = float(cli_dict["FK"]["freq_min"])
        settings_dict["waveforms_widget"]["filter_dict"]["lowpass"] = float(cli_dict["FK"]["freq_max"])

        # I still want to write the default values for signal/noise windows into the config file for later reading
        settings_dict["beamforming_widget"]["signal_start"] = None
        settings_dict["beamforming_widget"]["signal_end"] = None
        settings_dict["beamforming_widget"]["noise_start"] = None
        settings_dict["beamforming_widget"]["noise_end"] = None

        settings_dict["beamforming_widget"]["backAz_start"] = cli_dict["FK"]["back_az_min"]
        settings_dict["beamforming_widget"]["backAz_end"] = cli_dict["FK"]["back_az_max"]
        settings_dict["beamforming_widget"]["backAz_resolution"] = cli_dict["FK"]["back_az_step"] 
        settings_dict["beamforming_widget"]["traceV_min"] = cli_dict["FK"]["trace_vel_min"]
        settings_dict["beamforming_widget"]["traceV_max"] = cli_dict["FK"]["trace_vel_max"]
        settings_dict["beamforming_widget"]["traceV_resolution"] = cli_dict["FK"]["trace_vel_step"]
        settings_dict["beamforming_widget"]["method"] = cli_dict["FK"]["method"]
        settings_dict["beamforming_widget"]["win_length"] = cli_dict["FK"]["window_len"]
        settings_dict["beamforming_widget"]["sub_win_length"] = cli_dict["FK"]["sub_window_len"]
        settings_dict["beamforming_widget"]["win_step"] = cli_dict["FK"]["window_step"]
        #cli_template_dict["FK"]["cpu_cnt"] # GUI doesn't save cpu cnt

        # FD node
        # settings_dict["beamforming_widget"]["detector_settings"]["window_len"] = cli_dict["FD"]["window_len"]
        settings_dict["beamforming_widget"]["detector_settings"]["thresh_ceil"] = cli_dict["FD"]["thresh_ceil"]
        settings_dict["beamforming_widget"]["detector_settings"]["pval"] = cli_dict["FD"]["p_value"]
        settings_dict["beamforming_widget"]["detector_settings"]["min_peak_width"] = cli_dict["FD"]["min_duration"]
        settings_dict["beamforming_widget"]["detector_settings"]["back_az_limit"] = cli_dict["FD"]["back_az_width"]
        settings_dict["beamforming_widget"]["detector_settings"]["manual_level"] = cli_dict["FD"]["fixed_thresh"]
        settings_dict["beamforming_widget"]["detector_settings"]["merge"] = cli_dict["FD"]["merge_dets"]

        # SD node
        # Spectrogram v. spectrogram - check string sanitization
        settings_dict["spectral_widget"]["spec_type"] = cli_dict["SD"]["spectral_option"]
        settings_dict["spectral_widget"]["omega0"] = cli_dict["SD"]["morlet_omega0"]
        settings_dict["spectral_widget"]["fmin"] = cli_dict["SD"]["freq_min"]
        settings_dict["spectral_widget"]["fmax"] = cli_dict["SD"]["freq_max"]
        settings_dict["spectral_widget"]["adapt_win_len"] = cli_dict["SD"]["window_len"]  
        # cli_dict["SD"]["window_step"].  GUI autmatically calulates window_step from window length
        settings_dict["spectral_widget"]["pval"] = cli_dict["SD"]["p_value"] 
        settings_dict["spectral_widget"]["cluster_freq_scale"] = cli_dict["SD"]["freq_tm_factor"] 
        settings_dict["spectral_widget"]["cluster_eps"] = cli_dict["SD"]["cluster_eps"]
        settings_dict["spectral_widget"]["cluster_min_samples"] = cli_dict["SD"]["cluster_min_samples"]

        settings_dict["spectral_widget"]["cwt_fmin"] = cli_dict["SD-CWT"]["freq_min"]
        settings_dict["spectral_widget"]["cwt_fmax"] = cli_dict["SD-CWT"]["freq_max"]
        settings_dict["spectral_widget"]["cwt_adapt_win_len"] = cli_dict["SD-CWT"]["window_len"]
        # cli_dict["SD-CWT"]["window_step"] GUI automatically calculates window step from window length 
        settings_dict["spectral_widget"]["cwt_pval"] = cli_dict["SD-CWT"]["p_value"]
        settings_dict["spectral_widget"]["cwt_cluster_freq_scale"] = cli_dict["SD-CWT"]["freq_tm_factor"]
        settings_dict["spectral_widget"]["cwt_cluster_eps"] = cli_dict["SD-CWT"]["cluster_eps"]
        settings_dict["spectral_widget"]["cwt_cluster_min_samples"] = cli_dict["SD-CWT"]["cluster_min_samples"]

        # ASSOC node
        settings_dict["location_widget"]["bisl_dict"]["bisl_bm_width"] = cli_dict["ASSOC"]["back_az_width"]
        settings_dict["location_widget"]["bisl_dict"]["bisl_rng_max"] = cli_dict["ASSOC"]["range_max"]
        settings_dict["location_widget"]["bisl_dict"]["bisl_resolution"] = cli_dict["ASSOC"]["resolution"]
        # cli_dict["ASSOC"]["cpu_cnt"] GUI doesn't save cpu cnt

        # LOC node
        settings_dict["location_widget"]["bisl_dict"]["bisl_bm_width"] = cli_dict["LOC"]["back_az_width"]
        settings_dict["location_widget"]["bisl_dict"]["bisl_rng_max"] = cli_dict["LOC"]["range_max"]

        # GUI node  -- gui specific settings go here
        settings_dict['waveforms_widget']['filter_dict']['apply_filter'] = cli_dict['GUI']['wf_filter_apply']
        settings_dict['waveforms_widget']['filter_dict']['filter_type'] = cli_dict['GUI']['wf_filter_type']
        settings_dict['waveforms_widget']['filter_dict']['order'] = cli_dict['GUI']['wf_filter_order']
        settings_dict['waveforms_widget']['filter_dict']['zero_phase'] = cli_dict['GUI']['wf_filter_zerophase']
        settings_dict['waveforms_widget']['filter_dict']['show_unfiltered'] = cli_dict['GUI']['wf_filter_show_unfiltered']
        settings_dict['waveforms_widget']['psd_dict']['fft_n'] = cli_dict['GUI']['wf_psd_fftn']
        settings_dict['waveforms_widget']['psd_dict']['window'] = cli_dict['GUI']['wf_psd_window']

        settings_dict['beamforming_widget']['gui_bf_colormap'] = cli_dict["GUI"]["bf_colormap"]
        settings_dict['location_widget']['gui_borders'] = cli_dict["GUI"]["loc_borders"]
        settings_dict['location_widget']['gui_states'] = cli_dict["GUI"]["loc_states"]
        settings_dict['location_widget']['gui_lakes'] = cli_dict["GUI"]["loc_lakes"]
        settings_dict['location_widget']['gui_rivers'] = cli_dict["GUI"]["loc_rivers"]
        settings_dict['location_widget']['gui_coastline'] = cli_dict["GUI"]["loc_coastline"]

        settings_dict['location_widget']['gui_ocean_color'] = cli_dict['GUI']['loc_ocean_color']
        settings_dict['location_widget']['gui_land_color'] = cli_dict['GUI']['loc_land_color'] 
        settings_dict['location_widget']['gui_resolution'] = cli_dict['GUI']['loc_resolution']
        settings_dict['location_widget']['gui_use_background'] = cli_dict['GUI']['loc_use_background']
        settings_dict['location_widget']['gui_offline'] = cli_dict['GUI']['loc_offline']
        settings_dict['location_widget']['gui_offline_dir'] = cli_dict['GUI']['loc_offline_map_dir']
        settings_dict['location_widget']['extent_dict']['ll_lat'] = cli_dict['GUI']['loc_ext_ll_lat']
        settings_dict['location_widget']['extent_dict']['ll_lon'] = cli_dict['GUI']['loc_ext_ll_lon']
        settings_dict['location_widget']['extent_dict']['ur_lat'] = cli_dict['GUI']['loc_ext_ur_lat']
        settings_dict['location_widget']['extent_dict']['ur_lon'] = cli_dict['GUI']['loc_ext_ur_lon']

        settings_dict['spectral_widget']['gui_colormap'] = cli_dict['GUI']['spec_colormap']
        settings_dict['spectral_widget']['gui_colorbar'] = cli_dict['GUI']['spec_colorbar']

        return settings_dict

    def map_infraview_settings_to_cli(self, settings_dict):
        # Convert the settings widget dictionary to one that can be read into the command line config ini file.

        # load infrapy defaults, then if there's something we don't write over, it will automatically have the default value
        if self.default_config_path.exists():
            # print("Reading: {}".format(str(self.default_config_path)))
            cli_template_dict = self.ini_to_dict(str(self.default_config_path))
        else:
            raise FileNotFoundError

        # FK node
        cli_template_dict["FK"]["freq_min"] = settings_dict["waveforms_widget"]["filter_dict"]["highpass"]
        cli_template_dict["FK"]["freq_max"] = settings_dict["waveforms_widget"]["filter_dict"]["lowpass"]
        # cli_template_dict["FK"]["signal_start"] = settings_dict["beamforming_widget"]["signal_start"]
        # cli_template_dict["FK"]["signal_end"] = settings_dict["beamforming_widget"]["signal_end"]
        # cli_template_dict["FK"]["noise_start"] = settings_dict["beamforming_widget"]["noise_start"]
        # cli_template_dict["FK"]["noise_end"] = settings_dict["beamforming_widget"]["noise_end"]
        cli_template_dict["FK"]["back_az_min"] = settings_dict["beamforming_widget"]["backAz_start"]
        cli_template_dict["FK"]["back_az_max"] = settings_dict["beamforming_widget"]["backAz_end"]
        cli_template_dict["FK"]["back_az_step"] = settings_dict["beamforming_widget"]["backAz_resolution"]
        cli_template_dict["FK"]["trace_vel_min"] = settings_dict["beamforming_widget"]["traceV_min"]
        cli_template_dict["FK"]["trace_vel_max"] = settings_dict["beamforming_widget"]["traceV_max"]
        cli_template_dict["FK"]["trace_vel_step"] = settings_dict["beamforming_widget"]["traceV_resolution"]
        cli_template_dict["FK"]["method"] = settings_dict["beamforming_widget"]["method"]

        cli_template_dict["FK"]["window_len"] = settings_dict["beamforming_widget"]["win_length"]
        # cli_template_dict["FK"]["sub_window_len"] = settings_dict["beamforming_widget"]["sub_win_length"] not currently in gui
        cli_template_dict["FK"]["window_step"] = settings_dict["beamforming_widget"]["win_step"]
        cli_template_dict["FK"]["cpu_cnt"] = None  # # GUI doesn't save cpu cnt

        # FD node
        # cli_template_dict["FD"]["window_len"] = None  # Don't write the FD window_len from the gui at this point
        cli_template_dict["FD"]["p_value"] = settings_dict["beamforming_widget"]["detector_settings"]["pval"]
        cli_template_dict["FD"]["min_duration"] = settings_dict["beamforming_widget"]["detector_settings"]["min_peak_width"]
        cli_template_dict["FD"]["back_az_width"] = settings_dict["beamforming_widget"]["detector_settings"]["back_az_limit"]
        cli_template_dict["FD"]["fixed_thresh"] = settings_dict["beamforming_widget"]["detector_settings"]["manual_level"]
        cli_template_dict["FD"]["thresh_ceil"] = None  # Default
        cli_template_dict["FD"]["merge_dets"] = settings_dict["beamforming_widget"]["detector_settings"]["merge"]

        # SD node
        # Spectrogram v. spectrogram - check string sanitization
        cli_template_dict["SD"]["spectral_option"] = settings_dict["spectral_widget"]["spec_type"]
        cli_template_dict["SD"]["morlet_omega0"] = settings_dict["spectral_widget"]["omega0"]
        cli_template_dict["SD"]["freq_min"] = settings_dict["spectral_widget"]["fmin"]
        cli_template_dict["SD"]["freq_max"] = settings_dict["spectral_widget"]["fmax"]
        cli_template_dict["SD"]["window_len"] = settings_dict["spectral_widget"]["adapt_win_len"]               # 900 default
        cli_template_dict["SD"]["window_step"] = settings_dict["spectral_widget"]["adapt_win_len"] / 2.         # hard coding for now
        cli_template_dict["SD"]["p_value"] = settings_dict["spectral_widget"]["pval"]
        cli_template_dict["SD"]["freq_tm_factor"] = settings_dict["spectral_widget"]["cwt_cluster_freq_scale"]
        cli_template_dict["SD"]["cluster_eps"] = settings_dict["spectral_widget"]["cluster_eps"]
        cli_template_dict["SD"]["cluster_min_samples"] = settings_dict["spectral_widget"]["cluster_min_samples"]

        cli_template_dict["SD-CWT"]["freq_min"] = settings_dict["spectral_widget"]["cwt_fmin"]
        cli_template_dict["SD-CWT"]["freq_max"] = settings_dict["spectral_widget"]["cwt_fmax"]
        cli_template_dict["SD-CWT"]["window_len"] = settings_dict["spectral_widget"]["cwt_adapt_win_len"]       # 900 default
        cli_template_dict["SD-CWT"]["window_step"] = settings_dict["spectral_widget"]["adapt_win_len"] / 2.     # hard coding for now
        cli_template_dict["SD-CWT"]["p_value"] = settings_dict["spectral_widget"]["cwt_pval"]
        cli_template_dict["SD-CWT"]["freq_tm_factor"] = settings_dict["spectral_widget"]["cwt_cluster_freq_scale"]
        cli_template_dict["SD-CWT"]["cluster_eps"] = settings_dict["spectral_widget"]["cluster_eps"]
        cli_template_dict["SD-CWT"]["cluster_min_samples"] = settings_dict["spectral_widget"]["cluster_min_samples"]

        # ASSOC node
        cli_template_dict["ASSOC"]["back_az_width"] = settings_dict["location_widget"]["bisl_dict"]["bisl_bm_width"]
        cli_template_dict["ASSOC"]["range_max"] = settings_dict["location_widget"]["bisl_dict"]["bisl_rng_max"]
        cli_template_dict["ASSOC"]["resolution"] = settings_dict["location_widget"]["bisl_dict"]["bisl_resolution"]
        cli_template_dict["ASSOC"]["cpu_cnt"] = None  # GUI doesn't save cpu cnt

        # LOC node
        cli_template_dict["LOC"]["back_az_width"] = settings_dict["location_widget"]["bisl_dict"]["bisl_bm_width"]
        cli_template_dict["LOC"]["range_max"] = settings_dict["location_widget"]["bisl_dict"]["bisl_rng_max"]
        # cli_template_dict["LOC"]["latlon_resol"] = 
        # cli_template_dict["LOC"]["tm_resol"] =
        # cli_template_dict["LOC"]["src_est"] = 
        # cli_template_dict["pgm_model"] = 

        #GUI node -- GUI specific settings go here
        cli_template_dict['GUI']['wf_filter_apply'] = settings_dict['waveforms_widget']['filter_dict']['apply_filter']
        cli_template_dict['GUI']['wf_filter_type'] = settings_dict['waveforms_widget']['filter_dict']['filter_type']
        cli_template_dict['GUI']['wf_filter_order'] = settings_dict['waveforms_widget']['filter_dict']['order']
        cli_template_dict['GUI']['wf_filter_zerophase'] = settings_dict['waveforms_widget']['filter_dict']['zero_phase']
        cli_template_dict['GUI']['wf_filter_show_unfiltered'] = settings_dict['waveforms_widget']['filter_dict']['show_unfiltered']
        cli_template_dict['GUI']['wf_psd_fftn'] = settings_dict['waveforms_widget']['psd_dict']['fft_n']
        cli_template_dict['GUI']['wf_psd_window'] =settings_dict['waveforms_widget']['psd_dict']['window']

        cli_template_dict['GUI']["bf_colormap"] = settings_dict['beamforming_widget']['gui_bf_colormap']
        
        cli_template_dict['GUI']['loc_borders'] = settings_dict['location_widget']['gui_borders']
        cli_template_dict['GUI']['loc_states'] = settings_dict['location_widget']['gui_states']
        cli_template_dict['GUI']['loc_lakes'] = settings_dict['location_widget']['gui_lakes']
        cli_template_dict['GUI']['loc_rivers'] = settings_dict['location_widget']['gui_rivers']
        cli_template_dict['GUI']['loc_coastline'] = settings_dict['location_widget']['gui_coastline']
        
        cli_template_dict['GUI']['loc_ocean_color'] = settings_dict['location_widget']['gui_ocean_color']
        cli_template_dict['GUI']['loc_land_color'] = settings_dict['location_widget']['gui_land_color']
        cli_template_dict['GUI']['loc_resolution'] = settings_dict['location_widget']['gui_resolution']
        cli_template_dict['GUI']['loc_use_background'] = settings_dict['location_widget']['gui_use_background']
        cli_template_dict['GUI']['loc_offline'] = settings_dict['location_widget']['gui_offline']
        cli_template_dict['GUI']['loc_offline_map_dir'] = settings_dict['location_widget']['gui_offline_dir']

        cli_template_dict['GUI']['loc_ext_ll_lat'] = settings_dict['location_widget']['extent_dict']['ll_lat']
        cli_template_dict['GUI']['loc_ext_ll_lon'] = settings_dict['location_widget']['extent_dict']['ll_lon']
        cli_template_dict['GUI']['loc_ext_ur_lat'] = settings_dict['location_widget']['extent_dict']['ur_lat']
        cli_template_dict['GUI']['loc_ext_ur_lon'] = settings_dict['location_widget']['extent_dict']['ur_lon']

        cli_template_dict['GUI']['spec_colormap'] = settings_dict['spectral_widget']['gui_colormap']
        cli_template_dict['GUI']['spec_colorbar'] = settings_dict['spectral_widget']['gui_colorbar']

        return cli_template_dict