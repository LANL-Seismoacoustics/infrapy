import pyqtgraph as pg
import numpy as np
import os, platform, yaml
import qdarktheme, darkdetect
from pathlib import Path


# obspy includes
from obspy.core import read as obsRead
from obspy.core import UTCDateTime
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.core.stream import Stream

# PyQt5 includes
from PyQt5 import QtWidgets
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt, pyqtSignal, pyqtSlot, QSettings, QSize, QPoint, QDir
from PyQt5.QtWidgets import (QDialog, QFileDialog, QTabWidget, QStackedWidget,
                             QFormLayout, QHBoxLayout, QVBoxLayout, QLabel, QWidget,
                             QDialog, QDialogButtonBox, QDoubleSpinBox, QLineEdit)

# Infrapy includes

# Application includes
from InfraView.widgets import IPBeamformingWidget
from InfraView.widgets import IPDatabaseWidget
from InfraView.widgets import IPFDSNDialog
from InfraView.widgets import IPLocationWidget
from InfraView.widgets import IPMainMenu
from InfraView.widgets import IPProject
from InfraView.widgets import IPSaveAllDialog
from InfraView.widgets import IPSettingsManager
from InfraView.widgets import IPSingleSensorWidget
from InfraView.widgets import IPUtils
from InfraView.widgets import IPWaveformWidget


# multiprocessing modules
import pathos.multiprocessing as mp
from multiprocessing import cpu_count

# for tracking down warnings...
import warnings


class IPApplicationWindow(QtWidgets.QMainWindow):

    warnings.filterwarnings('ignore', message='Item already added to PlotItem, ignoring')
    
    sig_stream_changed = pyqtSignal(Stream)
    sig_inventory_changed = pyqtSignal(Inventory, str)

    # this tab will emit the name of the widget corresponding to the clicked action 
    sig_widget_changed = pyqtSignal(str)


    # variable to hold the reference of the loaded project object (if any)
    project = None

    # we will have one multiprocessing pool used by any/all widgets
    mp_pool = None

    # application settings
    is_fullscreen = False

    def __init__(self, qApp, progname, progversion):
        super().__init__()

        self.ipApp = qApp   # reference to the application

        pg.setConfigOption('antialias', True)
        pg.setConfigOption('foreground', pg.mkColor(128,128,128))

        self.progname = progname
        self.progversion = progversion

        # initialize the multiproccessor pool
        self.mp_pool = mp.ProcessingPool(cpu_count() - 1)

        font = self.font()
        font.setFamily('monospace')

        self.buildUI()

        # Initialize theme based on what's checked.   Sometimes works
        if platform.system() == 'Darwin':
            self.set_theme(darkdetect.theme().lower())
        else:
            if self.menuBar.light_action.isChecked():
                self.set_theme('light')
            elif self.menuBar.dark_action.isChecked():
                self.set_theme('dark')

    def buildUI(self):
        self.setWindowTitle("Infraview")
        self.main_widget = QWidget(self)

        # Create the settings manager
        # The settings manager creates the settings widgets
        self.settings_manager = IPSettingsManager.IPSettingsManager(self)

        # Create the main widgets
        self.beamformingWidget = IPBeamformingWidget.IPBeamformingWidget(self, self.mp_pool)
        self.singleSensorWidget = IPSingleSensorWidget.IPSingleSensorWidget(self, self.mp_pool)
        self.waveformWidget = IPWaveformWidget.IPWaveformWidget(self, self.mp_pool)
        self.locationWidget = IPLocationWidget.IPLocationWidget(self, self.mp_pool)
        self.databaseWidget = IPDatabaseWidget.IPDatabaseWidget(self)

        # For convenience, keep referencces to the widgets in a dict. 
        # to make things simple, have the keys be lower case strings of the Tab Names
        self.widget_dict = {'waveforms': self.waveformWidget,
                            'beamforming': self.beamformingWidget,
                            'location': self.locationWidget,
                            'database': self.databaseWidget,
                            'spectral': self.singleSensorWidget,
                            'app_window': self}
        
        # now that we have a settings manager and settings widgets along with the actual controled widgets
        # we need to give a reference of the controlled widgets to the respective settings widgets
        self.settings_manager.connect_widgets_and_settings(self.widget_dict)

        # Add widgets to the main stack
        self.mainStack = QStackedWidget()
        self.mainStack.addWidget(self.widget_dict['waveforms'])
        self.mainStack.addWidget(self.widget_dict['beamforming'])
        self.mainStack.addWidget(self.widget_dict['location'])
        self.mainStack.addWidget(self.widget_dict['database'])
        self.mainStack.addWidget(self.widget_dict['spectral'])

        # Put the settings above the tabs
        mainLayout = QVBoxLayout(self.main_widget)
        mainLayout.addWidget(self.settings_manager)        # add the main widgets to the application tabs
        mainLayout.addWidget(self.mainStack)
        mainLayout.setContentsMargins(0,0,0,0)

        # Load the default settings
        self.settings_manager.load_default_settings()

        # All menu items should be located in makeMenuBar method
        self.menuBar = IPMainMenu.IPMainMenuBar(self, self.widget_dict)
        self.setMenuBar(self.menuBar)
        self.menuBar.activate_waveforms(True)

        # Create Dialogs
        self.fdsnDialog = IPFDSNDialog.IPFDSNDialog(self)
        self.saveAllDialog = IPSaveAllDialog.IPSaveAllDialog(self)

        self.restoreWindowGeometrySettings()
        self.connectSignalsAndSlots()

        self.setCentralWidget(self.main_widget)
        self.main_widget.setContentsMargins(0,0,0,0)

        # Instantiate dialogs here
        self.fill_sta_info_dialog = IPFillStationInfoDialog()
        self.redundant_trace_dialog = IPRedundantTraceDialog()
        self.aboutDialog = IPAboutDialog(self.progname, self.progversion)

    def debug_trace(self):  # for debugging, you have to call pyqtRemoveInputHook before set_trace()
        from PyQt5.QtCore import pyqtRemoveInputHook
        from pdb import set_trace
        pyqtRemoveInputHook()
        set_trace()

    @pyqtSlot()
    def on_palette_change(self):
        # OSX is always auto set, if Linux or Windows, check to see if the theme is auto changed
        if platform.system() == 'Darwin':
            self.set_theme(darkdetect.theme().lower())
            return
        # For linux and windows, check to see if auto is selected
        if self.menuBar.auto_action.isChecked():
            if darkdetect.isDark():
                self.set_theme('dark')
            elif darkdetect.isLight():
                self.set_theme('light')
    
    @pyqtSlot(str)
    def set_theme(self, t):
        if t == 'auto':
            t = darkdetect.theme().lower()

        # linux and windows acknowledge the pyqtdark themes, mac does its own thing
        if platform.system() == 'Linux' or platform.system() == 'Windows':
            qdarktheme.setup_theme(t, corner_shape='sharp')
        # Still need to update the pyqtgraph backgrounds and other similar things for all OSes...
        self.widget_dict['waveforms'].update_theme(t)
        self.widget_dict['beamforming'].update_theme(t)
        self.widget_dict['spectral'].update_theme(t)
        self.widget_dict['location'].update_theme(t)

    @pyqtSlot()
    def toggle_settings(self):
         self.settings_manager.setVisible(self.settings_manager.isHidden())

    @pyqtSlot(str)
    def activate_widget(self, name):
        self.mainStack.setCurrentWidget(self.widget_dict[name])
        self.setWindowTitle(self.progname + " - " + name.title())
        self.sig_widget_changed.emit(name)

    def connectSignalsAndSlots(self):
        self.menuBar.sig_set_theme.connect(self.set_theme)

        self.ipApp.paletteChanged.connect(self.on_palette_change)

        self.fdsnDialog.fdsnWidget.sigTracesAppended.connect(self.waveformWidget.appendTraces)
        self.fdsnDialog.fdsnWidget.sigTracesReplaced.connect(self.waveformWidget.replaceTraces)

        self.sig_stream_changed.connect(self.waveformWidget.update_streams)
        self.sig_inventory_changed.connect(self.waveformWidget.stationViewer.merge_new_inventory)

        self.beamformingWidget.detectionWidget.signal_detections_changed.connect(self.locationWidget.update_detections)
        self.beamformingWidget.detectionWidget.signal_detections_cleared.connect(self.locationWidget.detections_cleared)

        self.databaseWidget.ipdatabase_query_results_table.signal_new_stream_from_db.connect(self.database_add_streams)
        self.databaseWidget.ipev_query_results_table.sig_origin_changed.connect(self.locationWidget.showgroundtruth.ev_widget.setEvent)

        # connect tabs to the settings manager so it can adjust when tabs are clicked
        self.sig_widget_changed.connect(self.settings_manager.widget_changed)

        self.waveformWidget.plotViewer.lr_settings_widget.noiseSpinsChanged.connect(self.beamformingWidget.bottomSettings.setNoiseValues)
        self.waveformWidget.plotViewer.lr_settings_widget.signalSpinsChanged.connect(self.beamformingWidget.bottomSettings.setSignalValues)
        self.waveformWidget.psdWidget.f1_Spin.valueChanged.connect(self.beamformingWidget.bottomSettings.setFmin)
        self.waveformWidget.psdWidget.f2_Spin.valueChanged.connect(self.beamformingWidget.bottomSettings.setFmax)
        self.waveformWidget.psdWidget.psdPlot.getFreqRegion().sigRegionChanged.connect(self.beamformingWidget.bottomSettings.setFreqValues)

        self.menuBar.sig_activate_widget.connect(self.activate_widget)
        

    def setStatus(self, s, ms=0):
        self.statusBar().showMessage(s, ms)

    def filemenu_NewProject(self):
        self.project = IPProject.IPProject()
        if self.project.makeNewProject():
            self.setWindowTitle(self.progname + ' - ' + self.project.getprojectName())
            self.setStatus(self.project.get_projectName() + ' successfully loaded', 5000)

            settings = QSettings('LANL', 'InfraView')
            settings.beginGroup('General')
            settings.setValue("last_baseProject_directory", str(self.project.get_basePath()))
            settings.setValue('last_project_directory', str(self.project.get_projectPath))
            settings.endGroup()
            
    def filemenu_LoadProject(self):
        newProject = IPProject.IPProject()
        if newProject.loadProject():
            self.project = newProject
            self.setWindowTitle(self.progname + ' - ' + self.project.get_projectName())
            self.setStatus(self.project.get_projectName() + ' successfully loaded', 5000)

            settings = QSettings('LANL', 'InfraView')
            settings.beginGroup('General')
            settings.setValue("last_baseProject_directory", str(self.project.get_basePath()))
            settings.setValue('last_project_directory', str(self.project.get_projectPath))
            settings.endGroup()

    def filemenu_CloseProject(self):
        self.setWindowTitle(self.progname)
        # self.consoleBox.append('Closed Project: ' + self.project.get_projectName())
        self.filemenu_ClearWaveforms()
        if self.project is not None:
            self.setStatus(self.project.get_projectName() + ' successfully closed', 5000)
            self.project = None

    def getProject(self):
        return self.project

    def get_earliest_start_time(self):
        return self._earliest_start_time

    def get_latest_end_time(self):
        return self._latest_end_time

    def filemenu_Quit(self):
        self.close()

    def filemenu_Open(self):
        # if this hasn't been run yet, I want to default to opening the /data directory.  This can be found relative to the 
        # current directory via...
        default_data_dir = os.path.join(os.path.dirname(__file__), '../../examples/data')
        if self.project is None:
            # force a new filename...
            settings = QSettings('LANL', 'InfraView')
            settings.beginGroup('General')
            previous_directory = settings.value("last_open_directory", default_data_dir)
            settings.endGroup()
        else:
            # There is an open project, so make the default save location
            # correspond to what the project wants
            previous_directory = str(self.project.get_dataPath())

        ifiles = QFileDialog.getOpenFileNames(self, 'Open File...', previous_directory)

        new_inventory = None
        new_stream = None

        if len(ifiles[0]) > 0:

            for ifile in ifiles[0]:
                current_trace_names = []
                if self.waveformWidget._sts is not None:
                    for trace in self.waveformWidget._sts:
                        current_trace_names.append(trace.id)

                ipath = os.path.dirname(ifile)
                if self.project is None:
                    settings = QSettings('LANL', 'InfraView')
                    settings.beginGroup('General')
                    settings.setValue("last_open_directory", ipath)
                    settings.endGroup()
                else:
                    self.project.set_dataPath(ipath)
                try:

                    new_stream = obsRead(ifile)
                    
                    for trace in new_stream:
                        
                        trace_name = trace.id

                        if trace_name in current_trace_names:
                            # redundant trace!
                            _, staid, _, _ = self.parseTraceName(trace_name)
                            self.redundant_trace_dialog.exec_(trace_name)

                            if  self.redundant_trace_dialog.get_result():
                                # if accepted, they want to use the new trace so first remove the old one
                                self.waveformWidget.remove_from_inventory(staid)
                                
                            else:
                                # if rejected, they want to keep the old trace, and ignore this one
                                #so remove the trace from new_stream, and continue to the next trace
                                new_stream.remove(trace)
                                continue

                        # do our best to generate new inventory from the new stream
                        if new_inventory is None:
                            new_inventory = self.trace_to_inventory(trace)
                        else:
                            new_inventory += self.trace_to_inventory(trace)

                        # for now we will remove dc offset when loading the file.  Maybe should be an option?
                        trace.data = trace.data - np.mean(trace.data)
                    
                    if self.waveformWidget._sts is not None:   
                        self.waveformWidget._sts += new_stream
                    else:
                        self.waveformWidget._sts = new_stream
                except Exception as e:
                    print("Unexpected error: ", str(e))
                    self.setStatus("Unexpected error: ", 5000)
                    continue

                self.waveformWidget._sts.merge(fill_value=0)
        else:
            # No files were chosen to open
            return

        # it's possible that self.waveformWidget._sts is still None, so if it is, bail out
        # if not populate the trace stats viewer and plot the traces

        if self.waveformWidget._sts is not None:
            self.beamformingWidget.setStreams(self.waveformWidget._sts)
            self.sig_stream_changed.emit(self.waveformWidget._sts)

            if new_inventory is not None:
                
                self.sig_inventory_changed.emit(new_inventory, 'PROMPT')

            # self.mainTabs.setCurrentIndex(0)
            self.mainStack.setCurrentWidget(self.widget_dict['waveforms'])
            self.menuBar.toggle_enable('waveforms')

        else:
            return

    def filemenu_import(self):
        if self.fdsnDialog.exec_():
            self.menuBar.activate_waveforms(True)

    def filemenu_ClearWaveforms(self):
        self.beamformingWidget.clearWaveformPlot()
        self.waveformWidget.clearWaveforms()
        self.singleSensorWidget.clearWaveformPlots()
        self.waveformWidget.stationViewer.clear_all()

    def filemenu_saveAllWaveforms(self):
        if self.waveformWidget._sts is None:
            IPUtils.errorPopup('Oh dear...no waveforms to save.')
            return

        if self.project is None:
            # force a new filename...
            settings = QSettings('LANL', 'InfraView')
            settings.beginGroup('General')
            previousDirectory = settings.value("last_open_directory", QDir.homePath())
            settings.endGroup()
        else:
            # There is an open project, so make the default save location
            # correspond to what the project wants
            previousDirectory = str(self.project.get_dataPath())

        accepted = self.saveAllDialog.exec_(self.waveformWidget._sts, previousDirectory)

        if accepted:

            saveDir = self.saveAllDialog.getSaveDirectory()
            whichData = self.saveAllDialog.getFileChoiceData()
            if whichData == 1:
                fileInfo = self.saveAllDialog.getFileInfo()
                for idx, trace in enumerate(self.waveformWidget.get_streams()):
                    trace.write(saveDir + '/' + fileInfo[idx]['fname'], format=fileInfo[idx]['format'])
            elif whichData == 2:
                # filter traces, then save them
                filteredfileInfo = self.saveAllDialog.getFilteredFileInfo()
                for idx, trace in enumerate(self.waveformWidget.get_filtered_streams()):
                    trace.write(saveDir + '/' + filteredfileInfo[idx]['fname'], format=filteredfileInfo[idx]['format'])

            elif whichData == 3:
                fileInfo = self.saveAllDialog.getFileInfo()
                filteredfileInfo = self.saveAllDialog.getFilteredFileInfo()
                for idx, trace in enumerate(self.waveformWidget.get_streams()):
                    trace.write(
                        saveDir + '/' + fileInfo[idx]['fname'], format=fileInfo[idx]['format'])
                # filter traces, then save them
                for idx, trace in enumerate(self.waveformWidget.get_filtered_streams()):
                    trace.write(saveDir + '/' + filteredfileInfo[idx]['fname'], format=filteredfileInfo[idx]['format'])

        return

    @pyqtSlot(Stream, bool)
    def database_add_streams(self, new_stream, append):
        current_trace_names = []
        new_inventory = None
        if self.waveformWidget._sts:
            for trace in self.waveformWidget._sts:
                current_trace_names.append(trace.id)
        for trace in new_stream:
            trace_name = trace.id
            if trace_name in current_trace_names:
                # redundant trace!
                _, staid, _, _ = self.parseTraceName(trace_name)
                self.redundant_trace_dialog.exec_(trace_name)

                if  self.redundant_trace_dialog.get_result():
                    # if accepted, they want to use the new trace so first remove the old one
                    self.waveformWidget.remove_from_inventory(staid)
                else:
                    # if rejected, they want to keep the old trace, and ignore this one
                    #so remove the trace from new_stream, and continue to the next trace
                    new_stream.remove(trace)
                    continue
            # do our best to generate new inventory from the new stream
            if new_inventory is None:
                new_inventory = self.trace_to_inventory(trace)
            else:
                new_inventory += self.trace_to_inventory(trace)

            # for now we will remove dc offset when loading the file.  Maybe should be an option?
            trace.data = trace.data - np.mean(trace.data)
        
        if append:
            self.waveformWidget.appendTraces(new_stream, new_inventory)
        else: # Replace
            self.waveformWidget.replaceTraces(new_stream, new_inventory)

    def parseTraceName(self, trace_name):
        bits = trace_name.split('.')
        return bits[0], bits[1], bits[2], bits[3]

    def viewmenu_toggle_fullscreen(self):
        if self.is_fullscreen:
            self.showNormal()
            self.is_fullscreen = False
        else:
            self.showFullScreen()
            self.is_fullscreen = True
                
    def trace_to_inventory(self, trace):
        # if sac files are opened, it's useful to extract inventory from their streams so that we can populate the 
        # stations tabs and the location widget
        new_inventory = None

        # The next bit is modified from the obspy webpage on building a stationxml site from scratch
        # https://docs.obspy.org/tutorial/code_snippets/stationxml_file_from_scratch.html
        #
        # We'll first create all the various objects. These strongly follow the
        # hierarchy of StationXML files.
        # initialize the lat/lon/ele 
        lat = 0.0
        lon = 0.0
        ele = 0.0

        network = trace.stats['network'] 
        station = trace.stats['station']
        channel = trace.stats['channel']
        location = trace.stats['location']

        # if the trace is from a sac file, the sac header might have some inventory information
        if trace.stats['_format'] == 'SAC':

            if 'stla' in trace.stats['sac']:
                lat = trace.stats['sac']['stla']

            if 'stlo' in trace.stats['sac']:
                lon = trace.stats['sac']['stlo']
        
            if 'stel' in trace.stats['sac']:
                ele = trace.stats['sac']['stel']
            else:
                ele = 0.0

        if lat == 0.0 or lon == 0.0 or ele < 0:
            if self.fill_sta_info_dialog.exec_(network, station, location, channel, lat, lon, ele):
                
                edited_values = self.fill_sta_info_dialog.get_values()
                
                lat = edited_values['lat']
                lon = edited_values['lon']
                ele = edited_values['ele']

                network = edited_values['net'] 
                station = edited_values['sta']
                location = edited_values['loc']
                channel = edited_values['cha']

                # (re)populate sac headers where possible
                if trace.stats['_format'] == 'SAC':
                    trace.stats['sac']['stla'] = lat
                    trace.stats['sac']['stlo'] = lon
                    trace.stats['sac']['stel'] = ele
                    trace.stats['sac']['knetwk'] = network
                    trace.stats['sac']['kstnm'] = station
                # (re)populate trace stats where possible
                trace.stats['network'] = network
                trace.stats['station'] = station
                trace.stats['location'] = location
                trace.stats['channel'] = channel
        try:            
            new_inventory = Inventory(
                # We'll add networks later.
                networks=[],
                # The source should be the id whoever create the file.
                source="InfraView")
            
            net = Network(
                # This is the network code according to the SEED standard.
                code=network,
                # A list of stations. We'll add one later.
                stations=[],
                # Description isn't something that's in the trace stats or SAC header, so lets set it to the network cod
                description=network,
                # Start-and end dates are optional.

                # Start and end dates for the network are not stored in the sac header so lets set it to 1/1/1900
                start_date=UTCDateTime(1900, 1, 1))

            sta = Station(
                # This is the station code according to the SEED standard.
                code=station,
                latitude=lat,
                longitude=lon,
                elevation=ele,
                # Creation_date is not saved in the trace stats or sac header
                creation_date=UTCDateTime(1900, 1, 1),
                # Site name is not in the trace stats or sac header, so set it to the site code
                site=Site(name=station))

            # This is the channel code according to the SEED standard.
            cha = Channel(code=channel,
                        # This is the location code according to the SEED standard.
                        location_code=location,
                        # Note that these coordinates can differ from the station coordinates.
                        latitude=lat,
                        longitude=lon,
                        elevation=ele,
                        depth=0.0)

            # Now tie it all together.
            # cha.response = response
            sta.channels.append(cha)
            net.stations.append(sta)
            new_inventory.networks.append(net)

            return new_inventory

        except ValueError:
            bad_values = ""
            if lon < -180 or lon > 180:
                bad_values = bad_values + "\tlon = " + str(lon) + "\n"
            if lat < -90 or lat > 90:
                bad_values = bad_values + "\tlat = " + str(lat)
            IPUtils.errorPopup("There seems to be a value error in "+ network + "." + station + "." + channel + "\nPossible bad value(s) are:\n" + bad_values)

    # ------------------------------------------------------------------------------
    # Settings methods

    def restoreWindowGeometrySettings(self):
        # Restore the widgets geometry settings
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('MainWindow')
        self.resize(settings.value("windowSize", QSize(1000, 900)))
        self.move(settings.value("windowPos", QPoint(200, 200)))
        settings.endGroup()

        self.beamformingWidget.restoreWindowGeometrySettings()
        self.locationWidget.restoreWindowGeometrySettings()
        self.waveformWidget.restoreWindowGeometrySettings()

    def saveWindowGeometrySettings(self):
        # save the widgets geometry settings
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('MainWindow')
        settings.setValue("windowSize", self.size())
        settings.setValue("windowPos", self.pos())
        settings.endGroup()

        self.beamformingWidget.saveWindowGeometrySettings()
        self.locationWidget.saveWindowGeometrySettings()
        self.waveformWidget.saveWindowGeometrySettings()

    # -------------------------------------------------------
    # Clean up

    def closeEvent(self, event):
        self.saveWindowGeometrySettings()
        self.mp_pool.close()
        event.accept()

    # -------------------------------------------------------
    # Obligatory about

    def helpmenu_about(self):
        self.aboutDialog.exec_()
        #QtWidgets.QMessageBox.about(self, "About", self.progname + "   " + self.progversion + "\n" + "Copyright 2018\nLos Alamos National Laboratory")


class IPAboutDialog(QDialog):
    def __init__(self, progname, progversion):
        super().__init__()
        self.buildUI(progname, progversion)

    def buildUI(self, name, version):
        self.setWindowTitle("InfraView: About")
        
        info_label = QLabel('\n' + name + '\nVersion: ' + version + '\nCopyright 2018 Los Alamos National Laboratory\n')
        label_font = info_label.font()
        label_font.setPixelSize(12)
        info_label.setFont(label_font)

        image_path = Path(Path(__file__).parent.parent.parent, 'infrapy', 'resources', 'PNG', 'LANL_Logo_Ultramarine.png')
        logo_pixmap = QPixmap(str(image_path))
        logo_label = QLabel(self)
        logo_label.setPixmap(logo_pixmap)

        cpu_label = QLabel('Total CPU Count: ' + str(cpu_count()) + '\nCPUs Used: ' + str(cpu_count() -1))
        cpu_label.setFont(label_font)
        
        layout = QVBoxLayout()
        layout.addWidget(logo_label)
        layout.addWidget(info_label)
        layout.addWidget(cpu_label)
        self.setLayout(layout)


class IPRedundantTraceDialog(QDialog):

    def __init__(self):
        super().__init__()
        self.buildUI()
        
    def exec_(self, trace_name):
        intro_text = "I appears there is already a trace loaded with the name " + trace_name +". Would you like to keep the one that is currently in memory, or replace it with this one?"
        self.intro_label.setText(intro_text)

        return super().exec_()

    def buildUI(self):
        self.setWindowTitle("InfraView: Redundant Trace")
        self.intro_label = QLabel("")
        self.intro_label.setWordWrap(True)

        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
                                   Qt.Horizontal,
                                   self)
        buttons.button(QDialogButtonBox.Ok).setText("Use New One")
        buttons.button(QDialogButtonBox.Cancel).setText("Keep Old One")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        vbox_layout = QVBoxLayout()
        vbox_layout.addWidget(self.intro_label)
        vbox_layout.addWidget(buttons)

        self.setLayout(vbox_layout)

    def accept(self):
        self.result = True
        super().accept()

    def reject(self):
        self.result = False
        super().reject()

    def get_result(self):
        return self.result


class IPFillStationInfoDialog(QDialog):

    def __init__(self):
        super().__init__()
        self.buildUI()
        

    def exec_(self, net, sta, loc, cha, lat, lon, ele):

        self.lat_spin.setValue(lat)
        self.lon_spin.setValue(lon)
        self.ele_spin.setValue(ele)

        self.net_edit.setText(net)
        self.sta_edit.setText(sta)
        self.loc_edit.setText(loc)
        self.cha_edit.setText(cha)
        
        self.update_fullname_label()

        return super().exec_()

    def buildUI(self):
        self.setWindowTitle('InfraView: Trace metadata')
        self.setMinimumWidth(300)
        self.setMaximumWidth(400)

        intro_text = """Some of the inventory information in this trace's stats appears to be absent or incorrect. This could cause problems later. 
        You can edit the information here.  If you have a stationxml file, you can skip this form and load your stationxml file from the Station tab later."""
        self.intro_label = QLabel(intro_text)
        self.intro_label.setWordWrap(True)

        self.lat_spin = QDoubleSpinBox()
        self.lat_spin.setRange(-90.1, 90.0)  # the -90.1 is used as the "unset" value
        self.lat_spin.setDecimals(8)
        self.lat_spin.setMaximumWidth(120)
        self.lat_spin.setSingleStep(0.1)

        self.lon_spin = QDoubleSpinBox()
        self.lon_spin.setRange(-180.1, 180.0)
        self.lon_spin.setDecimals(8)
        self.lon_spin.setMaximumWidth(120)
        self.lon_spin.setSingleStep(0.1)

        self.ele_spin = QDoubleSpinBox()
        self.ele_spin.setRange(-20000,20000)
        self.ele_spin.setMaximumWidth(100)
        self.ele_spin.setDecimals(2)

        self.full_name = QLabel('')

        caps_validator = IPUtils.CapsValidator(self)

        # Network selector
        self.net_edit = QLineEdit()
        self.net_edit.setToolTip('Wildcards OK \nCan be SEED network codes or data center defined codes. \nMultiple codes are comma-separated (e.g. "IU,TA").')
        self.net_edit.setMaximumWidth(60)
        self.net_edit.setValidator(caps_validator)

        self.sta_edit = QLineEdit()
        self.sta_edit.setMaximumWidth(60)
        self.sta_edit.setToolTip('Wildcards OK \nCan be SEED network codes or data center defined codes. \nMultiple codes are comma-separated (e.g. "IU,TA").')
        self.sta_edit.setValidator(caps_validator)

        self.loc_edit = QLineEdit()
        self.loc_edit.setMaximumWidth(60)
        self.loc_edit.setToolTip('Wildcards OK \nCan be SEED network codes or data center defined codes. \nMultiple codes are comma-separated (e.g. "IU,TA").')
        self.loc_edit.setValidator(caps_validator)

        self.cha_edit = QLineEdit()
        self.cha_edit.setMaximumWidth(60)
        self.cha_edit.setToolTip('Wildcards OK \nCan be SEED network codes or data center defined codes. \nMultiple codes are comma-separated (e.g. "IU,TA").')
        self.cha_edit.setValidator(caps_validator)
        
        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
                                   Qt.Horizontal,
                                   self)
        buttons.button(QDialogButtonBox.Ok).setText("Set Values")
        buttons.button(QDialogButtonBox.Cancel).setText("Skip")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        form_layout_col1 = QFormLayout()
        form_layout_col1.addRow("", self.full_name)
        form_layout_col1.addRow("Net: ", self.net_edit)
        form_layout_col1.addRow("Sta: ", self.sta_edit)
        form_layout_col1.addRow("Loc: ", self.loc_edit)
        form_layout_col1.addRow("Cha: ", self.cha_edit)

        form_layout_col2 = QFormLayout()
        form_layout_col2.addRow("Lat: ", self.lat_spin)
        form_layout_col2.addRow("Lon: ", self.lon_spin)
        form_layout_col2.addRow("Ele: ", self.ele_spin)

        hbox = QHBoxLayout()
        hbox.addLayout(form_layout_col1)
        hbox.addLayout(form_layout_col2)

        vbox = QVBoxLayout()
        vbox.addWidget(self.intro_label)
        vbox.addLayout(hbox)
        vbox.addWidget(buttons)

        self.setLayout(vbox)

        self.connectSignalsAndSlots()

    def connectSignalsAndSlots(self):
        self.net_edit.textChanged.connect(self.update_fullname_label)
        self.sta_edit.textChanged.connect(self.update_fullname_label)
        self.loc_edit.textChanged.connect(self.update_fullname_label)
        self.cha_edit.textChanged.connect(self.update_fullname_label)

    def get_values(self):
        values = {"net": self.net_edit.text(), 
                  "sta": self.sta_edit.text(),
                  "loc": self.loc_edit.text(),
                  "cha": self.cha_edit.text(),
                  "lat": self.lat_spin.value(),
                  "lon": self.lon_spin.value(),
                  "ele": self.ele_spin.value()}
        return values

    @pyqtSlot()
    def update_fullname_label(self):
        self.full_name.setText(self.net_edit.text() + '.' + 
                               self.sta_edit.text() + '.' + 
                               self.loc_edit.text() + '.' + 
                               self.cha_edit.text())


