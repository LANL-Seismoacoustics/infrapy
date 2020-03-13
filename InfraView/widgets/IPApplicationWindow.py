import pyqtgraph as pg
import platform
import os
import pdb

import numpy as np
import copy

# obspy includes
import obspy
from obspy.core import read as obsRead
from obspy.core import UTCDateTime
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.core.stream import Stream

# PyQt5 includes
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtGui import QKeySequence
from PyQt5.QtCore import Qt, pyqtSignal, pyqtSlot, QSettings, QSize, QPoint, QDir
from PyQt5.QtWidgets import (QAction, QFileDialog, QTabWidget, QGridLayout, QVBoxLayout,
                             QHBoxLayout, QInputDialog, QLabel, QTableWidgetItem, QMessageBox, QWidget,
                             QApplication, QDialog, QDoubleSpinBox)

# Infrapy includes

# Application includes
from InfraView.widgets import IPPlotWidget
from InfraView.widgets import IPDisplaySettingsWidget
from InfraView.widgets import IPBeamformingWidget
from InfraView.widgets import IPPSDWidget
from InfraView.widgets import IPFilterSettingsWidget
from InfraView.widgets import IPProject
from InfraView.widgets import IPStream
from InfraView.widgets import IPStationView
from InfraView.widgets import IPStatsView
from InfraView.widgets import IPFDSNDialog
from InfraView.widgets import IPLocationWidget
from InfraView.widgets import IPSaveAllDialog
from InfraView.widgets import IPSimpleLegend
from InfraView.widgets import IPWaveformWidget

# multiprocessing modules
import pathos.multiprocessing as mp
from multiprocessing import cpu_count


class IPApplicationWindow(QtWidgets.QMainWindow):

    sig_stream_changed = pyqtSignal(Stream)
    sig_inventory_changed = pyqtSignal(Inventory)

    # variable to hold the reference of the loaded project object (if any)
    _project = None

    # we will have one multiprocessing pool used by any/all widgets
    mp_pool = None

    # application settings
    is_fullscreen = False

    def __init__(self, qApp, progname, progversion):
        super().__init__()

        self.ipApp = qApp   # reference to the application

        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
        pg.setConfigOptions(antialias=False)

        self.settings = QSettings('LANL', 'InfraView')

        self.progname = progname
        self.progversion = progversion

        # initialize the multiproccessor pool
        self.mp_pool = mp.ProcessingPool(cpu_count() - 1)

        self.buildUI()
        self.restoreSettings()

    def buildUI(self):

        self.main_widget = QWidget(self)
        mainLayout = QGridLayout(self.main_widget)

        # All menu items should be located in makeMenuBar method
        self.makeMenuBar()

        # Create the display settings widget
        self.displaySettingsWidget = IPDisplaySettingsWidget.IPDisplaySettingsWidget()

        # Create the main widgets
        self.beamformingWidget = IPBeamformingWidget.IPBeamformingWidget(self, self.mp_pool)
        self.locationWidget = IPLocationWidget.IPLocationWidget(self, self.mp_pool)
        self.waveformWidget = IPWaveformWidget.IPWaveformWidget(self, self.mp_pool)

        # add the main widgets to the application tabs
        self.mainTabs = QTabWidget()
        self.mainTabs.addTab(self.waveformWidget, 'Waveforms')
        self.mainTabs.addTab(self.beamformingWidget, 'Beamforming')
        self.mainTabs.addTab(self.locationWidget, 'Location')

        mainLayout.addWidget(self.mainTabs)

        # Create Dialogs
        self.fdsnDialog = IPFDSNDialog.IPFDSNDialog(self)
        self.saveAllDialog = IPSaveAllDialog.IPSaveAllDialog(self)

        self.restoreSettings()
        self.connectSignalsAndSlots()

        self.setCentralWidget(self.main_widget)

    def debug_trace(self):  # for debugging, you have to call pyqtRemoveInputHook before set_trace()
        from PyQt5.QtCore import pyqtRemoveInputHook
        from pdb import set_trace
        pyqtRemoveInputHook()
        set_trace()

    def makeMenuBar(self):
        main_menu_bar = self.menuBar()

        if platform.system() == 'Darwin':
            main_menu_bar.setNativeMenuBar(False)  # This is because I couldn't get the normal mac menu to work...

        self.file_menu = QtWidgets.QMenu('&File', self)
        self.file_menu.addAction(self.tr(' New Project...'), self.filemenu_NewProject)
        self.file_menu.addAction(self.tr(' Load Project...'), self.filemenu_LoadProject)
        self.file_menu.addAction(self.tr(' Close Project'), self.filemenu_CloseProject)

        self.file_menu.addSeparator()

        self.file_menu.addAction(self.tr(' Load Waveform(s)...'), self.filemenu_Open)
        self.file_menu.addAction(self.tr(' Import Waveform(s)...'), self.filemenu_import)
        self.file_menu.addAction(self.tr(' Clear Waveform(s)'), self.filemenu_ClearWaveforms)
        self.file_menu.addAction(self.tr(' Save Waveform(s)...'), self.filemenu_saveAllWaveforms)

        self.file_menu.addSeparator()

        self.file_menu.addAction(self.tr(' &Exit'), self.filemenu_Quit, QtCore.Qt.CTRL + QtCore.Qt.Key_Q)

        action_show_fullscreen = QAction(self.tr('Toggle Fullscreen'), self)
        action_show_fullscreen.setShortcut(QKeySequence.FullScreen)
        action_show_fullscreen.triggered.connect(self.view_menu_toggle_fullscreen)

        self.view_menu = QtWidgets.QMenu('View', self)
        self.view_menu.addAction(action_show_fullscreen)

        self.help_menu = QtWidgets.QMenu('&Help', self)
        self.help_menu.addAction(self.tr("You're on your own buddy"))
        self.help_menu.addAction('&About', self.about)

        main_menu_bar.addMenu(self.file_menu)
        main_menu_bar.addMenu(self.view_menu)
        main_menu_bar.addSeparator()
        main_menu_bar.addMenu(self.help_menu)

    def connectSignalsAndSlots(self):

        self.fdsnDialog.fdsnWidget.sigTracesAppended.connect(self.waveformWidget.appendTraces)
        self.fdsnDialog.fdsnWidget.sigTracesReplaced.connect(self.waveformWidget.replaceTraces)

        # new connections
        self.sig_stream_changed.connect(self.waveformWidget.update_streams)
        self.sig_inventory_changed.connect(self.waveformWidget.update_inventory)

        self.beamformingWidget.waveformPlot.sigXRangeChanged.connect(self.waveformWidget.plotViewer.pl_widget.adjustSignalRegionRange)
        self.beamformingWidget.detectionWidget.signal_detections_changed.connect(self.locationWidget.update_detections)
        self.beamformingWidget.detectionWidget.signal_detections_cleared.connect(self.locationWidget.detections_cleared)

    def errorPopup(self, message):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setText(message)
        msgBox.setWindowTitle("Oops...")
        msgBox.exec_()

    def setStatus(self, s, ms=0):
        self.statusBar().showMessage(s, ms)

    def filemenu_NewProject(self):
        self._project = IPProject.IPProject()
        if self._project.makeNewProject():
            self.setWindowTitle(self.progname + ' - ' + self._project.get_projectName())
            self.setStatus(self._project.get_projectName() + ' successfully created', 5000)
            self.settings.setValue("last_baseProject_directory", str(self._project.get_basePath()))
            self.settings.setValue('last_project_directory', str(self._project.get_projectPath))

    def filemenu_LoadProject(self):
        newProject = IPProject.IPProject()
        if newProject.loadProject():
            self._project = newProject
            self.setWindowTitle(self.progname + ' - ' + self._project.get_projectName())
            self.setStatus(self._project.get_projectName() + ' successfully loaded', 5000)
            self.settings.setValue("last_baseProject_directory", str(self._project.get_basePath()))
            self.settings.setValue('last_project_directory', str(self._project.get_projectPath))

    def filemenu_CloseProject(self):
        self.setWindowTitle(self.progname)
        # self.consoleBox.append('Closed Project: ' + self._project.get_projectName())
        self.filemenu_ClearWaveforms()
        if self._project is not None:
            self.setStatus(self._project.get_projectName() + ' successfully closed', 5000)
            self._project = None

    def getProject(self):
        return self._project

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
        if self._project is None:
            # force a new filename...
            previous_directory = self.settings.value("last_open_directory", default_data_dir)
        else:
            # There is an open project, so make the default save location
            # correspond to what the project wants
            previous_directory = str(self._project.get_dataPath())

        ifiles = QFileDialog.getOpenFileNames(self, 'Open File...', previous_directory)

        if len(ifiles[0]) > 0:

            for ifile in ifiles[0]:
                ipath = os.path.dirname(ifile)
                if self._project is None:
                    self.settings.setValue("last_open_directory", ipath)
                else:
                    self._project.set_dataPath(ipath)
                try:
                    if self.waveformWidget._sts is not None:
                        self.waveformWidget._sts += obsRead(ifile)
                    else:
                        self.waveformWidget._sts = obsRead(ifile)
                except Exception:
                    self.setStatus("File Read Error", 5000)
                    continue

                for trace in self.waveformWidget._sts:
                    trace.data = trace.data - np.mean(trace.data)

                self.waveformWidget._sts.merge(fill_value=0)
        else:
            # No files were chosen to open
            return

        # it's possible, if the open failed, that self.waveformWidget._sts is still None, so if it is, bail out
        # if not populate the trace stats viewer and plot the traces
        if self.waveformWidget._sts is not None:
            self.beamformingWidget.setStreams(self.waveformWidget._sts)

            new_inventory = None
            for trace in self.waveformWidget._sts:
                if trace.stats['_format'] == 'SAC':
                    if new_inventory is None:
                        new_inventory = self.sac_trace_to_inventory(trace)
                    else:
                        new_inventory += self.sac_trace_to_inventory(trace)
                elif trace.stats['_format'] == 'MSEED':
                    # miniseed files have no metadata, so we need to deal with
                    # the inventory seperately

                    # First, there's a chance that the inventory data has been
                    # loaded already, so lets check the current inventory, and if it
                    # has been, leave it alone.  We do need to remove inventory that does
                    pass

            print(new_inventory)
            self.sig_stream_changed.emit(self.waveformWidget._sts)
            if new_inventory is not None:
                self.sig_inventory_changed.emit(new_inventory)

            self.mainTabs.setCurrentIndex(0)

            self.setStatus("Ready", 5000)
        else:
            return

    def filemenu_import(self):
        if self.fdsnDialog.exec_():
            self.mainTabs.setCurrentIndex(0)

    def filemenu_saveAllWaveforms(self):
        if self.waveformWidget._sts is None:
            self.errorPopup('Oops... No waveforms to save')
            return

        if self._project is None:
            # force a new filename...
            previousDirectory = self.settings.value("last_open_directory", QDir.homePath())
        else:
            # There is an open project, so make the default save location
            # correspond to what the project wants
            previousDirectory = str(self._project.get_dataPath())

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

    def view_menu_toggle_fullscreen(self):

        if self.is_fullscreen:
            self.showNormal()
            self.is_fullscreen = False
        else:
            self.showFullScreen()
            self.is_fullscreen = True

    def filemenu_ClearWaveforms(self):

        self.beamformingWidget.clearWaveformPlot()
        self.waveformWidget.clearWaveforms()

    def sac_trace_to_inventory(self, trace):
        # if sac files are opened, it's useful to extract inventory from their streams so that we can populate the stations tabs and the location widget
        new_inventory = None

        # The next bit is cribbed from the obspy webpage on building a stationxml site from scratch
        # https://docs.obspy.org/tutorial/code_snippets/stationxml_file_from_scratch.html
        #
        # We'll first create all the various objects. These strongly follow the
        # hierarchy of StationXML files.
        new_inventory = Inventory(
            # We'll add networks later.
            networks=[],
            # The source should be the id whoever create the file.
            source="InfraView")

        # Attempt to retrieve it from the sac header, if not found, set it to '---'
        _network = trace.stats['network']

        if _network == '':
            _network = '###'

        net = Network(
            # This is the network code according to the SEED standard.
            code=_network,
            # A list of stations. We'll add one later.
            stations=[],
            # Description isn't something that's in the SAC header, so lets set it to the network cod
            description=_network,
            # Start-and end dates are optional.

            # Start and end dates for the network are not stored in the sac header so lets set it to 1/1/1900
            start_date=UTCDateTime(1900, 1, 1))

        _station = trace.stats['station']
        if _station == '':
            _station = '###'

        if 'stla' in trace.stats['sac']:
            _lat = trace.stats['sac']['stla']
        else:
            self.errorPopup("SAC header doesn't contain latitude information")
            _lat = 10.0
        if 'stlo' in trace.stats['sac']:
            _lon = trace.stats['sac']['stlo']
        else:
            self.errorPopup("SAC header doesn't contain longitude information")
            _lon = 10.0
        if 'stel' in trace.stats['sac']:
            _ele = trace.stats['sac']['stel']
        else:
            _ele = -999.9

        sta = Station(
            # This is the station code according to the SEED standard.

            code=_station,
            latitude=_lat,
            longitude=_lon,
            elevation=_ele,
            # Creation_date is not saved in the sac header
            creation_date=UTCDateTime(1900, 1, 1),
            # Site name is not in the sac header, so set it to the site code
            site=Site(name=_station))

        _channel = trace.stats['channel']
        _location = trace.stats['location']
        # This is the channel code according to the SEED standard.
        cha = Channel(code=_channel,
                      # This is the location code according to the SEED standard.
                      location_code=_location,
                      # Note that these coordinates can differ from the station coordinates.
                      latitude=_lat,
                      longitude=_lon,
                      elevation=_ele,
                      depth=0.0)

        # Now tie it all together.
        # cha.response = response
        sta.channels.append(cha)
        net.stations.append(sta)
        new_inventory.networks.append(net)

        return new_inventory

    # ------------------------------------------------------------------------------
    # Settings methods

    def restoreSettings(self):
        # Restore settings

        self.settings.beginGroup('MainWindow')
        self.resize(self.settings.value("windowSize", QSize(1000, 900)))
        self.move(self.settings.value("windowPos", QPoint(200, 200)))
        self.settings.endGroup()

        self.beamformingWidget.restoreSettings()
        self.locationWidget.restoreSettings()
        self.waveformWidget.restoreSettings()

    def saveSettings(self):
        self.settings.beginGroup('MainWindow')
        self.settings.setValue("windowSize", self.size())
        self.settings.setValue("windowPos", self.pos())
        self.settings.endGroup()

        self.beamformingWidget.saveSettings()
        self.locationWidget.saveSettings()
        self.waveformWidget.saveSettings()

    # -------------------------------------------------------
    # Clean up

    def closeEvent(self, ce):
        self.saveSettings()

    # Obligatory about
    def about(self):
        QtWidgets.QMessageBox.about(self, "About", self.progname + "   " + self.progversion + "\n" + "Copyright 2018\nLos Alamos National Laboratory")
