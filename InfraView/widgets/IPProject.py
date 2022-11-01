import sys

from pathlib import Path

from PyQt5.QtCore import Qt, QSettings
from PyQt5.QtWidgets import (QDialog, QGridLayout, QLabel, QLineEdit,
                             QPushButton, QFileDialog, QWidget, QVBoxLayout,
                             QDialogButtonBox, QSizePolicy)


class IPProject:

    basePath = None                 # base path projects are saved in
    projectPath = None              # path of actual project
    dataPath = None                 # where locally stored data (can) be saved
    beamformingResutsPath = None    # path to csv files holding fstat, ba, and tracev data
    detectionsPath = None           # where picks will be saved
    customFilterPath = None         # where custom filters will be saved
    homePath = None                 # user's home directory
    stationsPath = None             # where station xml files will be saved
    eventPath = None                # where event json files will be saved

    projectName = None
    projectFileName = None

    def __init__(self):
        self.globalSettings = QSettings('LANL', 'InfraView')

    def makeNewProject(self):
        newDialog = IPNewProjectDialog(self)
        if newDialog.exec_():
            self.basePath, self.projectName = newDialog.getBasePathAndProjectName()
            self.projectPath = Path(str(self.basePath) + '/' + self.projectName)
            self.dataPath = Path(str(self.projectPath) + '/data')
            self.detectionsPath = Path(str(self.projectPath) + '/detections')
            self.stationsPath = Path(str(self.projectPath) + '/stations')
            self.customFilterPath = Path(str(self.projectPath) + '/customFilters')
            self.beamformingResutsPath = Path(str(self.projectPath) + '/beamformingResults')
            self.eventPath = Path(str(self.projectPath) + '/events')

            # Create the project directories
            self.projectPath.mkdir(parents=True, exist_ok=True)
            self.dataPath.mkdir(parents=True, exist_ok=True)
            self.detectionsPath.mkdir(parents=True, exist_ok=True)
            self.stationsPath.mkdir(parents=True, exist_ok=True)
            self.customFilterPath.mkdir(parents=True, exist_ok=True)
            self.beamformingResutsPath.mkdir(parents=True, exist_ok=True)
            self.eventPath.mkdir(parents=True, exist_ok=True)

            # Create a settings object/file for the new project and populate it with the directories
            self.projectFileName = self.projectName + '.ipprj'
            self.projectSettings = QSettings(str(self.projectPath) + '/' + self.projectFileName, QSettings.IniFormat)
            self.projectSettings.beginGroup('Main')
            self.projectSettings.setValue('projectName', str(self.projectName))
            self.projectSettings.endGroup()

            self.projectSettings.beginGroup('PathNames')
            self.projectSettings.setValue('basePathName', str(self.basePath))
            self.projectSettings.setValue('projectPathName', str(self.projectPath))
            self.projectSettings.setValue('dataPathName', str(self.dataPath))
            self.projectSettings.setValue('detectionsPathName', str(self.detectionsPath))
            self.projectSettings.setValue('stationsPathName', str(self.stationsPath))
            self.projectSettings.setValue('customFilterPathName', str(self.customFilterPath))
            self.projectSettings.setValue('beamformingResultsPath', str(self.beamformingResutsPath))
            self.projectSettings.setValue('eventPath', str(self.eventPath))
            self.projectSettings.endGroup()

            # From here down are settings that are the default for infrapy's CLI interface.  These should 
            # be kept up-to-date with what is in that file.  Not all of these are used by infraview
            self.projectSettings.beginGroup('FK')
            self.projectSettings.setValue('freq_min', 0.5)
            self.projectSettings.setValue('freq_max', 5.0)
            self.projectSettings.setValue('back_az_min', -180.0)
            self.projectSettings.setValue('back_az_max', 180.0)
            self.projectSettings.setValue('back_az_step', 2.0)
            self.projectSettings.setValue('trace_vel_min', 300.0)
            self.projectSettings.setValue('trace_vel_max', 600.0)
            self.projectSettings.setValue('trace_vel_step', 2.5)
            self.projectSettings.setValue('method', 'bartlett')
            self.projectSettings.setValue('signal_start', 'None')
            self.projectSettings.setValue('signal_end', 'None')
            self.projectSettings.setValue('noise_start', 'None')
            self.projectSettings.setValue('noise_end', 'None')
            self.projectSettings.setValue('window_len', 10)
            self.projectSettings.setValue('sub_window_len', 'None')
            self.projectSettings.setValue('window_step', 5)
            self.projectSettings.setValue('cpu_cnt', 'None')
            self.projectSettings.endGroup()

            self.projectSettings.beginGroup('FD')
            self.projectSettings.setValue('window_len', 3600)
            self.projectSettings.setValue('p_value', 0.01)
            self.projectSettings.setValue('min_duration', 10)
            self.projectSettings.setValue('back_az_width', 15.0)
            self.projectSettings.setValue('fixed_thresh', 'None')
            self.projectSettings.setValue('thresh_ceil', 'None')
            self.projectSettings.setValue('return_thresh', False)
            self.projectSettings.setValue('merge_dets', False)
            self.projectSettings.endGroup()

            self.projectSettings.beginGroup('ASSOC')
            self.projectSettings.setValue('back_az_width', 10.0)
            self.projectSettings.setValue('range_max', 2000.0)
            self.projectSettings.setValue('resolution', 180)
            self.projectSettings.setValue('distance_matrix_max', 8.0)
            self.projectSettings.setValue('cluster_linkage', 'weighted')
            self.projectSettings.setValue('cluster_threshold', 5.0)
            self.projectSettings.setValue('trimming_threshold', 3.8)
            self.projectSettings.setValue('event_population_min', 3)
            self.projectSettings.setValue('event_station_min', 2)
            self.projectSettings.setValue('multithread', False)
            self.projectSettings.setValue('cpu_cnt', 'None')
            self.projectSettings.endGroup()

            self.projectSettings.beginGroup('LOC')
            self.projectSettings.setValue('back_az_width', 10.0)
            self.projectSettings.setValue('range_max', 2000.0)
            self.projectSettings.setValue('resolution', 180)
            self.projectSettings.setValue('src_est', 'None')
            self.projectSettings.setValue('pgm_model', 'None')
            self.projectSettings.endGroup()

            self.projectSettings.beginGroup('YIELD')
            self.projectSettings.setValue('source_loc', [30.0, -105])
            self.projectSettings.setValue('freq_min', 0.25)
            self.projectSettings.setValue('freq_max', 1.0)
            self.projectSettings.setValue('yld_min', 1)
            self.projectSettings.setValue('yld_max', 1000)
            self.projectSettings.setValue('ref_rng', 1.0)
            self.projectSettings.setValue('resolution', 200)
            self.projectSettings.setValue('noise_option', 'post')
            self.projectSettings.setValue('window_buffer', 0.2)
            self.projectSettings.setValue('amp_press', 101.325)
            self.projectSettings.setValue('amp_temp', 288.15)
            self.projectSettings.setValue('grnd_burst', True)
            self.projectSettings.setValue('exp_type', 'chemical')
            self.projectSettings.endGroup()

            self.projectSettings.beginGroup('VISUALIZATION')
            self.projectSettings.setValue('offline_maps_dir', '')
            self.projectSettings.endGroup()

            return True

        else:

            return False

    def loadProject(self):
        mydirectory = self.globalSettings.value('last_baseProject_directory', self.homePath)
        ipprjPathname, _ = QFileDialog.getOpenFileName(caption='Open InfraView Project', directory=mydirectory, filter='InfraView Project Files (*.ipprj)')
        if ipprjPathname:
            self.projectSettings = QSettings(ipprjPathname, QSettings.IniFormat)
            self.projectSettings.beginGroup('Main')
            self.projectName = self.projectSettings.value('projectName')
            self.projectFileName = self.projectName + '.ipprj'
            self.projectSettings.endGroup()
            self.projectSettings.beginGroup('PathNames')
            self.basePath = Path(self.projectSettings.value('basePathName'))
            self.projectPath = Path(self.projectSettings.value('projectPathName'))
            self.dataPath = Path(self.projectSettings.value('dataPathName'))
            self.detectionsPath = Path(self.projectSettings.value('detectionsPathName'))
            self.stationsPath = Path(self.projectSettings.value('stationsPathName'))
            self.customFilterPath = Path(self.projectSettings.value('customFilterPathName'))
            # when opening old projects, newer settings might not be present
            if self.projectSettings.value('beamformingResultsPath') is None:
                self.beamformingResutsPath = Path(str(self.projectPath) + '/beamformingResults')
                # the directory might not exist at this point, so create it
                Path(self.beamformingResutsPath).mkdir(parents=True, exist_ok=True)
            else:
                self.beamformingResutsPath = Path(self.projectSettings.value('beamformingResultsPath'))
            if self.projectSettings.value('eventPath') is None:
                self.eventPath = Path(str(self.projectPath) + '/events')
                # the directory might not exist at this point, so create it
                Path(self.eventPath).mkdir(parents=True, exist_ok=True)
            else:
                self.eventPath = Path(self.projectSettings.value('eventPath'))

            self.projectSettings.endGroup()
            return True
        else:
            return False

    def get_basePath(self):
        return self.basePath

    def get_projectPath(self):
        return self.projectPath

    def get_dataPath(self):
        return self.dataPath

    def set_dataPath(self, path):
        self.dataPath = path

    def get_detectionsPath(self):
        return self.detectionsPath

    def get_stationsPath(self):
        return self.stationsPath

    def get_customFilterPath(self):
        return self.customFilterPath

    def get_projectName(self):
        return self.projectName

    def get_projectFileName(self):
        return self.projectFileName

    def get_beamformResultsPath(self):
        return self.beamformingResutsPath

    def get_eventPath(self):
        return self.eventPath

    def clear(self):
        self.basePath = None                # base path projects are saved in
        self.projectPath = None             # path of actual project
        self.dataPath = None                # where locally stored data (can) be saved
        self.detectionsPath = None          # where picks will be saved
        self.stationsPath = None            # where exported picks will be saved
        self.customFilterPath = None        # where custom filters will be saved
        self.homePath = None                # user's home directory
        self.projectName = None
        self.projectFileName = None
        self.beamformingResutsPath = None   # beamforming results directory
        self.eventPath = None


class IPNewProjectDialog(QDialog):

    basePath = None
    projectName = None

    def __init__(self, parent):
        super().__init__()

        homePath = Path.home()
        self.basePath = Path(homePath, 'IPProjects')

        self.buildUI()

    def buildUI(self):
        self.setWindowTitle('InfraView: Create a New Project')
        label_projectName = QLabel(self.tr('Project Name: '))
        self.lineEdit_projectName = QLineEdit()
        self.lineEdit_projectName.textChanged.connect(self.updateProjectPath)

        label_basePath = QLabel(self.tr('Base Directory: '))
        self.lineEdit_basePath = QLineEdit(str(self.basePath))
        self.lineEdit_basePath.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.lineEdit_basePath.setMinimumWidth(400)
        self.lineEdit_basePath.textChanged.connect(self.updateProjectPath)
        button_basePathEdit = QPushButton('Edit...')
        button_basePathEdit.clicked.connect(self.directoryDialog)

        self.label_projectDirectory = QLabel('Project Directory: ')
        self.label_projectDirectory_value = QLabel(str(self.basePath) + '/' + self.lineEdit_projectName.text())

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        gridWidget = QWidget()
        gridlayout = QGridLayout()
        gridlayout.addWidget(label_projectName, 0, 0)
        gridlayout.addWidget(self.lineEdit_projectName, 0, 1)
        gridlayout.addWidget(label_basePath, 1, 0)
        gridlayout.addWidget(self.lineEdit_basePath, 1, 1)
        gridlayout.addWidget(button_basePathEdit, 1, 2)
        gridlayout.addWidget(self.label_projectDirectory, 2, 0)
        gridlayout.addWidget(self.label_projectDirectory_value, 2, 1)
        gridWidget.setLayout(gridlayout)

        mainLayout = QVBoxLayout()
        mainLayout.addWidget(gridWidget)
        mainLayout.addWidget(buttons)

        self.setLayout(mainLayout)

    def updateProjectPath(self):
        self.basePath = self.lineEdit_basePath.text()
        self.projectName = self.lineEdit_projectName.text()
        self.label_projectDirectory_value.setText(self.lineEdit_basePath.text() + '/' + self.lineEdit_projectName.text())

    def directoryDialog(self):
        newBasePathName = QFileDialog.getExistingDirectory(
            self, "Choose a Directory", str(self.basePath), QFileDialog.ShowDirsOnly)
        if newBasePathName != '':
            # self.settings.setValue("last_projectbase_directory", newBasePathName)
            self.lineEdit_basePath.setText(newBasePathName)
            self.basePath = Path(newBasePathName)

    def getBasePathAndProjectName(self):
        return self.basePath, self.projectName
