import sys

from pathlib import Path

from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import Qt, QObject, QSettings, QDir
from PyQt5.QtWidgets import (QApplication, QDialog, QGridLayout, QLabel, QLineEdit,
                             QPushButton, QFileDialog, QWidget, QGroupBox, QVBoxLayout,
                             QDialogButtonBox, QSizePolicy)


class IPProject:

    __basePath = None           # base path projects are saved in
    __projectPath = None        # path of actual project
    __dataPath = None           # where locally stored data (can) be saved
    __detectionsPath = None           # where picks will be saved
    __pickExportsPath = None    # where exported picks will be saved
    __customFilterPath = None   # where custom filters will be saved
    __homePath = None           # user's home directory
    __stationsPath = None       # where station xml files will be saved

    __projectName = None
    __projectFileName = None

    def __init__(self):
        self.__globalSettings = QSettings('LANL', 'InfraView')
        print(self.__globalSettings)

    def makeNewProject(self):
        newDialog = IPNewProjectDialog(self)
        if newDialog.exec_():
            self.__basePath, self.__projectName = newDialog.getBasePathAndProjectName()
            self.__projectPath = Path(str(self.__basePath) + '/' + self.__projectName)
            self.__dataPath = Path(str(self.__projectPath) + '/data')
            self.__detectionsPath = Path(str(self.__projectPath) + '/detections')
            self.__stationsPath = Path(str(self.__projectPath) + '/stations')
            self.__customFilterPath = Path(str(self.__projectPath) + '/customFilters')

            # Create the project directories
            self.__projectPath.mkdir(parents=True, exist_ok=True)
            self.__dataPath.mkdir(parents=True, exist_ok=True)
            self.__detectionsPath.mkdir(parents=True, exist_ok=True)
            self.__stationsPath.mkdir(parents=True, exist_ok=True)
            self.__customFilterPath.mkdir(parents=True, exist_ok=True)

            # Create a settings object/file for the new project and populate it with the directories
            self.__projectFileName = self.__projectName + '.ipprj'
            self.projectSettings = QSettings(str(self.__projectPath) + '/' + self.__projectFileName, QSettings.IniFormat)
            self.projectSettings.beginGroup('Main')
            self.projectSettings.setValue('projectName', str(self.__projectName))
            self.projectSettings.endGroup()
            self.projectSettings.beginGroup('PathNames')
            self.projectSettings.setValue('basePathName', str(self.__basePath))
            self.projectSettings.setValue('projectPathName', str(self.__projectPath))
            self.projectSettings.setValue('dataPathName', str(self.__dataPath))
            self.projectSettings.setValue('detectionsPathName', str(self.__detectionsPath))
            self.projectSettings.setValue('stationsPathName', str(self.__stationsPath))
            self.projectSettings.setValue('customFilterPathName', str(self.__customFilterPath))

            self.projectSettings.endGroup()

            return True

        else:

            return False

    def loadProject(self):
        mydirectory = self.__globalSettings.value('last_baseProject_directory', self.__homePath)
        ipprjPathname, _ = QFileDialog.getOpenFileName(
            caption='Open InfraView Project', directory=mydirectory, filter='InfraView Project Files (*.ipprj)')
        if ipprjPathname:
            self.projectSettings = QSettings(ipprjPathname, QSettings.IniFormat)
            self.projectSettings.beginGroup('Main')
            self.__projectName = self.projectSettings.value('projectName')
            self.__projectFileName = self.__projectName + '.ipprj'
            self.projectSettings.endGroup()
            self.projectSettings.beginGroup('PathNames')
            self.__basePath = Path(self.projectSettings.value('basePathName'))
            self.__projectPath = Path(self.projectSettings.value('projectPathName'))
            self.__dataPath = Path(self.projectSettings.value('dataPathName'))
            self.__detectionsPath = Path(self.projectSettings.value('detectionsPathName'))
            self.__stationsPath = Path(self.projectSettings.value('stationsPathName'))
            self.__customFilterPath = Path(self.projectSettings.value('customFilterPathName'))
            self.projectSettings.endGroup()
            return True
        else:
            return False

    def get_basePath(self):
        return self.__basePath

    def get_projectPath(self):
        return self.__projectPath

    def get_dataPath(self):
        return self.__dataPath

    def set_dataPath(self, path):
        self.__dataPath = path

    def get_detectionsPath(self):
        return self.__detectionsPath

    def get_stationsPath(self):
        return self.__stationsPath

    def get_customFilterPath(self):
        return self.__customFilterPath

    def get_projectName(self):
        return self.__projectName

    def get_projectFileName(self):
        return self.__projectFileName

    def clear():
        __basePath = None           # base path projects are saved in
        __projectPath = None        # path of actual project
        __dataPath = None           # where locally stored data (can) be saved
        __detectionsPath = None     # where picks will be saved
        __stationsPath = None       # where exported picks will be saved
        __customFilterPath = None   # where custom filters will be saved
        __homePath = None           # user's home directory

        __projectName = None
        __projectFileName = None


class IPNewProjectDialog(QDialog):

    __basePath = None
    __projectName = None

    def __init__(self, parent):
        super().__init__()

        __homePath = Path.home()
        self.__basePath = Path(__homePath, 'IPProjects')

        self.buildUI()

    def buildUI(self):
        self.setWindowTitle('Create a New Project')
        label_projectName = QLabel(self.tr('Project Name: '))
        self.lineEdit_projectName = QLineEdit()
        self.lineEdit_projectName.textChanged.connect(self.updateProjectPath)

        label_basePath = QLabel(self.tr('Base Directory: '))
        self.lineEdit_basePath = QLineEdit(str(self.__basePath))
        self.lineEdit_basePath.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.lineEdit_basePath.setMinimumWidth(400)
        self.lineEdit_basePath.textChanged.connect(self.updateProjectPath)
        button_basePathEdit = QPushButton('Edit...')
        button_basePathEdit.clicked.connect(self.directoryDialog)

        self.label_projectDirectory = QLabel('Project Directory: ')
        self.label_projectDirectory_value = QLabel(str(self.__basePath) + '/' + self.lineEdit_projectName.text())

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
        self.__basePath = self.lineEdit_basePath.text()
        self.__projectName = self.lineEdit_projectName.text()
        self.label_projectDirectory_value.setText(self.lineEdit_basePath.text() + '/' + self.lineEdit_projectName.text())

    def directoryDialog(self):
        newBasePathName = QFileDialog.getExistingDirectory(
            self, "Choose a Directory", str(self.__basePath), QFileDialog.ShowDirsOnly)
        if newBasePathName != '':
            # self.settings.setValue("last_projectbase_directory", newBasePathName)
            self.lineEdit_basePath.setText(newBasePathName)
            self.__basePath = Path(newBasePathName)

    def getBasePathAndProjectName(self):
        return self.__basePath, self.__projectName
