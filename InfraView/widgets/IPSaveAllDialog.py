import sys

from pathlib import Path

from PyQt5 import QtCore, QtGui

from PyQt5.QtCore import Qt, QSettings, QDir

from PyQt5.QtWidgets import (QApplication, QDialog, QGridLayout, QFormLayout, QLabel, QLayout, QLineEdit,
                             QPushButton, QFileDialog, QWidget, QGroupBox, QVBoxLayout, QHBoxLayout,
                             QDialogButtonBox, QSizePolicy, QCheckBox, QMessageBox, QScrollArea)


class IPSaveAllDialog(QDialog):

    __streams = []
    __fileInfo = []
    __filtered_fileInfo = []
    __fileEdits = []
    __filtered_fileEdits = []

    def __init__(self, parent):
        super().__init__()
        self.__parent = parent
        self.__buildUI__()

    def exec_(self, streams='', directoryPath=None):

        self.__streams = streams

        self.settings = QSettings('LANL', 'InfraView')

        if directoryPath is None:
            self.__directoryName = self.settings.value("last_save_directory", QDir.homePath())
        else:
            self.__directoryName = str(directoryPath)

        self.lineEdit_Directory.setText(self.__directoryName)

        # manually call this slot to make sure things are populated correctly
        self.checkBoxClicked()

        return super().exec_()

    def __buildUI__(self):

        self.setWindowTitle('Save Data')
        self.setMinimumWidth(500)

        label_directory = QLabel(self.tr('Directory: '))
        self.lineEdit_Directory = QLineEdit()
        button_directory = QPushButton('Edit...')
        button_directory.clicked.connect(self.directoryDialog)

        pathWidget = QWidget()
        pathLayout = QGridLayout()
        pathLayout.addWidget(label_directory, 0, 0)
        pathLayout.addWidget(self.lineEdit_Directory, 0, 1)
        pathLayout.addWidget(button_directory, 0, 2)
        pathWidget.setLayout(pathLayout)

        label_saveFiltered = QLabel('Save Filtered Data: ')
        self.saveFiltered_check = QCheckBox()
        self.saveFiltered_check.clicked.connect(self.checkBoxClicked)
        label_saveOriginal = QLabel('Save Original Data: ')
        self.saveOriginal_check = QCheckBox()
        self.saveOriginal_check.clicked.connect(self.checkBoxClicked)

        hlayout = QHBoxLayout()
        hlayout.addWidget(label_saveOriginal)
        hlayout.addWidget(self.saveOriginal_check)
        hlayout.addWidget(label_saveFiltered)
        hlayout.addWidget(self.saveFiltered_check)

        self.fileWidget = QWidget()
        self.gridlayout1 = QGridLayout()
        self.gridlayout1.addWidget(self.fileWidget)

        self.fileGridLayout = QGridLayout()
        self.fileWidget.setLayout(self.fileGridLayout)

        self.fileGroupBox = QGroupBox('File Names')
        self.fileGroupBox.setLayout(self.gridlayout1)

        self.fileScrollArea = QScrollArea()
        # self.fileScrollArea.setWidget(self.fileWidget)
        self.fileScrollArea.setWidgetResizable(True)

        self.filteredGridLayout = QGridLayout()
        self.filteredGroupBox = QGroupBox('Filtered File Names')
        self.filteredGroupBox.setLayout(self.filteredGridLayout)
        self.filteredGroupBox.setVisible(False)

        filteredScrollArea = QScrollArea()
        # filteredScrollArea.setWidget(self.filteredGroupBox)
        filteredScrollArea.setWidgetResizable(True)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel, Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        mainLayout = QVBoxLayout()
        mainLayout.setSizeConstraint(QLayout.SetFixedSize)
        mainLayout.addWidget(pathWidget)
        mainLayout.addLayout(hlayout)
        mainLayout.addWidget(self.fileGroupBox)
        mainLayout.addWidget(self.filteredGroupBox)
        mainLayout.addStretch()
        mainLayout.addWidget(buttons)
        self.setLayout(mainLayout)

        # this needs to be called after the grid layouts are made...putting it at the end to make sure
        self.saveOriginal_check.setChecked(True)

    def generateOriginalDataFileInfo(self):
        for idx, trace in enumerate(self.__streams):
            stats = trace.stats
            fileFormat = stats['_format']

            basename = trace.id
            if basename[0] == '.':
                basename = basename[1:]
            filename = basename + '.' + fileFormat

            self.__fileInfo.append({'fname': filename, 'format': fileFormat, 'directory': self.__directoryName})

            self.__fileEdits.append(QLineEdit(filename))
            self.__fileEdits[-1].textChanged.connect(self.fileEditsChanged)
            self.fileGridLayout.addWidget(self.__fileEdits[idx], idx, 0)

    def generateFilteredDataFileInfo(self):
        for idx, trace in enumerate(self.__streams):
            stats = trace.stats
            fileFormat = stats['_format']

            basename = trace.id
            if basename[0] == '.':
                basename = basename[1:]
            filename = 'filtered.' + basename + '.' + fileFormat

            self.__filtered_fileInfo.append({'fname': filename, 'format': fileFormat, 'directory': self.__directoryName})

            self.__filtered_fileEdits.append(QLineEdit(filename))
            self.__filtered_fileEdits[idx].textChanged.connect(self.filtered_fileEditsChanged)
            self.filteredGridLayout.addWidget(self.__filtered_fileEdits[idx], idx, 0)

    def directoryDialog(self):
        self.__directoryName = QFileDialog.getExistingDirectory(self, "Choose a Directory", self.__directoryName, QtGui.QFileDialog.ShowDirsOnly)
        if self.__directoryName != '':
            self.settings.setValue("last_save_directory", self.__directoryName)
            self.lineEdit_Directory.setText(self.__directoryName)

    def getSaveDirectory(self):
        return self.lineEdit_Directory.text()

    def getFileInfo(self):
        return self.__fileInfo

    def getFilteredFileInfo(self):
        return self.__filtered_fileInfo

    def getFileChoiceData(self):
        # This is the for the checkboxes for whether to save the original data, filtered data, or both
        if self.saveOriginal_check.isChecked() and not self.saveFiltered_check.isChecked():
            # Save just the original data
            return 1
        elif self.saveFiltered_check.isChecked() and not self.saveOriginal_check.isChecked():
            # Save just the filtered data
            return 2
        elif self.saveOriginal_check.isChecked() and self.saveFiltered_check.isChecked():
            # Save both
            return 3
        else:
            # Peculiar case where neither is checked
            return 0

    def errorPopup(self, message):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setText(message)
        msgBox.setWindowTitle("Oops...")
        msgBox.exec_()

    @QtCore.pyqtSlot()
    def fileEditsChanged(self):
        # slot function called when a fileEdit box is edited
        for idx, newFileName in enumerate(self.__fileEdits):
            self.__fileInfo[idx]['fname'] = self.__fileEdits[idx].text()

    @QtCore.pyqtSlot()
    def filtered_fileEditsChanged(self):
        # slot function called when a fileEdit box is edited
        for idx, newFileName in enumerate(self.__filtered_fileEdits):
            self.__filtered_fileInfo[idx]['fname'] = self.__filtered_fileEdits[idx].text()

    @QtCore.pyqtSlot()
    def checkBoxClicked(self):

        if self.saveFiltered_check.isChecked():
            filterDisplaySettings = self.__parent.waveformWidget.filterSettingsWidget.get_filter_display_settings()

            if filterDisplaySettings['apply'] is False:
                # self.saveFiltered_check.blockSignals(True)
                self.saveFiltered_check.setChecked(False)
                # self.saveFiltered_check.blockSignals(False)

                self.errorPopup('Filter is not currently applied to any data')
                return

        # clear out previous file info
        self.__fileInfo.clear()

        # Clear out the layouts
        for i in reversed(range(self.fileGridLayout.count())):
            self.fileGridLayout.itemAt(i).widget().setParent(None)

        for i in reversed(range(self.filteredGridLayout.count())):
            self.filteredGridLayout.itemAt(i).widget().setParent(None)

        # Repopulate the fileInfo and make the new lineedits
        if self.saveOriginal_check.isChecked():
            self.fileGroupBox.setVisible(True)
            self.generateOriginalDataFileInfo()
        else:
            self.fileGroupBox.setVisible(False)

        if self.saveFiltered_check.isChecked():
            self.filteredGroupBox.setVisible(True)
            self.generateFilteredDataFileInfo()
        else:
            self.filteredGroupBox.setVisible(False)
