import os, json

from PyQt5.QtWidgets import (QDateEdit, QPushButton, QLabel, QLineEdit, QFileDialog, QHBoxLayout,
                             QFormLayout, QDoubleSpinBox, QSpinBox, QVBoxLayout, QTimeEdit, QWidget,
                             QCheckBox)
from PyQt5.QtCore import QDir, QTime, QDate, QSettings, pyqtSignal
from PyQt5.QtGui import QIcon

from InfraView.widgets import IPEventBrowser

class IPEventWidget(QWidget):

    sigEventWidgetLoaded = pyqtSignal()
    sigEventWidgetChanged = pyqtSignal(dict)
    sigEventCleared = pyqtSignal()
    
    savefile = None
    parent = None

    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self.parent = parent 

        self.settings = QSettings('LANL', 'InfraView')

        self.buildUI()

    def buildUI(self):

        self.buildIcons()

        formWidget = QWidget()
        formLayout = QFormLayout()

        self.displayEvent_cb = QCheckBox(self.tr('Show on waveform plots'))
        self.displayEvent_cb.setChecked(False)
        self.displayEvent_cb.stateChanged.connect(self.eventChanged)

        label_event_name = QLabel(self.tr('Event ID: '))
        self.event_name_edit = QLineEdit()
        self.event_name_edit.textChanged.connect(self.eventChanged)

        label_latitude = QLabel(self.tr('Latitude'))
        self.event_lat_edit = QDoubleSpinBox()
        self.event_lat_edit.setRange(-90.1,90.0)
        self.event_lat_edit.setDecimals(8)
        self.event_lat_edit.setSpecialValueText('--')
        self.event_lat_edit.setSingleStep(0.1)
        self.event_lat_edit.setValue(self.event_lat_edit.minimum())
        self.event_lat_edit.valueChanged.connect(self.eventChanged)

        label_longitude = QLabel(self.tr('Longitude'))
        self.event_lon_edit = QDoubleSpinBox()
        self.event_lon_edit.setRange(-180.1, 180.0)
        self.event_lon_edit.setDecimals(8)
        self.event_lon_edit.setSpecialValueText('--')
        self.event_lon_edit.setSingleStep(0.1)
        self.event_lon_edit.setValue(self.event_lon_edit.minimum())
        self.event_lon_edit.valueChanged.connect(self.eventChanged)

        label_event_date = QLabel(self.tr('Date (UTC):'))
        self.event_date_edit = QDateEdit()
        self.event_date_edit.setDisplayFormat("yyyy-MM-dd")
        self.event_date_edit.setDate(self.event_date_edit.minimumDate())
        self.event_date_edit.setSpecialValueText("yyyy-MM-dd")
        self.event_date_edit.dateChanged.connect(self.eventChanged)

        label_event_time = QLabel(self.tr('Time (UTC):'))
        self.event_time_edit = QTimeEdit()
        self.event_time_edit.setDisplayFormat('HH:mm:ss.zzz')
        self.event_time_edit.setTime(self.event_time_edit.minimumTime())
        self.event_time_edit.setSpecialValueText('HH:mm:ss.zzz')
        self.event_time_edit.timeChanged.connect(self.eventChanged)

        label_event_elev = QLabel(self.tr('Elevation'))
        self.event_elev_edit = QSpinBox()
        self.event_elev_edit.setRange(-1000000, 1000000)
        self.event_elev_edit.setSpecialValueText('--')
        self.event_elev_edit.setValue(self.event_elev_edit.minimum())
        self.event_elev_edit.setSuffix(' m')
        self.event_elev_edit.valueChanged.connect(self.eventChanged)

        label_event_depth = QLabel(self.tr('Depth'))
        self.event_depth_edit = QSpinBox()
        self.event_depth_edit.setRange(-1000000, 1000000)
        self.event_depth_edit.setSpecialValueText('--')
        self.event_depth_edit.setValue(self.event_depth_edit.minimum())
        self.event_depth_edit.setSuffix(' m')
        self.event_depth_edit.valueChanged.connect(self.eventChanged)

        label_event_mag = QLabel(self.tr('Magnitude'))
        self.event_mag_edit = QDoubleSpinBox()
        self.event_mag_edit.setRange(0.0, 1000.0)
        self.event_mag_edit.setSpecialValueText('--')
        self.event_mag_edit.setSingleStep(0.1)
        self.event_mag_edit.setValue(self.event_mag_edit.minimum())
        self.event_mag_edit.valueChanged.connect(self.eventChanged)

        self.load_button = QPushButton(self.tr(' Load Event...'))
        self.load_button.setMaximumWidth(100)
        self.load_button.clicked.connect(self.loadEvent)
        self.load_button.setIcon(self.openIcon)

        self.update_button = QPushButton(self.tr('Update'))
        self.update_button.setMaximumWidth(100)
        self.update_button.clicked.connect(self.eventChanged)

        self.save_button = QPushButton(self.tr(' Save Event...'))
        self.save_button.setMaximumWidth(100)
        self.save_button.clicked.connect(self.saveEventAs)
        self.save_button.setIcon(self.saveAsIcon)

        self.clear_button = QPushButton(self.tr(' Clear'))
        self.save_button.setMaximumWidth(100)
        self.clear_button.clicked.connect(self.clear)
        self.clear_button.setIcon(self.clearIcon)
        
        self.browse_button = QPushButton(self.tr(' IRIS Event Browser...'))
        self.browse_button.setMaximumWidth(200)
        self.browse_button.clicked.connect(self.browse)

        show_layout = QHBoxLayout()
        show_layout.addStretch()
        show_layout.addWidget(self.displayEvent_cb)
        show_layout.addStretch()

        formLayout.addRow(label_event_name, self.event_name_edit)
        formLayout.addRow(label_longitude, self.event_lon_edit)
        formLayout.addRow(label_latitude, self.event_lat_edit)
        formLayout.addRow(label_event_date, self.event_date_edit)
        formLayout.addRow(label_event_time, self.event_time_edit)
        formLayout.addRow(label_event_elev, self.event_elev_edit)
        formLayout.addRow(label_event_depth, self.event_depth_edit)
        formLayout.addRow(label_event_mag, self.event_mag_edit)

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(self.load_button)
        buttonLayout.addWidget(self.save_button)
        buttonLayout.addWidget(self.clear_button)
        buttonLayout.addWidget(self.update_button)
        buttonLayout.addStretch()

        browseLayout = QHBoxLayout()
        browseLayout.addStretch()
        browseLayout.addWidget(self.browse_button)
        browseLayout.addStretch()

        formWidget.setLayout(formLayout)

        verticalLayout = QVBoxLayout()
        # verticalLayout.addLayout(show_layout)
        verticalLayout.addWidget(formWidget)
        verticalLayout.addLayout(buttonLayout)
        verticalLayout.addLayout(browseLayout)
        verticalLayout.addStretch()

        self.setLayout(verticalLayout)

        self.eventBrowser = IPEventBrowser.IPEventDialog()

        self.show()

    def buildIcons(self):
        self.clearIcon = QIcon.fromTheme("edit-clear")
        self.openIcon = QIcon.fromTheme("document-open")
        self.saveIcon = QIcon.fromTheme("document-save")
        self.saveAsIcon = QIcon.fromTheme("document-save-as")

    def getID(self):
        return self.event_name_edit.text()

    def getLat(self):
        return self.event_lat_edit.value()

    def getLon(self):
        return self.event_lon_edit.value()

    def getDepth(self):
        return self.event_depth_edit.value()

    def getElevation(self):
        return self.event_elev_edit.value()

    def getMag(self):
        return self.event_mag_edit.value()

    def getUTCDateTimeString(self):
        date = self.event_date_edit.date().toPyDate()
        time = self.event_time_edit.time().toPyTime()
        utcString = str(date) + 'T' + str(time)
        return utcString
        
    def setUTCDateTime(self, newUTCTime_iso):
        datetime = newUTCTime_iso.isoformat()

        date = datetime[0:10]
        self.event_date_edit.setDate(QDate.fromString(date, 'yyyy-MM-dd'))

        time = datetime[11:]
        self.event_time_edit.setTime(QTime.fromString(time[0:12], "hh:mm:ss.zzz"))

    def hasValidEvent(self):
        #returns true if the minimum info needed for an event is set... Date, Lat, and Lon
        if (self.event_lat_edit.value() == self.event_lat_edit.minimum() or self.event_lon_edit.value() == self.event_lon_edit.minimum() or self.event_date_edit.date() == self.event_date_edit.minimumDate()):
            return False
        else:
            return True

    def eventChanged(self):
        # whenever any widget changes, this signal is emitted with the new event dictionary
        self.sigEventWidgetChanged.emit(self.Dict())

    def Dict(self):
        # Returns a dictionary representation of the data in the widget
        if self.event_lat_edit.value() == self.event_lat_edit.minimum():
            lat = None
        else:
            lat = self.event_lat_edit.value()

        if self.event_lon_edit.value() == self.event_lon_edit.minimum():
            lon = None
        else:
            lon = self.event_lon_edit.value()

        if self.event_elev_edit.value() == self.event_elev_edit.minimum():
            elev = None
        else:
            elev = self.event_elev_edit.value()

        if self.event_depth_edit.value() == self.event_depth_edit.minimum():
            depth = None
        else:
            depth = self.event_depth_edit.value()

        if self.event_mag_edit.value() == self.event_mag_edit.minimum():
            mag = None
        else:
            mag = self.event_mag_edit.value()

        eventdict = {'Name':self.event_name_edit.text(), 
                    'UTC Date':str(self.event_date_edit.date().toPyDate()), 
                    'UTC Time':str(self.event_time_edit.time().toPyTime()), 
                    'Longitude':lon, 
                    'Latitude':lat, 
                    'Elevation':elev, 
                    'Depth':depth,
                    'Magnitude':mag }
        return eventdict

    def saveEventAs(self):
        # pop up a save file dialog, default to project directory if a project is open, otherwise use the last used directory
        if self.parent.getProject() is None:
            # force a new filename...
            savePath=self.settings.value("last_eventfile_directory", QDir.homePath())
        else:
            # There is an open project, so make the default save location correspond to what the project wants
            savePath = str(self.parent.getProject().get_eventPath())

        self.savefile = QFileDialog.getSaveFileName(self, 'Save File', savePath)
        if self.savefile[0]:
            with open(self.savefile[0], 'w') as of:
                json.dump(self.Dict(), of, indent=4)

                if self.parent.getProject() is None:
                    # if there is no open project, update the global settings 
                    self.settings.setValue("last_eventfile_directory", os.path.dirname(self.savefile[0]))

    def loadEvent(self):
        if self.parent.getProject() is None:
            loadPath=self.settings.value("last_eventfile_directory", QDir.homePath())
        else:
            # There is an open project, so make the default save location correspond to what the project wants
            loadPath = str(self.parent.getProject().get_eventPath())

        self.__openfile = QFileDialog.getOpenFileName(self, 'Open File', loadPath)
        if self.__openfile[0]:
            with open(self.__openfile[0], 'r') as infile:
                new_event = json.load(infile)

            if new_event['Name'] is not None:
                self.event_name_edit.setText(new_event['Name'])
            else:
                self.event_name_edit.setText('')

            if new_event['Longitude'] is not None:
                self.event_lon_edit.setValue(float(new_event['Longitude']))
            else:
                self.event_lon_edit.setValue(self.event_lon_edit.minimum())

            if new_event['Latitude'] is not None:
                self.event_lat_edit.setValue(float(new_event['Latitude']))
            else:
                self.event_lat_edit.setValue(self.event_lat_edit.minimum())

            if new_event['Elevation'] is not None:
                self.event_elev_edit.setValue(int(new_event['Elevation']))
            else:
                self.event_elev_edit.setValue(self.event_elev_edit.minimum())

            if new_event['Depth'] is not None:
                self.event_depth_edit.setValue(int(new_event['Depth']))
            else:
                self.event_depth_edit.setValue(self.event_depth_edit.minimum())

            if new_event['UTC Time'] is not None:
                if len(new_event['UTC Time']) == 8:
                    self.event_time_edit.setTime(QTime.fromString(new_event['UTC Time'][0:8], "hh:mm:ss"))
                elif len(new_event['UTC Time']) >=12:
                    self.event_time_edit.setTime(QTime.fromString(new_event['UTC Time'][0:12], "hh:mm:ss.zzz"))
            else:
                self.event_time_edit.setTime(self.event_time_edit.minimumTime())

            if new_event['UTC Date'] is not None:
                self.event_date_edit.setDate(QDate.fromString(new_event['UTC Date'], 'yyyy-MM-dd'))
            else:
                self.event_date_edit.setDate(self.event_date_edit.minimumDate())

            if new_event['Magnitude'] is not None:
                self.event_mag_edit.setValue(float(new_event['Magnitude']))

            if self.parent.getProject() is None:
                # if there is no open project, update the global settings 
                self.settings.setValue("last_eventfile_directory", os.path.dirname(self.__openfile[0]))

            self.sigEventWidgetLoaded.emit()

    def browse(self):
        self.eventDialog = IPEventBrowser.IPEventDialog()

        if self.eventDialog.exec_():
            event = self.eventDialog.getEvent()

            self.event_name_edit.setText(event['Name'])
            self.event_lon_edit.setValue(event['Longitude'])
            self.event_lat_edit.setValue(event['Latitude'])
            self.event_depth_edit.setValue(event['Depth'])
            self.event_elev_edit.setValue(event['Elevation'])
            self.event_mag_edit.setValue(event['Magnitude'])
            self.event_date_edit.setDate(QDate.fromString(event['UTC Date'], 'yyyy-MM-dd'))
            self.event_time_edit.setTime(QTime.fromString(event['UTC Time'], 'hh:mm:ss.zzz'))
            self.sigEventWidgetLoaded.emit()

    def clear(self):
        self.event_name_edit.setText('')
        self.event_lon_edit.setValue(self.event_lon_edit.minimum())
        self.event_lat_edit.setValue(self.event_lat_edit.minimum())
        self.event_depth_edit.setValue(self.event_depth_edit.minimum())
        self.event_elev_edit.setValue(self.event_elev_edit.minimum())
        self.event_mag_edit.setValue(self.event_mag_edit.minimum())
        self.event_date_edit.setDate(self.event_date_edit.minimumDate())
        self.event_time_edit.setTime(self.event_time_edit.minimumTime())
        self.sigEventCleared.emit()
        #self.event_time_edit.
