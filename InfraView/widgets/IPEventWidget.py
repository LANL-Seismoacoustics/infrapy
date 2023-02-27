import os, json

from PyQt5.QtWidgets import (QDateEdit, QPushButton, QLabel, QLineEdit, QFileDialog, QHBoxLayout,
                             QFormLayout, QDoubleSpinBox, QVBoxLayout, QTimeEdit, QWidget,
                             QCheckBox)
from PyQt5.QtCore import QDir, QTime, QDate, QSettings, pyqtSignal
from PyQt5.QtGui import QIcon

from InfraView.widgets import IPEventBrowser

from obspy.core import UTCDateTime
import pyproj


class IPEventWidget(QWidget):

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

        formLayout = QFormLayout()
        self.showGT_cb = QCheckBox("Show Event on Map")
        self.showGT_cb.setChecked(False)

        self.displayEvent_cb = QCheckBox(self.tr('Show event line on waveform plots'))
        self.displayEvent_cb.setChecked(False)

        self.displayArrivals_cb = QCheckBox(self.tr('Show arrival estimations on waveform plots'))
        self.displayArrivals_cb.setChecked(False)
        self.displayArrivals_cb.setEnabled(self.displayEvent_cb.isChecked())
        
        label_event_name = QLabel(self.tr('Event ID: '))
        self.event_name_edit = QLineEdit()
        self.event_name_edit.setFixedWidth(100)

        label_latitude = QLabel(self.tr('Latitude (deg)'))
        self.event_lat_edit = QDoubleSpinBox()
        self.event_lat_edit.setFixedWidth(125)
        self.event_lat_edit.setRange(-90.1,90.0)
        self.event_lat_edit.setDecimals(8)
        self.event_lat_edit.setSingleStep(0.1)
        self.event_lat_edit.setValue(0.0)
        
        label_longitude = QLabel(self.tr('Longitude (deg)'))
        self.event_lon_edit = QDoubleSpinBox()
        self.event_lon_edit.setFixedWidth(125)
        self.event_lon_edit.setRange(-180.0, 180.0)
        self.event_lon_edit.setDecimals(8)
        self.event_lon_edit.setSingleStep(0.1)
        self.event_lon_edit.setValue(0.0)

        label_event_date = QLabel(self.tr('Date (UTC):'))
        self.event_date_edit = QDateEdit()
        self.event_date_edit.setFixedWidth(125)
        self.event_date_edit.setDisplayFormat("yyyy-MM-dd")

        label_event_time = QLabel(self.tr('Time (UTC):'))
        self.event_time_edit = QTimeEdit()
        self.event_time_edit.setFixedWidth(125)
        self.event_time_edit.setDisplayFormat('HH:mm:ss.zzz')        

        self.load_button = QPushButton(self.tr(' Load Event...'))
        self.load_button.setMaximumWidth(100)
        self.load_button.setIcon(self.openIcon)

        self.save_button = QPushButton(self.tr(' Save Event...'))
        self.save_button.setMaximumWidth(100)
        self.save_button.setIcon(self.saveAsIcon)

        self.clear_button = QPushButton(self.tr(' Clear'))
        self.save_button.setMaximumWidth(100)
        self.clear_button.setIcon(self.clearIcon)
        
        self.browse_button = QPushButton(self.tr(' IRIS Event Browser...'))
        self.browse_button.setMaximumWidth(200)
        
        show_layout = QVBoxLayout()
        show_layout.addWidget(self.showGT_cb)
        show_layout.addWidget(self.displayEvent_cb)
        show_layout.addWidget(self.displayArrivals_cb)

        show_layout_horiz = QHBoxLayout()
        show_layout_horiz.addStretch()
        show_layout_horiz.addLayout(show_layout)
        show_layout_horiz.addStretch()

        formLayout.addRow(label_event_name, self.event_name_edit)
        formLayout.addRow(label_longitude, self.event_lon_edit)
        formLayout.addRow(label_latitude, self.event_lat_edit)
        formLayout.addRow(label_event_time, self.event_time_edit)
        formLayout.addRow(label_event_date, self.event_date_edit)           
        
        form_hbox = QHBoxLayout()
        form_hbox.addStretch()
        form_hbox.addLayout(formLayout)
        form_hbox.addStretch()

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(self.load_button)
        buttonLayout.addWidget(self.save_button)
        buttonLayout.addWidget(self.clear_button)
        buttonLayout.addStretch()

        browseLayout = QHBoxLayout()
        browseLayout.addStretch()
        browseLayout.addWidget(self.browse_button)
        browseLayout.addStretch()

        verticalLayout = QVBoxLayout()
        verticalLayout.addLayout(show_layout_horiz)
        verticalLayout.addLayout(form_hbox)
        verticalLayout.addLayout(buttonLayout)
        verticalLayout.addLayout(browseLayout)
        verticalLayout.addStretch()

        self.setLayout(verticalLayout)

        self.eventBrowser = IPEventBrowser.IPEventDialog()

        self. connect_signals_and_slots()

        self.show()

    def connect_signals_and_slots(self):
        self.displayEvent_cb.stateChanged.connect(self.displayArrivals_cb.setEnabled)
        self.displayEvent_cb.stateChanged.connect(self.eventChanged)
        
        self.event_name_edit.textChanged.connect(self.eventChanged)

        self.event_lat_edit.valueChanged.connect(self.eventChanged)
        self.event_lon_edit.valueChanged.connect(self.eventChanged)

        self.event_date_edit.dateChanged.connect(self.eventChanged)
        self.event_time_edit.timeChanged.connect(self.eventChanged)

        self.load_button.clicked.connect(self.loadEvent)
        self.save_button.clicked.connect(self.saveEventAs)
        self.clear_button.clicked.connect(self.clear)
        self.browse_button.clicked.connect(self.browse)

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
        #self.displayArrivals_cb.setEnabled(self.displayEvent_cb.isChecked())
        self.sigEventWidgetChanged.emit(self.Dict())

    def Dict(self):
        # Returns a dictionary representation of the data in the widget
        lat = self.event_lat_edit.value()

        lon = self.event_lon_edit.value()

        eventdict = {'Name':self.event_name_edit.text(), 
                    'UTC Date':str(self.event_date_edit.date().toPyDate()), 
                    'UTC Time':str(self.event_time_edit.time().toPyTime()), 
                    'Longitude':lon, 
                    'Latitude':lat}

        return eventdict

    def saveEventAs(self):
        # pop up a save file dialog, default to project directory if a project is open, otherwise use the last used directory
        if self.window().getProject() is None:
            # force a new filename...
            savePath=self.settings.value("last_eventfile_directory", QDir.homePath())
        else:
            # There is an open project, so make the default save location correspond to what the project wants
            savePath = str(self.window().getProject().get_eventPath())

        self.savefile = QFileDialog.getSaveFileName(self, 'Save File', savePath)
        if self.savefile[0]:
            with open(self.savefile[0], 'w') as of:
                json.dump(self.Dict(), of, indent=4)

                if self.window().getProject() is None:
                    # if there is no open project, update the global settings 
                    self.settings.setValue("last_eventfile_directory", os.path.dirname(self.savefile[0]))

    def setEvent(self, event):
        # event is a dictionary containing the relevant information
        self.event_name_edit.setText(str(event['Name']))
        d = QDate(event['UTC Date'].year, event['UTC Date'].month, event['UTC Date'].day)
        self.event_date_edit.setDate(d)
        t = QTime(event['UTC Time'].hour, event['UTC Time'].minute, event['UTC Time'].second)
        self.event_time_edit.setTime(t)
        self.event_lat_edit.setValue(event['Latitude'])
        self.event_lon_edit.setValue(event['Longitude'])
        self.sigEventWidgetChanged.emit(event)

    def loadEvent(self):
        if self.window().getProject() is None:
            loadPath=self.settings.value("last_eventfile_directory", QDir.homePath())
        else:
            # There is an open project, so make the default save location correspond to what the project wants
            loadPath = str(self.window().getProject().get_eventPath())

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

            if self.window().getProject() is None:
                # if there is no open project, update the global settings 
                self.settings.setValue("last_eventfile_directory", os.path.dirname(self.__openfile[0]))

            self.sigEventWidgetChanged.emit(self.Dict())

    def calculate_arrival_travel_times(self, receiver_coord):
        # maybe this belongs in a different widget?
        if self.hasValidEvent():

            # first get event day/time in usable form
            eventTime = UTCDateTime(self.getUTCDateTimeString())
            
            # calculate the travel time in seconds for a pressure wave to go from event to receiver
            # receiver_coord is a tuple containing the (lat, lon) of the receiver
            apparent_vels = {'thermospheric': 250.0, 'stratospheric': 290.0, 'tropospheric': 340.0}      # all in m/s

            geod = pyproj.Geod(ellps='WGS84')
            _, _, distance = geod.inv(self.event_lon_edit.value(), self.event_lat_edit.value(), receiver_coord[1], receiver_coord[0])

            result_dict = {}
            for arrival, velocity in apparent_vels:
                result_dict[arrival] = distance/velocity
            return result_dict
        else:
            return None

    def browse(self):
        #self.eventDialog = IPEventBrowser.IPEventDialog()
        
        if self.eventBrowser.exec_():
            event = self.eventBrowser.getEvent()

            self.event_name_edit.setText(event['Name'])
            self.event_lon_edit.setValue(event['Longitude'])
            self.event_lat_edit.setValue(event['Latitude'])
            self.event_date_edit.setDate(QDate.fromString(event['UTC Date'], 'yyyy-MM-dd'))
            self.event_time_edit.setTime(QTime.fromString(event['UTC Time'], 'hh:mm:ss.zzz'))
            self.sigEventWidgetChanged.emit(self.Dict())

    def clear(self):
        self.event_name_edit.setText('')
        self.event_lon_edit.setValue(0.0)
        self.event_lat_edit.setValue(0.0)
        self.event_date_edit.setDate(QDate(2000,1,1))
        self.event_time_edit.setTime(QTime(00,00,00))
        self.sigEventCleared.emit()
