import  sys

from PyQt5.QtWidgets import (QApplication, QDateTimeEdit, QDialogButtonBox, QMessageBox, QPushButton, 
                             QLabel, QLineEdit, QGridLayout, QHBoxLayout, QDoubleSpinBox, 
                             QVBoxLayout, QWidget, QFrame, QAbstractItemView,
                             QTabWidget, QTableWidget, QTableWidgetItem, QDialog)

from PyQt5.QtCore import QDateTime, pyqtSlot, Qt

from InfraView.widgets import IPUtils

import obspy
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client

class IPEventBrowser(QWidget):

    def __init__(self):
        super().__init__()

        self.buildUI()


    def buildUI(self):

        self.eventID_edit = QLineEdit()

        self.startDateTime_edit = QDateTimeEdit()
        self.startDateTime_edit.setDisplayFormat('yyyy-MM-ddTHH:mm:ss.zzz')
        self.startDateTime_edit.setDateTime(self.startDateTime_edit.minimumDateTime())
        # self.startDateTime_edit.setSpecialValueText('yyyy-MM-ddTHH:mm:ss.zzz')

        self.endDateTime_edit = QDateTimeEdit()
        self.endDateTime_edit.setDisplayFormat('yyyy-MM-ddTHH:mm:ss.zzz')
        self.endDateTime_edit.setDateTime(self.endDateTime_edit.minimumDateTime())
        # self.endDateTime_edit.setSpecialValueText('yyyy-MM-ddTHH:mm:ss.zzz')

        self.lat_edit = QDoubleSpinBox()
        self.lat_edit.setRange(-90.1,90.0)  # the -90.1 is used as the "unset" value 
        self.lat_edit.setDecimals(8)
        self.lat_edit.setSpecialValueText('--')
        self.lat_edit.setSingleStep(0.1)
        self.lat_edit.setValue(self.lat_edit.minimum())

        self.lon_edit = QDoubleSpinBox()
        self.lon_edit.setRange(-180.1, 180.0)
        self.lon_edit.setDecimals(8)
        self.lon_edit.setSpecialValueText('--')
        self.lon_edit.setSingleStep(0.1)
        self.lon_edit.setValue(self.lon_edit.minimum())

        self.minLat_edit = QDoubleSpinBox()
        self.minLat_edit.setRange(-90.1,90.0)
        self.minLat_edit.setDecimals(8)
        self.minLat_edit.setSpecialValueText('--')
        self.minLat_edit.setSingleStep(0.1)
        self.minLat_edit.setValue(self.minLat_edit.minimum())

        self.maxLat_edit = QDoubleSpinBox()
        self.maxLat_edit.setRange(-90.1,90.0)
        self.maxLat_edit.setDecimals(8)
        self.maxLat_edit.setSpecialValueText('--')
        self.maxLat_edit.setSingleStep(0.1)
        self.maxLat_edit.setValue(self.maxLat_edit.minimum())

        self.minLon_edit = QDoubleSpinBox()
        self.minLon_edit.setRange(-180.1, 180.0)
        self.minLon_edit.setDecimals(8)
        self.minLon_edit.setSpecialValueText('--')
        self.minLon_edit.setSingleStep(0.1)
        self.minLon_edit.setValue(self.minLon_edit.minimum())

        self.maxLon_edit = QDoubleSpinBox()
        self.maxLon_edit.setRange(-180.1, 180.0)
        self.maxLon_edit.setDecimals(8)
        self.maxLon_edit.setSpecialValueText('--')
        self.maxLon_edit.setSingleStep(0.1)
        self.maxLon_edit.setValue(self.maxLon_edit.minimum())

        self.minRadius_edit = QDoubleSpinBox()
        self.minRadius_edit.setMinimum(0.0)
        self.minRadius_edit.setMaximum(180.0)
        self.minRadius_edit.setSpecialValueText('--')
        self.minRadius_edit.setSuffix(' deg')
        self.minRadius_edit.setValue(self.minRadius_edit.minimum())
        self.minRadius_edit.setEnabled(False)   #Need lat and lon before this means anything
        
        self.maxRadius_edit = QDoubleSpinBox()
        self.maxRadius_edit.setMinimum(0.0)
        self.maxRadius_edit.setMaximum(180.0)
        self.maxRadius_edit.setSpecialValueText('--')
        self.maxRadius_edit.setSuffix(' deg')
        self.maxRadius_edit.setValue(self.maxRadius_edit.minimum())
        self.maxRadius_edit.setEnabled(False)  #Need lat and lon before this means anything


        self.minDepth_edit = QDoubleSpinBox()
        self.minDepth_edit.setMinimum(0.0)
        self.minDepth_edit.setSpecialValueText('--')
        self.minDepth_edit.setSuffix(' km')
        self.minDepth_edit.setValue(self.minDepth_edit.minimum())

        self.maxDepth_edit = QDoubleSpinBox()
        self.maxDepth_edit.setMinimum(0.0)
        self.maxDepth_edit.setSpecialValueText('--')
        self.maxDepth_edit.setSuffix(' km')
        self.maxDepth_edit.setValue(self.maxDepth_edit.minimum())

        self.minMag_edit = QDoubleSpinBox()
        self.minMag_edit.setMinimum(0.0)
        self.minMag_edit.setSpecialValueText('--')
        self.minMag_edit.setValue(self.minMag_edit.minimum())

        self.maxMag_edit = QDoubleSpinBox()
        self.maxMag_edit.setMinimum(0.0)
        self.maxMag_edit.setSpecialValueText('--')
        self.maxMag_edit.setValue(self.maxMag_edit.minimum())

        self.searchButton = QPushButton('Search')
        self.searchButton.clicked.connect(self.downloadEvents)

        self.resetButton = QPushButton('Reset Inputs')
        self.resetButton.clicked.connect(self.resetWidget)

        self.resultsTable = QTableWidget()
        self.resultsTable.setSelectionBehavior(QAbstractItemView.SelectRows)

        main_layout = QVBoxLayout()

        layout = QGridLayout()

        layout.addWidget(QLabel(self.tr('Event ID: ')), 0 , 0, alignment=Qt.AlignRight)
        layout.addWidget(self.eventID_edit, 0, 1)

        layout.addWidget(self.HLine(), 1, 0, 1, 4) # This is a horizontal line

        layout.addWidget(QLabel(self.tr('Start Date/Time (UTC)')), 2, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.startDateTime_edit, 2, 1, 1, 2)

        layout.addWidget(QLabel(self.tr('End Date/Time (UTC)')), 3, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.endDateTime_edit, 3, 1, 1, 2)

        layout.addWidget(self.HLine(), 4, 0, 1, 4)  # This is a horizontal line

        layout.addWidget(QLabel(self.tr('Latitude: ')), 5, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.lat_edit, 5, 1)
        layout.addWidget(QLabel(self.tr('Longitude: ')), 5, 2, alignment=Qt.AlignRight)
        layout.addWidget(self.lon_edit, 5, 3)

        layout.addWidget(QLabel(self.tr('Minimum Latitude: ')), 6, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.minLat_edit, 6, 1)
        layout.addWidget(QLabel(self.tr('Maximum Latitude: ')), 6, 2, alignment=Qt.AlignRight)
        layout.addWidget(self.maxLat_edit, 6, 3)

        layout.addWidget(QLabel(self.tr('Minimum Longitude: ')), 7, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.minLon_edit, 7, 1)
        layout.addWidget(QLabel(self.tr('Maximum Longitude: ')), 7, 2, alignment=Qt.AlignRight)
        layout.addWidget(self.maxLon_edit, 7, 3)

        layout.addWidget(self.HLine(), 8, 0, 1, 4)  # This is a horizontal line

        layout.addWidget(QLabel(self.tr('Minimum Radius: ')), 9, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.minRadius_edit, 9, 1)
        layout.addWidget(QLabel(self.tr('Maximum Radius: ')), 9, 2, alignment=Qt.AlignRight)
        layout.addWidget(self.maxRadius_edit, 9, 3)

        layout.addWidget(QLabel(self.tr('Minimum Depth: ')), 10, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.minDepth_edit, 10, 1)
        layout.addWidget(QLabel(self.tr('Maximum Depth: ')), 10, 2, alignment=Qt.AlignRight)
        layout.addWidget(self.maxDepth_edit, 10, 3)

        layout.addWidget(QLabel(self.tr('Minimum Magnitude')), 11, 0, alignment=Qt.AlignRight)
        layout.addWidget(self.minMag_edit, 11, 1)
        layout.addWidget(QLabel(self.tr('Maximum Magnitude')), 11, 2, alignment=Qt.AlignRight)
        layout.addWidget(self.maxMag_edit, 11, 3)

        layout.addWidget(self.HLine(), 12, 0, 1, 4) # This is a horizontal line

        #self.spinner = QLabel()
        #self.spinner_movie = QMovie('../graphics/loader_small.gif')
        #self.spinner.setMovie(self.spinner_movie)
        #self.spinner.hide()

        layout.addWidget(self.resetButton, 13, 0)
        layout.addWidget(self.searchButton, 13, 2, 1, 2)

        centering_layout = QHBoxLayout()
        centering_layout.addStretch()
        centering_layout.addLayout(layout)
        centering_layout.addStretch()

        # placeholder for mapwidget, whenever I get around to making one
        self.mapWidget = QWidget()

        self.tabs = QTabWidget()
        self.tabs.addTab(self.resultsTable, 'Results')
        # self.tabs.addTab(self.mapWidget, 'Map')

        main_layout.addLayout(centering_layout)
        main_layout.addWidget(self.tabs)


        self.setLayout(main_layout)

        self.connectSignalsAndSlots()

        self.show()

    def connectSignalsAndSlots(self):
        self.lat_edit.valueChanged.connect(self.updateWidget)
        self.lon_edit.valueChanged.connect(self.updateWidget)

        self.startDateTime_edit.dateTimeChanged.connect(self.adjust_end_datetime)

    @pyqtSlot(QDateTime)
    def adjust_end_datetime(self, start_datetime: UTCDateTime):
        end_datetime = self.endDateTime_edit.dateTime()
        if start_datetime > end_datetime:
            self.endDateTime_edit.setDateTime(start_datetime)

    def getUTCDateTimeString(self, date: UTCDateTime, time: UTCDateTime) -> str:
        pdate = date.toPyDate()
        ptime = time.toPyTime()
        utcString = str(pdate) + 'T' + str(ptime)
        return utcString

    def downloadEvents(self): 

        # Iris is the default client, not sure if it would be useful to include others or not
        # Events come from the NEIC and the ISC
        try:
            client = Client("IRIS")
        except:
            IPUtils.errorPopup("failed to connect to client...proxy issue?")
            return None
        ####
        # laboriously gather the inputs
        if self.eventID_edit.text() == '':
            evt_id = None
        else:
            evt_id = self.eventID_edit.text()

        if self.startDateTime_edit.date() == self.startDateTime_edit.minimumDateTime():
            startDateTime = None
        else:
            date = self.startDateTime_edit.date()
            time = self.startDateTime_edit.time()
            sdatetime = self.getUTCDateTimeString(date,time)
            startDateTime = UTCDateTime(sdatetime)

        if self.endDateTime_edit.date() == self.endDateTime_edit.minimumDateTime():
            startDateTime = None
        else:
            date = self.endDateTime_edit.date()
            time = self.endDateTime_edit.time()
            sdatetime = self.getUTCDateTimeString(date, time)
            endDateTime = UTCDateTime(sdatetime)

        if self.lat_edit.value() == self.lat_edit.minimum():
            lat = None
        else:
            lat = self.lat_edit.value()

        if self.lon_edit.value() == self.lon_edit.minimum():
            lon = None
        else:
            lon = self.lon_edit.value()

        if self.minLat_edit.value() == self.minLat_edit.minimum():
            min_lat = None
        else:
            min_lat = self.minLat_edit.value()

        if self.maxLat_edit.value() == self.maxLat_edit.minimum():
            max_lat = None
        else:
            max_lat = self.maxLat_edit.value()

        if self.minLon_edit.value() == self.minLon_edit.minimum():
            min_lon = None
        else:
            min_lon = self.minLon_edit.value()

        if self.maxLon_edit.value() == self.maxLon_edit.minimum():
            max_lon = None
        else:
            max_lon = self.maxLon_edit.value()

        if self.lat_edit.value() != self.lat_edit.minimum() and self.lon_edit.value() != self.lon_edit.minimum():

            if self.minRadius_edit.value() == self.minRadius_edit.minimum():
                min_rad = None
            else:
                min_rad = self.minRadius_edit.value()

            if self.maxRadius_edit.value() == self.maxRadius_edit.minimum():
                max_rad = None
            else:
                max_rad = self.maxRadius_edit.value()
        else:
            min_rad = None
            max_rad = None


        if self.minDepth_edit.value() == self.minDepth_edit.minimum():
            minDep = None
        else:
            minDep = self.minDepth_edit.value()

        if self.maxDepth_edit.value() == self.maxDepth_edit.minimum():
            max_dep = None
        else:
            max_dep = self.maxDepth_edit.value()

        if self.minMag_edit.value() == self.minMag_edit.minimum():
            min_mag = None
        else:
            min_mag = self.minMag_edit.value()

        if self.maxMag_edit.value() == self.maxMag_edit.minimum():
            max_mag = None
        else:
            max_mag = self.maxMag_edit.value()
        
        try:
            self.cat = client.get_events(eventid=evt_id,
                                     starttime=startDateTime,
                                     endtime=endDateTime,
                                     latitude=lat,
                                     longitude=lon,
                                     minlatitude=min_lat,
                                     maxlatitude=max_lat,
                                     minlongitude=min_lon,
                                     maxlongitude=max_lon,
                                     minradius=min_rad,
                                     maxradius=max_rad,
                                     minmagnitude=min_mag,
                                     maxmagnitude=max_mag)
        except obspy.clients.fdsn.header.FDSNNoDataException:
            IPUtils.errorPopup('No data available')
            return

        self.populateresultsTable()
        
    def populateresultsTable(self):
        
        self.resultsTable.setRowCount(len(self.cat))
        self.resultsTable.setColumnCount(7)

        headers = ['Magnitude','UTC Date/Time', 'Latitude', 'Longitude', 'Depth', '']
        self.resultsTable.setHorizontalHeaderLabels(headers)

        for idx, item in enumerate(self.cat):
            for mag in item.magnitudes:
                if mag.magnitude_type is None:
                    mag.magnitude_type = ''
                self.resultsTable.setItem(idx, 0, QTableWidgetItem(str(mag.mag)+' '+mag.magnitude_type))
            for origin in item.origins:
                self.resultsTable.setItem(idx, 1, QTableWidgetItem(str(origin.time.isoformat())))
                self.resultsTable.setItem(idx, 2, QTableWidgetItem(str(origin.latitude)))
                self.resultsTable.setItem(idx, 3, QTableWidgetItem(str(origin.longitude)))
                self.resultsTable.setItem(idx, 4, QTableWidgetItem(str(origin.depth)))
                self.resultsTable.setItem(idx, 5, QTableWidgetItem(item.resource_id.id))
            #self.resultsTable.setData(item.origins[0].latitude)

        self.resultsTable.resizeColumnsToContents()

    def updateWidget(self):
        # We need lat and lon values if we are going to use the radius inputs
        enableRadius = self.lat_edit.value() != self.lat_edit.minimum() and self.lon_edit.value() != self.lon_edit.minimum()
        self.maxRadius_edit.setEnabled(enableRadius)
        self.minRadius_edit.setEnabled(enableRadius)

    def resetWidget(self):
        self.startDateTime_edit.setDateTime(self.startDateTime_edit.minimumDateTime())
        self.endDateTime_edit.setDateTime(self.endDateTime_edit.minimumDateTime())
        self.lat_edit.setValue(self.lat_edit.minimum())
        self.lon_edit.setValue(self.lon_edit.minimum())
        self.minLat_edit.setValue(self.minLat_edit.minimum())
        self.maxLat_edit.setValue(self.maxLat_edit.minimum())
        self.minLon_edit.setValue(self.minLon_edit.minimum())
        self.maxLon_edit.setValue(self.maxLon_edit.minimum())
        self.minRadius_edit.setValue(self.minRadius_edit.minimum())
        self.maxRadius_edit.setValue(self.maxRadius_edit.minimum())
        self.minDepth_edit.setValue(self.minDepth_edit.minimum())
        self.maxDepth_edit.setValue(self.maxDepth_edit.minimum())
        self.minMag_edit.setValue(self.minMag_edit.minimum())
        self.maxMag_edit.setValue(self.maxMag_edit.minimum())
        self.eventID_edit.setText('')

    def clearWidget(self):
        self.resetWidget()
        self.resultsTable.clearContents()

    def HLine(self):
        hl = QFrame()
        hl.setFrameShape(QFrame.HLine)
        hl.setFrameShadow(QFrame.Sunken)
        return hl

    def getEvent(self):
        event = {'Name':None, 
                 'UTC Date':None, 
                 'UTC Time':None, 
                 'Longitude':-180.1, 
                 'Latitude':-90.1, 
                 'Elevation':-1000000, 
                 'Depth':-1000000,
                 'Magnitude':0.0 }

        # returns a dictionary with the contents of the selected row
        row = self.resultsTable.currentRow()
        if row == -1:
            return None
        else:
            event['Magnitude'] = float(self.resultsTable.item(row, 0).text().split()[0])
            event['UTC Date'] = self.resultsTable.item(row, 1).text()[0:10]
            event['UTC Time'] = self.resultsTable.item(row, 1).text()[11:23]
            event['Latitude'] = float(self.resultsTable.item(row, 2).text())
            event['Longitude'] = float(self.resultsTable.item(row, 3).text())
            event['Depth'] = float(self.resultsTable.item(row, 4).text())
            event['Name'] = self.resultsTable.item(row, 5).text().split('=')[-1]
            event['Elevation'] = 0.0
        return event


class IPEventDialog(QDialog):

    def __init__(self, parent=None):
        super(IPEventDialog, self).__init__(parent)

        self.setWindowTitle('InfraView: Event Browser')
        self.eventBrowser = IPEventBrowser()
        layout = QVBoxLayout()
        layout.addWidget(self.eventBrowser)

        button_cancel = QPushButton(self.tr('Cancel'))
        button_ok = QPushButton('Export Selected')

        buttons = QDialogButtonBox()
        buttons.addButton(button_ok, QDialogButtonBox.AcceptRole)
        buttons.addButton(button_cancel, QDialogButtonBox.RejectRole)
        buttons.rejected.connect(self.reject)
        buttons.accepted.connect(self.myaccept)

        layout.addWidget(buttons)

        self.setLayout(layout)

    def myaccept(self):
        event = self.eventBrowser.getEvent()
        if event is None:
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Warning)
            msgBox.setText('No event selected.')
            msgBox.setWindowTitle("Oops...")
            msgBox.exec_()
            return
        else:
            super().accept()

    
    def getEvent(self):
        event = self.eventBrowser.getEvent()
        return event




if __name__ == '__main__':
    app = QApplication(sys.argv)

    widget = IPEventBrowser()
    sys.exit(app.exec_())


