import obspy
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import URL_MAPPINGS

from PyQt5 import QtGui, QtCore

from PyQt5.QtCore import Qt, QDate
from PyQt5.QtWidgets import (QApplication, QWidget,
                             QVBoxLayout, QLineEdit,
                             QGridLayout,
                             QLabel, QComboBox, QDialog,
                             QDialogButtonBox, 
                             QPushButton, QDateEdit,
                             QHBoxLayout, QDoubleSpinBox,
                             QListWidget, QAbstractItemView,
                             QTabWidget, QFrame)

from InfraView.widgets import IPUtils


class IPStationBrowser(QWidget):

    inventory = None

    def __init__(self, service=None, network=None):
        super().__init__()

        # these are for managing the event info in form elements
        self.eventInfoPopulated = False
        self.currentEvent = None

        self.__buildUI__(service, network)
        self.show()

    def __buildUI__(self, service, network):
        vertLayout = QVBoxLayout()
        self.setLayout(vertLayout)

        gridLayout = QGridLayout()

        # First lets populate the client drop down with all available services
        self.service_cb = QComboBox()
        # self.service_cb.addItem("choose...")
        fdsn_dictionary = URL_MAPPINGS
        fdsn_dictionary.update({'RaspShake':'https://fdsnws.rasberryshakedata.com'})

        for key in sorted(fdsn_dictionary.keys()):
            self.service_cb.addItem(key)
            self.service_cb.setCurrentText('IRIS')
        if service is not None:
            self.service_cb.setCurrentText(service)

        validator = IPUtils.CapsValidator(self)

        # Network selector
        self.network_edit = QLineEdit(network)
        self.network_edit.setToolTip('Wildcards OK \nCan be SEED network codes or data center defined codes. \nMultiple codes are comma-separated (e.g. "IU,TA").')
        self.network_edit.setValidator(validator)
        if network is not None:
            self.network_edit.setText(network)

        # Station selector
        self.station_edit = QLineEdit()
        self.station_edit.setToolTip('Wildcards OK \nOne or more SEED station codes. \nMultiple codes are comma-separated (e.g. "ANMO,PFO")')
        self.station_edit.setValidator(validator)

        # Location selector
        self.location_edit = QLineEdit()
        self.location_edit.setToolTip('One or more SEED location identifiers. \nMultiple identifiers are comma-separated (e.g. "00,01"). \nAs a special case “--“ (two dashes) will be translated to a string of two space characters to match blank location IDs.')
        self.location_edit.setText('')

        # Channel selector
        self.channel_edit = QLineEdit()
        self.channel_edit.setToolTip('Wildcards OK \nOne or more SEED channel codes. \nMultiple codes are comma-separated (e.g. "BHZ,HHZ")')
        self.channel_edit.setValidator(validator)

        # Date time editors
        minimumDate = QDate(1900, 1, 1)
        self.startDate_edit = QDateEdit()
        self.startDate_edit.setDisplayFormat('yyyy-MM-dd')
        self.startDate_edit.setMinimumDate(minimumDate)
        self.startDate_edit.setDate(self.startDate_edit.minimumDate())

        self.endDate_edit = QDateEdit()
        self.endDate_edit.setDisplayFormat('yyyy-MM-dd')
        self.endDate_edit.setMinimumDate(minimumDate)
        self.endDate_edit.setDate(self.endDate_edit.minimumDate())

        self.lat_edit = QDoubleSpinBox()
        self.lat_edit.setRange(-90.1, 90.0)  # the -90.1 is used as the "unset" value
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

        self.evt_import_button = QPushButton('Populate with Event Info')

        self.minLat_edit = QDoubleSpinBox()
        self.minLat_edit.setRange(-90.1, 90.0)
        self.minLat_edit.setDecimals(8)
        self.minLat_edit.setSpecialValueText('--')
        self.minLat_edit.setSingleStep(0.1)
        self.minLat_edit.setValue(self.minLat_edit.minimum())

        self.maxLat_edit = QDoubleSpinBox()
        self.maxLat_edit.setRange(-90.1, 90.0)
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
        self.minRadius_edit.setEnabled(False)   # Need lat and lon before this means anything

        self.maxRadius_edit = QDoubleSpinBox()
        self.maxRadius_edit.setMinimum(0.0)
        self.maxRadius_edit.setMaximum(180.0)
        self.maxRadius_edit.setSpecialValueText('--')
        self.maxRadius_edit.setSuffix(' deg')
        self.maxRadius_edit.setValue(self.maxRadius_edit.minimum())
        self.maxRadius_edit.setEnabled(False)  # Need lat and lon before this means anything

        # Set up the search button here
        self.searchButton = QPushButton('Search')

        # Reset Button
        self.resetButton = QPushButton('Reset Inputs')

        # assemble the grid layout
        gridLayout.addWidget(QLabel(self.tr('Service: ')), 0, 0, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.service_cb, 0, 1)

        gridLayout.addWidget(QLabel(self.tr('Network(s): ')), 1, 0, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.network_edit, 1, 1)
        gridLayout.addWidget(QLabel(self.tr('Station(s): ')), 1, 2, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.station_edit, 1, 3)

        gridLayout.addWidget(QLabel(self.tr('Location(s): ')), 2, 0, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.location_edit, 2, 1)
        gridLayout.addWidget(QLabel(self.tr('Channel(s): ')), 2, 2, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.channel_edit, 2, 3)

        gridLayout.addWidget(self.HLine(), 3, 0, 1, 4)  # This is a horizontal line

        gridLayout.addWidget(QLabel(self.tr('Start Date (UTC)')), 4, 0, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.startDate_edit, 4, 1)

        gridLayout.addWidget(QLabel(self.tr('End Date (UTC)')), 4, 2, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.endDate_edit, 4, 3)

        gridLayout.addWidget(QLabel(self.tr('Latitude: ')), 5, 0, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.lat_edit, 5, 1)
        gridLayout.addWidget(QLabel(self.tr('Longitude: ')), 5, 2, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.lon_edit, 5, 3)

        gridLayout.addWidget(self.evt_import_button, 6, 1, 1, 2)

        gridLayout.addWidget(self.HLine(), 7, 0, 1, 4)  # This is a horizontal line

        gridLayout.addWidget(QLabel(self.tr('Minimum Latitude: ')), 8, 0, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.minLat_edit, 8, 1)
        gridLayout.addWidget(QLabel(self.tr('Maximum Latitude: ')), 8, 2, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.maxLat_edit, 8, 3)

        gridLayout.addWidget(QLabel(self.tr('Minimum Longitude: ')), 9, 0, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.minLon_edit, 9, 1)
        gridLayout.addWidget(QLabel(self.tr('Maximum Longitude: ')), 9, 2, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.maxLon_edit, 9, 3)

        gridLayout.addWidget(self.HLine(), 10, 0, 1, 4)  # This is a horizontal line

        gridLayout.addWidget(QLabel(self.tr('Minimum Radius: ')), 11, 0, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.minRadius_edit, 11, 1)
        gridLayout.addWidget(QLabel(self.tr('Maximum Radius: ')), 11, 2, alignment=Qt.AlignRight)
        gridLayout.addWidget(self.maxRadius_edit, 11, 3)

        gridLayout.addWidget(self.HLine(), 12, 0, 1, 4)  # This is a horizontal line

        gridLayout.addWidget(self.resetButton, 13, 0)
        gridLayout.addWidget(self.searchButton, 13, 2, 1, 2)

        # Center the form in the widget
        centeringLayout = QHBoxLayout()
        centeringLayout.addStretch()
        centeringLayout.addLayout(gridLayout)
        centeringLayout.addStretch()

        self.stationListWidget = QListWidget()
        # self.stationListWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        # for multiple selection, uncomment out the next line, and comment out the above line
        self.stationListWidget.setSelectionMode(QAbstractItemView.ExtendedSelection)

        # create a tab widget to hold the stationList and the station map
        stationTabs = QTabWidget()
        stationTabs.addTab(self.stationListWidget, "Station List")

        # Placeholder for map widget, whenever I get around to making one
        self.mapWidget = QWidget()
        stationTabs.addTab(self.mapWidget, 'Map')

        vertLayout.addLayout(centeringLayout)

        vertLayout.addWidget(stationTabs)

        self.connectSignalsAndSlots()

    def connectSignalsAndSlots(self):

        self.searchButton.clicked.connect(self.onActivated_search)
        self.resetButton.clicked.connect(self.onActivated_reset)
        self.evt_import_button.clicked.connect(self.populateEventInfo)

        self.stationListWidget.itemClicked.connect(self.getSelected)

        self.service_cb.currentIndexChanged[str].connect(self.serviceClicked)

        self.lat_edit.valueChanged.connect(self.updateWidget)
        self.lon_edit.valueChanged.connect(self.updateWidget)

        self.startDate_edit.dateChanged.connect(self.adjust_end_date)

    @QtCore.pyqtSlot(QDate)
    def adjust_end_date(self, start_date):
        if self.endDate_edit.date() == self.endDate_edit.minimumDate():
            self.endDate_edit.setDate(start_date)
        elif self.endDate_edit.date() < start_date:
            self.endDate_edit.setDate(start_date)

    def populateStationInfo(self):

        service = self.service_cb.currentText()

        if self.network_edit.text() == '':
            network = None
        else:
            network = self.network_edit.text()

        if self.station_edit.text() == '':
            station = None
        else:
            station = self.station_edit.text()

        if self.location_edit.text() == '':
            location = None
        else:
            location = self.location_edit.text()

        if self.channel_edit.text() == '':
            channel = None
        else:
            channel = self.channel_edit.text()

        if self.lat_edit.value() == self.lat_edit.minimum():
            lat = None
        else:
            lat = self.lat_edit.value()

        if self.lon_edit.value() == self.lon_edit.minimum():
            lon = None
        else:
            lon = self.lon_edit.value()

        if self.minLat_edit.value() == self.minLat_edit.minimum():
            minLat = None
        else:
            minLat = self.minLat_edit.value()

        if self.maxLat_edit.value() == self.maxLat_edit.minimum():
            maxLat = None
        else:
            maxLat = self.maxLat_edit.value()

        if self.minLon_edit.value() == self.minLon_edit.minimum():
            minLon = None
        else:
            minLon = self.minLon_edit.value()

        if self.maxLon_edit.value() == self.maxLon_edit.minimum():
            maxLon = None
        else:
            maxLon = self.maxLon_edit.value()

        if self.lat_edit.value() != self.lat_edit.minimum() and self.lon_edit.value() != self.lon_edit.minimum():

            if self.minRadius_edit.value() == self.minRadius_edit.minimum():
                minRad = None
            else:
                minRad = self.minRadius_edit.value()

            if self.maxRadius_edit.value() == self.maxRadius_edit.minimum():
                maxRad = None
            else:
                maxRad = self.maxRadius_edit.value()
        else:
            minRad = None
            maxRad = None

        if service != '':
            try:
                client = Client(service)
            except obspy.clients.fdsn.header.FDSNException as e:
                IPUtils.errorPopup(str(e))
                return

            startDate = self.startDate_edit.date().toString("yyyy-MM-dd")
            endDate = self.endDate_edit.date().toString("yyyy-MM-dd")

            try:
                self.inventory = client.get_stations(network=network,
                                                 station=station,
                                                 location=location,
                                                 channel=channel,
                                                 starttime=startDate,
                                                 endtime=endDate,
                                                 latitude=lat,
                                                 longitude=lon,
                                                 minlongitude=minLon,
                                                 maxlongitude=maxLon,
                                                 minlatitude=minLat,
                                                 maxlatitude=maxLon,
                                                 minradius=minRad,
                                                 maxradius=maxRad,
                                                 format='xml')
            except:
                emptyMessage = ["...No Data Found..."]
                self.stationListWidget.clear()
                self.stationListWidget.addItems(emptyMessage)
                return

            self.stationListWidget.clear()

            inv_contents = self.inventory.get_contents()

            stationList = []

            # fill up the stationList with stations
            for item in inv_contents['stations']:
                stationList.append(item.strip())
            self.stationListWidget.clear()
            self.stationListWidget.addItems(stationList)

    def onActivated_search(self):
        self.stationListWidget.clear()
        QtGui.QGuiApplication.processEvents()
        QApplication.processEvents()
        msg = ['...Searching...']
        self.stationListWidget.addItems(msg)
        QApplication.processEvents()
        QtGui.QGuiApplication.processEvents()
        self.populateStationInfo()
        return

    def onActivated_reset(self):
        self.inventory = None
        self.network_edit.setText('')
        self.station_edit.setText('')
        self.location_edit.setText('')
        self.channel_edit.setText('')
        self.startDate_edit.setDate(self.startDate_edit.minimumDate())
        self.endDate_edit.setDate(self.endDate_edit.minimumDate())
        self.lat_edit.setValue(self.lat_edit.minimum())
        self.lon_edit.setValue(self.lon_edit.minimum())
        self.minLat_edit.setValue(self.minLat_edit.minimum())
        self.minLon_edit.setValue(self.minLon_edit.minimum())
        self.maxLat_edit.setValue(self.maxLat_edit.minimum())
        self.maxLon_edit.setValue(self.maxLon_edit.minimum())
        self.minRadius_edit.setValue(self.minRadius_edit.minimum())
        self.maxRadius_edit.setValue(self.maxRadius_edit.minimum())
        self.stationListWidget.clear()
        self.eventInfoPopulated = False

    def populateEventInfo(self):
        # find out if there is a valid event in the eventwidget, if there is, proceed
        #lol...there's got to be a better way
        event_widget = self.parentWidget().parentWidget().parentWidget().parentWidget().eventWidget
        if event_widget.hasValidEvent():
            event_dict = event_widget.Dict()
            date = event_dict['UTC Date']
            lat = event_dict['Latitude']
            lon = event_dict['Longitude']
            # populate the form elements
            self.lat_edit.setValue(lat)
            self.lon_edit.setValue(lon)
            datedate = QDate.fromString(date, 'yyyy-MM-dd')
            self.startDate_edit.setDate(datedate)
            self.endDate_edit.setDate(datedate.addDays(1))
            self.eventInfoPopulated = True
        else:
            IPUtils.errorPopup('There is not a valid event in the Event Tab')

    # def setEvent(self, event):
    #     self.currentEvent = event

    #     # Update the event info if it has been imported
    #     if self.eventInfoPopulated:
    #         self.populateEventInfo()

    def updateWidget(self):
        # We need lat and lon values if we are going to use the radius inputs
        enableRadius = self.lat_edit.value() != self.lat_edit.minimum() and self.lon_edit.value() != self.lon_edit.minimum()
        self.maxRadius_edit.setEnabled(enableRadius)
        self.minRadius_edit.setEnabled(enableRadius)

    def getSelected(self):
        selectedText = self.stationListWidget.currentItem().text()

        # first split the string at the first space, this will return just the
        # network and the station
        name = selectedText.split(" ")
        fullID = name[0]
        # next split that at the period, to seperate the network and station
        networkID = fullID.split(".")[0]
        stationID = fullID.split(".")[1]

        result = {"networkID": networkID, "stationID": stationID}

        return result

    def get_current_center(self):
        # this method will calculate the center of the current inventory and will return a [lat,lon]

        # TODO: This is not really setup right now to handle the (very rare) case where an array straddles the
        # international date line

        lat, lon, cnt = 0, 0, 0

        for network in self.inventory:
            for station in network:
                lat += station.latitude
                lon += station.longitude
                cnt += 1

        return [lat / cnt, lon / cnt]

    def getInventory(self):
        return self.inventory

    def getSelectedInventory(self):
        reducedInv = []
        for item in self.inventory:
            if item.isSelected():
                reducedInv.append(item)
        return reducedInv

    def serviceClicked(self):
        self.stationListWidget.clear()
        self.network_edit.clear()

    def HLine(self):
        hl = QFrame()
        hl.setFrameShape(QFrame.HLine)
        hl.setFrameShadow(QFrame.Sunken)
        return hl


class IPStationDialog(QDialog):

    # I pass event to the Dialog in case the user wants to populate the date/lat/lon fields with its info
    def __init__(self, parent=None):
        super(IPStationDialog, self).__init__(parent)

        self.setWindowTitle('InfraView: Station Browser')
        self.stationBrowser = IPStationBrowser()
        layout = QVBoxLayout()
        layout.addWidget(self.stationBrowser)

        button_cancel = QPushButton(self.tr('Cancel'))
        button_ok = QPushButton('Export List')

        buttons = QDialogButtonBox()
        buttons.addButton(button_ok, QDialogButtonBox.AcceptRole)
        buttons.addButton(button_cancel, QDialogButtonBox.RejectRole)
        buttons.rejected.connect(self.reject)
        buttons.accepted.connect(self.myaccept)

        layout.addWidget(buttons)

        self.setLayout(layout)

    def exec_(self, start_date=None, network=None, station=None, location=None, channel=None):
        if start_date is not None:
            self.stationBrowser.startDate_edit.setDate(start_date)

        if network is not None:
            self.stationBrowser.network_edit.setText(network)

        if station is not None:
            self.stationBrowser.station_edit.setText(station)

        if location is not None:
            self.stationBrowser.location_edit.setText(location)

        if channel is not None:
            self.stationBrowser.channel_edit.setText(channel)

        return super().exec_()

    def myaccept(self):
        inv = self.stationBrowser.getInventory()
        if inv is None:
            IPUtils.errorPopup('No station(s) selected.')
            return
        else:
            super().accept()

    # this pipes the inventory out of the browser widget
    def getInventory(self):
        return self.stationBrowser.getInventory()

    # this pipes the start date out of the browser widget
    def getStartDate(self):
        if self.stationBrowser.startDate_edit.date() == self.stationBrowser.startDate_edit.minimumDate():
            return None
        else:
            return self.stationBrowser.startDate_edit.date()

    # this just pipes an event into the browser widget
    # @QtCore.pyqtSlot()
    # def setEvent(self, event):
    #     self.stationBrowser.setEvent(event)

