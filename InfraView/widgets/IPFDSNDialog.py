import obspy
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.clients.fdsn.header import URL_MAPPINGS

import numpy as np

from PyQt5.QtWidgets import (QDialog, QDialogButtonBox, QWidget, QAbstractItemView, QLineEdit, QFormLayout,
                             QComboBox, QLabel, QVBoxLayout, QHBoxLayout,
                             QGroupBox, QPushButton, QDateEdit, QTimeEdit,
                             QSizePolicy, QSpinBox, QListWidget)

from PyQt5.QtCore import pyqtSignal, pyqtSlot, QDate, Qt

from InfraView.widgets import IPStationBrowser
from InfraView.widgets import IPUtils

class IPFDSNDialog(QDialog):

    fdsnWidget = None

    def __init__(self, parent):
        super(IPFDSNDialog, self).__init__(parent)
        self.buildUI()

    def buildUI(self):
        self.setWindowTitle(self.tr('InfraView: FDSN Import'))

        self.fdsnWidget = IPFDSNWidget()

        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(self.fdsnWidget)
        layout.addWidget(buttons)

        self.setLayout(layout)

        self.resize(400, self.height())

    def getStreams(self):
        # pass through to get the stream info
        return self.fdsnWidget.getStreams()

    def getInventory(self):
        # pass through to get the inventory info out
        return self.fdsnWidget.getInventory()


class IPNewFDSNDialog(QDialog):
    def __init__(self, parent):
        super().__init__()
        self.buildUI()

    def buildUI(self):
        self.setWindowTitle("InfraView: Add FDSN Service")
        form_layout = QFormLayout()
        name_label = QLabel("Service Name")
        self.service_name_edit = QLineEdit()
        url_label = QLabel("Service URL")
        self.service_url_edit = QLineEdit()
        self.service_url_edit.setPlaceholderText("ex: http://service.iris.edu")
        self.service_url_edit.setMinimumWidth(220)

        form_layout.addRow(name_label, self.service_name_edit)
        form_layout.addRow(url_label, self.service_url_edit)

         # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        main_layout = QVBoxLayout()
        main_layout.addLayout(form_layout)
        main_layout.addWidget(buttons)

        self.setLayout(main_layout)

    def get_service(self):
        return self.service_name_edit.text(), self.service_url_edit.text()
        

class IPFDSNWidget(QWidget):

    sigTracesReplaced = pyqtSignal(obspy.core.stream.Stream, obspy.core.inventory.inventory.Inventory)
    sigTracesAppended = pyqtSignal(obspy.core.stream.Stream, obspy.core.inventory.inventory.Inventory)

    stream = None
    inventory = None

    def __init__(self, parent=None):
        super().__init__()
        self.parent = parent
        self.buildUI()

    def buildUI(self):

        # Put together the options container
        formLayout = QFormLayout()
        optionsContainer = QWidget()
        optionsContainer.setLayout(formLayout)

        # in order to allow for custom fdsn servers, we have to make our own fdsn dictionary that we can append to
        self.fdsn_dictionary = URL_MAPPINGS
        #self.fdsn_dictionary.update({'BEER':'https://fdsnws.ilikebeer.com'})

        # First lets populate the client drop down
        self.cb = QComboBox()
        self.cb.setMinimumWidth(150)
        self.cb.currentTextChanged.connect(self.service_changed)
        label_service_name = QLabel(self.tr('Service:'))

        self.cb.addItems(self.fdsn_dictionary.keys())
        self.cb.setCurrentText('IRIS')
        self.cb.setToolTip(self.fdsn_dictionary['IRIS'])
        self.cb.currentIndexChanged[str].connect(self.onActivated_cb)

        # add button for new fdsn service
        self.new_service_button = QPushButton("+")
        self.new_service_button.setToolTip("Add an FDSN service")
        self.new_service_button.clicked.connect(self.add_service)

        service_layout = QHBoxLayout()
        service_layout.addWidget(self.cb)
        service_layout.addWidget(self.new_service_button)

        validator = IPUtils.CapsValidator(self)
        
        label_network_name = QLabel(self.tr('Network: '))
        self.networkNameBox = QLineEdit()
        self.networkNameBox.setMinimumWidth(170)
        self.networkNameBox.setToolTip('Wildcards OK \nCan be SEED network codes or data center defined codes. \nMultiple codes are comma-separated (e.g. "IU,TA").')
        self.networkNameBox.setValidator(validator)

        label_station_name = QLabel(self.tr('Station: '))
        self.stationNameBox = QLineEdit()
        self.stationNameBox.setMinimumWidth(170)
        self.stationNameBox.setToolTip('Wildcards OK \nOne or more SEED station codes. \nMultiple codes are comma-separated (e.g. "ANMO,PFO")')
        self.stationNameBox.setValidator(validator)

        label_location_str = QLabel(self.tr('Location:'))
        self.location_Box = QLineEdit('*')
        self.location_Box.setMinimumWidth(170)
        self.location_Box.setToolTip('Wildcards OK \nOne or more SEED location identifiers. \nMultiple identifiers are comma-separated (e.g. "00,01"). \nAs a special case “--“ (two dashes) will be translated to a string of two space characters to match blank location IDs.')
        self.location_Box.setValidator(validator)

        label_channel_str = QLabel(self.tr('Channel:'))
        self.channel_Box = QLineEdit('BDF')
        self.channel_Box.setMinimumWidth(170)
        self.channel_Box.setToolTip('Wildcards OK \nOne or more SEED channel codes. \nMultiple codes are comma-separated (e.g. "BHZ,HHZ")')
        self.channel_Box.setValidator(validator)

        label_startDate = QLabel(self.tr('Start Date (UTC):'))
        self.startDate_edit = QDateEdit()
        self.startDate_edit.setMinimumWidth(170)
        self.startDate_edit.setMinimumDate(QDate(1900, 1, 1))
        self.startDate_edit.setDisplayFormat('yyyy-MM-dd')
        self.startDate_edit.setDate(self.startDate_edit.minimumDate())

        label_startTime = QLabel(self.tr('Start Time (UTC):'))
        self.startTime_edit = QTimeEdit()
        self.startTime_edit.setMinimumWidth(170)
        self.startTime_edit.setDisplayFormat('HH:mm:ss.zzz')

        label_traceLength = QLabel(self.tr('Trace Length (s)'))
        self.traceLength_t = QSpinBox()
        self.traceLength_t.setMinimumWidth(170)
        self.traceLength_t.setMinimum(1)
        self.traceLength_t.setMaximum(999999999)
        self.traceLength_t.setValue(3600)

        replaceWaveButton = QPushButton('Replace')
        replaceWaveButton.clicked.connect(self.onClicked_replace)
        appendWaveButton = QPushButton('Append')
        appendWaveButton.clicked.connect(self.onClicked_append)

        self.stationListWidget = QListWidget()
        self.stationListWidget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.stationListWidget.itemSelectionChanged.connect(self.populateStationInfoFromStationList)

        self.browserButton = QPushButton('Station Browser')
        self.browserButton.clicked.connect(self.onClicked_browserButton)

        formLayout.addRow(label_service_name, service_layout)

        horizontalLineWidget = QWidget()
        horizontalLineWidget.setFixedHeight(2);
        horizontalLineWidget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed);
        horizontalLineWidget.setStyleSheet("background-color: #c0c0c0;");
        formLayout.addWidget(horizontalLineWidget)

        formLayout.addRow(label_network_name, self.networkNameBox)
        formLayout.addRow(label_station_name, self.stationNameBox)
        formLayout.addRow(label_location_str, self.location_Box)
        formLayout.addRow(label_channel_str, self.channel_Box)
        formLayout.addRow(label_startDate, self.startDate_edit)
        formLayout.addRow(label_startTime, self.startTime_edit)
        formLayout.addRow(label_traceLength, self.traceLength_t)

        horzLayout = QHBoxLayout(self)
        horzLayout.addWidget(replaceWaveButton)
        horzLayout.addWidget(appendWaveButton)
        addGroupBox = QGroupBox("Get Waveform(s)")
        addGroupBox.setLayout(horzLayout)

        vertlayout = QVBoxLayout(self)
        vertlayout.addWidget(optionsContainer)
        vertlayout.addWidget(addGroupBox)
        vertlayout.addWidget(self.stationListWidget)
        vertlayout.addWidget(self.browserButton)

        self.setLayout(vertlayout)

        # create dialogs here so that you only create it once, from here on you just run exec_() to make it pop up
        self.stationDialog = IPStationBrowser.IPStationDialog(self)
        self.add_serviceDialog = IPNewFDSNDialog(self)

    def add_service(self):
        if self.add_serviceDialog.exec_():
            name, url = self.add_serviceDialog.get_service()
            self.fdsn_dictionary[name] = url
            self.cb.addItem(name)
            self.cb.setCurrentText(name)
            self.cb.setToolTip(url)

    @pyqtSlot(str)
    def service_changed(self, name):
        self.cb.setToolTip(self.fdsn_dictionary[name])

    def onClicked_browserButton(self):

        if self.stationDialog.exec_(self.startDate_edit.date(),
                                    network=self.networkNameBox.text(),
                                    station=self.stationNameBox.text(),
                                    location=self.location_Box.text(),
                                    channel=self.channel_Box.text()):
            self.inventory = self.stationDialog.getInventory()

            inv_contents = self.inventory.get_contents()
            stationList = []

            # fill up the stationList with stations
            for item in inv_contents['stations']:
                stationList.append(item.strip())
            self.stationListWidget.clear()
            self.stationListWidget.addItems(stationList)

            if self.stationDialog.getStartDate() is not None:
                self.startDate_edit.setDate(self.stationDialog.getStartDate())

    # def populateWithEventInfo(self):
    #     if self.parent.eventWidget.hasValidEvent():
    #         self.currentEvent = self.parent.eventWidget.Dict()
    #     else:
    #         msgBox = QMessageBox()
    #         msgBox.setIcon(QMessageBox.Warning)
    #         msgBox.setText('There is not a valid event in the Event Tab')
    #         msgBox.setWindowTitle("Oops...")
    #         msgBox.exec_()
    #         return

    #     if self.currentEvent is not None:
    #         date = self.currentEvent['UTC Date']
    #         time = self.currentEvent['UTC Time'][0:5]

    #     qdate = QDate.fromString(date, 'yyyy-MM-dd')
    #     qtime = QTime.fromString(time, 'HH:mm')
    #     qtime.addSecs(-5*60)  # start about 5 minutes before event

    #     self.startDate_edit.setDate(qdate)
    #     self.startTime_edit.setTime(qtime)

    #     self.eventInfoPopulated = True

    # # if someone edits the event info, reflect the changes here
    # @QtCore.pyqtSlot(dict)
    # def updateEventInfo(self, event):
    #     if not self.eventInfoPopulated:
    #         return

    #     if self.parent.eventWidget.hasValidEvent():
    #         self.currentEvent = event

    #     if self.currentEvent is not None:
    #         date = event['UTC Date']
    #         time = event['UTC Time'][0:5]

    #     qdate = QDate.fromString(date, 'yyyy-MM-dd')
    #     qtime = QTime.fromString(time, 'HH:mm')
    #     qtime.addSecs(-5*60)  # start about 5 minutes before event

    #     self.startDate_edit.setDate(qdate)
    #     self.startTime_edit.setTime(qtime)

    def populateStationInfoFromStationList(self):
        items = self.stationListWidget.selectedItems()
        if len(items) < 1:
            return  # nothing to do

        netList = []
        staList = []

        for item in items:
            text = item.text()
            text = text.split(' ')[0]
            netSta = text.split('.')

            netList.append(netSta[0])
            staList.append(netSta[1])

            netList = list(set(netList))
            staList = list(set(staList))

            netString = ''
            for net in netList:
                netString = netString + net + ', '
            staString = ''
            for sta in staList:
                staString = staString + sta + ', '

        self.networkNameBox.setText(netString[:-2])
        self.stationNameBox.setText(staString[:-2])

    def onActivated_cb(self, key):
        if(key != 'choose...'):
            url = URL_MAPPINGS[key]

    def onClicked_replace(self):
        self.downloadWaveforms()
        if self.stream is not None and self.inventory is not None:
            self.sigTracesReplaced.emit(self.stream, self.inventory)

    def onClicked_append(self):
        self.downloadWaveforms()
        if self.stream is not None and self.inventory is not None:
            self.sigTracesAppended.emit(self.stream, self.inventory)

    @pyqtSlot(str, str)
    def add_custom_fdsn(self, name, url):
        self.fdsn_dictionary.update({name : url})
        # for brevity, lets just clear the combobox and repopulate it
        self.cb.clear()
        for key in sorted(self.fdsn_dictionary.keys()):
            self.cb.addItem(key)
        self.cb.setCurrentText(name)

    # get waveform button was clicked
    def downloadWaveforms(self):

        # get the inputs...inputs
        service = self.cb.currentText()
        if(service == 'choose...'):
            IPUtils.errorPopup('Please select a service to search')
            return
        print(self.fdsn_dictionary[service])
        client = Client(self.fdsn_dictionary[service])

        # Clear old streams because we don't need them anymore
        self.clearWaveforms()

        network = self.networkNameBox.text().upper().replace(' ', '')
        self.networkNameBox.setText(network)
        station = self.stationNameBox.text().upper().replace(' ', '')
        self.stationNameBox.setText(station)
        location = self.location_Box.text().upper().replace(' ', '')
        self.location_Box.setText(location)
        channel = self.channel_Box.text().upper().replace(' ', '')
        self.channel_Box.setText(channel)
        date = self.startDate_edit.date().toPyDate()
        time = self.startTime_edit.time().toPyTime()
        traceLength = self.traceLength_t.value()
        utcString = str(date) + 'T' + str(time)
        startTime = UTCDateTime(utcString)
        endTime = startTime + traceLength

        # Check for unfilled boxes
        if (network == '' or station == '' or channel == ''):
            IPUtils.errorPopup('You are missing some important info...\nNetwork, Station, Location, and Channel are all required data.')
            return



        # self.parent.setStatus('downloading Waveforms...')
        try:
            self.stream = client.get_waveforms(network, station, location, channel, startTime, endTime)

        except Exception:
            IPUtils.errorPopup('Failure loading waveform. \nDouble check that the values you entered are valid and the time and date are appropriate.')
            return

        for trace in self.stream:
            trace.data = trace.data - np.mean(trace.data)
        self.stream.merge(fill_value=0)

        # self.parent.setStatus('downloading Inventory...')
        # self.parent.consoleBox.append( 'Downloaded waveforms from '+service )

        # Now get the corresponding stations
        try:
            self.inventory = client.get_stations(network=network, station=station)
        except:
            IPUtils.errorPopup('Failure loading Inventory.  \nDouble check that the values you entered are valid and the time and date are appropriate.')
            return

    def getStreams(self):
        return self.stream

    def getInventory(self):
        return self.inventory

    def getService(self):
        return self.cb.currentText()

    def clearWaveforms(self):
        self.stream = None
        self.inventory = None

    def clear(self):

        self.clearWaveforms()

        self.stationListWidget.clear()
        self.cb.setCurrentText('IRIS')  # default to IRIS because why the hell not?
        self.networkNameBox.clear()
        self.stationNameBox.clear()
        self.location_Box.setText('*')
        self.channel_Box.setText('*')
        self.startDate_edit.setDate(self.startDate_edit.minimumDate())
        self.startTime_edit.setTime(self.startTime_edit.minimumTime())
        self.traceLength_t.setValue(3600)
