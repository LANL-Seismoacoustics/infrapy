import obspy
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.clients.fdsn.header import URL_MAPPINGS

import numpy as np

from PyQt5 import QtGui
from PyQt5.QtWidgets import (QWidget, QAbstractItemView, QLineEdit, QGridLayout,
                             QComboBox, QLabel, QVBoxLayout, QHBoxLayout,
                             QGroupBox, QPushButton, QDateEdit, QTimeEdit,
                             QSpinBox, QMessageBox, QListWidget)

from PyQt5.QtCore import pyqtSignal, pyqtSlot, QDate

from InfraView.widgets import IPStationBrowser

# since most of the fields will require capitalized values only, here is a validator for the
# lineEdits


class IPValidator(QtGui.QValidator):
    def validate(self, string, pos):
        return QtGui.QValidator.Acceptable, string.upper(), pos


class IPFDSNWidget(QWidget):

    sigTracesReplaced = pyqtSignal(obspy.core.stream.Stream, obspy.core.inventory.inventory.Inventory)
    sigTracesAppended = pyqtSignal(obspy.core.stream.Stream, obspy.core.inventory.inventory.Inventory)

    stream = None
    inventory = None

    def __init__(self, parent=None):
        super().__init__()
        self.parent = parent

        # this is for managing the event populated form elements
        # self.eventInfoPopulated = False
        # self.currentEvent = None

        self.__buildUI__()
        self.show()

    def __buildUI__(self):

        # Put together the options container
        gridLayout = QGridLayout()
        optionsContainer = QWidget()
        optionsContainer.setLayout(gridLayout)

        # First lets populate the client drop down
        self.cb = QComboBox()
        label_service_name = QLabel(self.tr('Service:'))

        fdsn_dictionary = URL_MAPPINGS
        fdsn_dictionary.update({'RaspShake':'https://fdsnws.rasberryshakedata.com'})

        for key in sorted(URL_MAPPINGS.keys()):
            self.cb.addItem(key)
        self.cb.setCurrentText('IRIS')
        self.cb.currentIndexChanged[str].connect(self.onActivated_cb)

        validator = IPValidator(self)
        label_network_name = QLabel(self.tr('Network: '))
        self.networkNameBox = QLineEdit()
        self.networkNameBox.setToolTip('Wildcards OK \nCan be SEED network codes or data center defined codes. \nMultiple codes are comma-separated (e.g. "IU,TA").')
        self.networkNameBox.setValidator(validator)

        label_station_name = QLabel(self.tr('Station: '))
        self.stationNameBox = QLineEdit()
        self.stationNameBox.setToolTip('Wildcards OK \nOne or more SEED station codes. \nMultiple codes are comma-separated (e.g. "ANMO,PFO")')
        self.stationNameBox.setValidator(validator)

        label_location_str = QLabel(self.tr('Location:'))
        self.location_Box = QLineEdit('*')
        self.location_Box.setToolTip('Wildcards OK \nOne or more SEED location identifiers. \nMultiple identifiers are comma-separated (e.g. "00,01"). \nAs a special case “--“ (two dashes) will be translated to a string of two space characters to match blank location IDs.')
        self.location_Box.setValidator(validator)

        label_channel_str = QLabel(self.tr('Channel:'))
        self.channel_Box = QLineEdit('*')
        self.channel_Box.setToolTip('Wildcards OK \nOne or more SEED channel codes. \nMultiple codes are comma-separated (e.g. "BHZ,HHZ")')
        self.channel_Box.setValidator(validator)

        label_startDate = QLabel(self.tr('Start Date (UTC):'))
        self.startDate_edit = QDateEdit()
        self.startDate_edit.setMinimumDate(QDate(1900, 1, 1))
        self.startDate_edit.setDisplayFormat('yyyy-MM-dd')
        self.startDate_edit.setDate(self.startDate_edit.minimumDate())

        label_startTime = QLabel(self.tr('Start Time (UTC):'))
        self.startTime_edit = QTimeEdit()
        self.startTime_edit.setDisplayFormat('HH:mm:ss.zzz')

        label_traceLength = QLabel(self.tr('Trace Length (s)'))
        self.traceLength_t = QSpinBox()
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

        gridLayout.addWidget(label_service_name, 0, 0)
        gridLayout.addWidget(self.cb, 0, 1)
        gridLayout.addWidget(label_network_name, 1, 0)
        gridLayout.addWidget(self.networkNameBox, 1, 1)
        gridLayout.addWidget(label_station_name, 2, 0)
        gridLayout.addWidget(self.stationNameBox, 2, 1)
        gridLayout.addWidget(label_location_str, 3, 0)
        gridLayout.addWidget(self.location_Box, 3, 1)
        gridLayout.addWidget(label_channel_str, 4, 0)
        gridLayout.addWidget(self.channel_Box, 4, 1)
        gridLayout.addWidget(label_startDate, 5, 0)
        gridLayout.addWidget(self.startDate_edit, 5, 1)
        gridLayout.addWidget(label_startTime, 6, 0)
        gridLayout.addWidget(self.startTime_edit, 6, 1)
        gridLayout.addWidget(label_traceLength, 7, 0)
        gridLayout.addWidget(self.traceLength_t, 7, 1)
        # gridLayout.addWidget(importEventButton, 8, 1, 1, 2)

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

        # create stationdialog here so that you only create it once, from here on you just run exec_() to make it pop up
        self.stationDialog = IPStationBrowser.IPStationDialog()

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

    # get waveform button was clicked
    def downloadWaveforms(self):
        # first check to make sure the boxes are filled...addWidget
        # TODO!!!

        # get the inputs...inputs
        service = self.cb.currentText()
        if(service == 'choose...'):
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText('Please select a service to search')
            msg.setWindowTitle('oops...')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
            return

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
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText('You are missing some important info...')
            msg.setWindowTitle('oops...')
            msg.setDetailedText("Network, Station, Location, and Channel are all required data.")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
            return

        # Download the waveforms
        # self.parent.setStatus('connecting to '+service)
        client = Client(service)

        # self.parent.setStatus('downloading Waveforms...')
        try:
            self.stream = client.get_waveforms(network, station, location, channel, startTime, endTime)

        except Exception:
            # self.parent.setStatus('')
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText('Failure loading waveform')
            msg.setWindowTitle('oops...')
            msg.setDetailedText("Double check that the values you entered are valid and the time and date are appropriate.")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
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
            # self.parent.setStatus('')
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText('Failure loading Inventory')
            msg.setWindowTitle('oops...')
            msg.setDetailedText("Double check that the values you entered are valid and the time and date are appropriate.")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()

            return
        # self.parent.consoleBox.append( 'Downloaded stations from '+service)
        # self.parent.setStatus('Finished...', 3000)

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
