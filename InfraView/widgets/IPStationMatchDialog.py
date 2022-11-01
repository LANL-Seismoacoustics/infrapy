from PyQt5.QtWidgets import (QApplication, QComboBox, QDialog, QLabel,
                             QListWidget, QListWidgetItem, QPushButton,
                             QVBoxLayout, QDialogButtonBox)
from PyQt5.QtCore import Qt

from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.clients.fdsn.header import URL_MAPPINGS


class IPStationMatchDialog(QDialog):

    # This dialog is used to ATTEMPT to download stations (that aren't already downloaded) for
    # currently visible waveforms.

    inv = None

    def __init__(self, parent=None):
        super(IPStationMatchDialog, self).__init__(parent)

        self.__buildUI__()

    def __buildUI__(self):

        self.setWindowTitle('InfraView: Reconcile Stations')

        blurb = QLabel(self.tr('This will ATTEMPT to download the station info for the following Stations.'))
        blurb.setWordWrap(True)

        # First lets populate the client drop down
        self.cb = QComboBox()
        label_service_name = QLabel(self.tr('Service: '))
        for key in sorted(URL_MAPPINGS.keys()):
            self.cb.addItem(key)
        self.cb.setCurrentText('IRIS')

        self.stationListEdit = QListWidget()
        self.stationListEdit.setMinimumWidth(300)

        self.statusLine = QLabel(' ')

        self.attemptButton = QPushButton(self.tr('Attempt to Download'))

        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
                                   Qt.Horizontal, self)
        buttons.button(QDialogButtonBox.Ok).setText('Add to Station List')
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(blurb)
        layout.addWidget(label_service_name)
        layout.addWidget(self.cb)
        layout.addWidget(self.stationListEdit)
        layout.addWidget(self.attemptButton)
        layout.addWidget(self.statusLine)
        layout.addWidget(buttons)

        self.setLayout(layout)

        self.connectSignalsandSlots()

    def exec_(self, stationIDs, timeRange):

        self.stationListEdit.clear()
        self.statusLine.setText('')

        for stationID in stationIDs:
            self.stationListEdit.addItem(IPListItem(stationID))

        self.timeRange = timeRange

        return super().exec_()

    def connectSignalsandSlots(self):
        self.attemptButton.clicked.connect(self.download)

    def download(self):
        self.statusLine.setText('Connecting...')
        QApplication.instance().processEvents()

        service = self.cb.currentText()
        try:
            client = Client(service)
            self.statusLine.setText('Connected')
        except:
            self.statusLine.setText('Service unreachable')
        foundCount = 0

        for i in range(0, self.stationListEdit.count()):
            listItem = self.stationListEdit.item(i)
            traceID = listItem.text()

            self.statusLine.setText('Attempting to download...' + traceID)
            QApplication.instance().processEvents()

            s = traceID.split('.')
            network = s[0]
            station = s[1]
            # if s[2] == '':
            location = None
            # else:
            #     location = s[2]

            try:
                inv = client.get_stations(network=network, station=station, location=location, startbefore=UTCDateTime(self.timeRange[0]), endafter=UTCDateTime(self.timeRange[1]))
                listItem.setForeground(Qt.darkGreen)
                if self.inv is None:
                    self.inv = inv
                else:
                    self.inv += inv
                foundCount = foundCount + 1
            except:
                listItem.setForeground(Qt.red)

        self.statusLine.setText('Found ' + str(foundCount) + ' of ' + str(self.stationListEdit.count()) + ' stations.')

    def getInventory(self):
        return self.inv


class IPListItem(QListWidgetItem):

    def __init__(self, text, parent=None):
        super(IPListItem, self).__init__(parent)
        self.setText(text)
