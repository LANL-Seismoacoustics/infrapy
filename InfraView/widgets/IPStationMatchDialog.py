from PyQt5.QtWidgets import (QApplication, QComboBox, QDialog, QLabel,
                             QListWidget, QListWidgetItem, QPushButton,
                             QVBoxLayout, QDialogButtonBox, QCheckBox)
from PyQt5.QtCore import Qt

from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.clients.fdsn.header import URL_MAPPINGS


class IPStationMatchDialog(QDialog):

    # This dialog is used to ATTEMPT to download stations for
    # current waveforms.

    new_inv = None

    def __init__(self, parent):
        super(IPStationMatchDialog, self).__init__(parent)
        self.parent = parent
        self.__buildUI__()

    def __buildUI__(self):

        self.setWindowTitle('InfraView: Reconcile Stations')

        blurb = QLabel(self.tr("This will ATTEMPT to download from an FDSN service the station info for the selected Stations. A '*' by the label indicates there is alread a station loaded. Uncheck any that you don't want overwritten."))
        blurb.setWordWrap(True)

        # First lets populate the client drop down
        self.cb = QComboBox()
        label_service_name = QLabel(self.tr('Service: '))
        for key in sorted(URL_MAPPINGS.keys()):
            self.cb.addItem(key)
        self.cb.setCurrentText('IRIS')

        self.statusLine = QLabel(' ')

        self.attemptButton = QPushButton(self.tr('Attempt to Download'))

        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
                                   Qt.Horizontal, self)
        buttons.button(QDialogButtonBox.Ok).setText('Add to Station List')
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        self.checkbox_layout = QVBoxLayout()

        layout = QVBoxLayout()
        layout.addWidget(blurb)
        layout.addWidget(label_service_name)
        layout.addWidget(self.cb)
        layout.addLayout(self.checkbox_layout)
        layout.addWidget(self.attemptButton)
        layout.addWidget(self.statusLine)
        layout.addWidget(buttons)

        self.setLayout(layout)

        self.connectSignalsandSlots()

    def exec_(self, channel_ids, existing_channels, timeRange):

        self.new_inv = None
        
        self.statusLine.setText('')

        # clear any existing checkboxes, and load new ones
        self.clear_checkboxes()

        for chan_id in channel_ids:
            if chan_id in existing_channels:
                chan_id += '*'

            self.checkbox_list.append(QCheckBox(chan_id, parent=self))
        
        for checkbox in self.checkbox_list:
            text = checkbox.text()
            if '*' in text:
                checkbox.setChecked(False)
            else:
                checkbox.setChecked(True)
            self.checkbox_layout.addWidget(checkbox)

        self.timeRange = timeRange

        return super().exec_()

    def connectSignalsandSlots(self):
        self.attemptButton.clicked.connect(self.download)

    def clear_checkboxes(self):
        self.checkbox_list = []
        # first lets clear out the label layout, and redraw it
        for i in reversed(range(self.checkbox_layout.count())): 
            self.checkbox_layout.removeWidget(self.checkbox_layout.itemAt(i).widget())

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
        for checkbox in self.checkbox_list:
            if checkbox.isChecked():
                traceID = checkbox.text()
                traceID.replace('*', '')    # remove the asterisk if it's there

                self.statusLine.setText('Attempting to download...' + traceID)
                QApplication.instance().processEvents()

                s = traceID.split('.')
                network = s[0]
                station = s[1]
                location = s[2]
                channel = s[3]

                print("net = {}   sta = {}   loc = {}   chan = {}".format(network, station, location, channel))

                inv = None
                try:
                    inv = client.get_stations(network=network, station=station, channel=channel, level='channel', startbefore=UTCDateTime(self.timeRange[0]), endafter=UTCDateTime(self.timeRange[0]))
                    foundCount += 1
                    checkbox.setStyleSheet('QCheckBox {color: green}')

                    print(inv)
                    
                except Exception as e:
                    checkbox.setStyleSheet('QCheckBox {color: red}')

                
                if inv is not None:
                    if self.new_inv is None:
                        self.new_inv = inv
                    else:
                        self.new_inv += inv

                
                self.statusLine.setText('Found ' + str(foundCount) + ' stations.')

    def getInventory(self):
        return self.new_inv
    
    



class IPListItem(QListWidgetItem):

    def __init__(self, text, parent=None):
        super(IPListItem, self).__init__(parent)
        self.setText(text)
