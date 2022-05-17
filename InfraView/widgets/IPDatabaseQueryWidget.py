from PyQt5.QtWidgets import (QDateEdit, QFormLayout, QFrame, QHBoxLayout, QLabel, QLineEdit, 
                             QMessageBox, QPushButton, QSpinBox, QTimeEdit, QVBoxLayout,
                             QPlainTextEdit, QSizePolicy)
from PyQt5.QtCore import QDate, pyqtSlot

from obspy.core.utcdatetime import UTCDateTime

from infrapy.utils import database

from InfraView.widgets import IPFDSNWidget

class IPEventQueryWidget(QFrame):
    def __init__(self, parent):
        super().__init__()
        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        size_policy = self.sizePolicy()
        size_policy.setHorizontalPolicy(QSizePolicy.Fixed)
        self.setSizePolicy(size_policy)
        self.parent = parent
        self.session = None
        self.query_string = ""
        self.buildUI()

    def buildUI(self):
        validator = IPFDSNWidget.IPValidator()

        title_label = QLabel("\tEvent Query")
        title_label.setStyleSheet("QLabel {font-weight:bold; color: white; background-color: black}")

        self.evid_edit = QLineEdit()
        self.evid_edit.setMaximumWidth(150)
        self.evid_edit.setToolTip('')

        self.event_name_edit = QLineEdit()
        self.event_name_edit.setMaximumWidth(150)
        self.event_name_edit.setToolTip("Wildcards ok, case sensitive")

        form_layout1 = QFormLayout()
        form_layout1.addRow("EVID: ", self.evid_edit)
        form_layout1.addRow("Event Name: ", self.event_name_edit)

        self.query_textEdit = QPlainTextEdit()
        self.query_textEdit.setMaximumHeight(120)
        self.query_textEdit.setPlaceholderText("Query string...")

        self.clear_button = QPushButton("Clear")
        self.query_button = QPushButton("Query Events")

        # this bit centers the query button...
        horiz_layout = QHBoxLayout()
        horiz_layout.addStretch()
        horiz_layout.addWidget(self.clear_button)
        horiz_layout.addWidget(self.query_button)
        horiz_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addWidget(title_label)
        main_layout.addLayout(form_layout1)
        main_layout.addStretch()
        main_layout.addWidget(self.query_textEdit)
        main_layout.addLayout(horiz_layout)
        self.setLayout(main_layout)


class IPDatabaseQueryWidget(QFrame):

    def __init__(self, parent):
        super().__init__()
        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        size_policy = self.sizePolicy()
        size_policy.setHorizontalPolicy(QSizePolicy.Fixed)
        self.setSizePolicy(size_policy)
        self.parent = parent 
        self.session = None
        self.query_string = ""
        self.buildUI()

    def buildUI(self):

        validator = IPFDSNWidget.IPValidator()

        title_label = QLabel("\tWaveform Query")
        title_label.setStyleSheet("QLabel {font-weight:bold; color: white; background-color: black}")

        self.net_edit = QLineEdit()
        self.net_edit.setMaximumWidth(150)
        self.net_edit.setToolTip('Wildcards OK \nCan be SEED network codes or data center defined codes. \nMultiple codes are comma-separated (e.g. "IU,TA").')
        self.net_edit.setValidator(validator)

        self.sta_edit = QLineEdit()
        self.sta_edit.setMaximumWidth(150)
        self.sta_edit.setToolTip('Wildcards OK \nOne or more SEED station codes. \nMultiple codes are comma-separated (e.g. "ANMO,PFO")')
        self.sta_edit.setValidator(validator)

        self.cha_edit = QLineEdit()
        self.cha_edit.setMaximumWidth(150)
        self.cha_edit.setToolTip('Wildcards OK \nOne or more SEED channel codes. \nMultiple codes are comma-separated (e.g. "BHZ,HHZ")')
        self.cha_edit.setValidator(validator)

        self.start_date_edit = QDateEdit()
        self.start_date_edit.setMaximumWidth(150)
        self.start_date_edit.setDisplayFormat('yyyy-MM-dd')
        self.start_date_edit.setMinimumDate(QDate(1900, 1, 1))
        self.start_date_edit.setDate(self.start_date_edit.minimumDate())

        self.start_time_edit = QTimeEdit()
        self.start_time_edit.setMaximumWidth(150)
        self.start_time_edit.setDisplayFormat("HH:mm:ss")

        self.duration_edit = QSpinBox()
        self.duration_edit.setMinimum(1)
        self.duration_edit.setMaximum(999999)
        self.duration_edit.setMaximumWidth(150)
        self.duration_edit.setValue(600)

        form_layout1 = QFormLayout()
        form_layout1.addRow("Network: ", self.net_edit)
        form_layout1.addRow("Station: ", self.sta_edit)
        form_layout1.addRow("Channel: ", self.cha_edit)
        form_layout1.addRow("Start date: ", self.start_date_edit)
        form_layout1.addRow("Start time(UTC): ", self.start_time_edit)
        form_layout1.addRow("Duration (s): ", self.duration_edit)

        self.query_textEdit = QPlainTextEdit()
        self.query_textEdit.setMaximumHeight(120)
        self.query_textEdit.setPlaceholderText("Query string...")

        self.clear_button = QPushButton("Clear")
        self.query_button = QPushButton("Query Database")

        # this bit centers the query button...
        horiz_layout = QHBoxLayout()
        horiz_layout.addStretch()
        horiz_layout.addWidget(self.clear_button)
        horiz_layout.addWidget(self.query_button)
        horiz_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addWidget(title_label)
        main_layout.addLayout(form_layout1)
        main_layout.addStretch()
        main_layout.addWidget(self.query_textEdit)
        main_layout.addLayout(horiz_layout)
        self.setLayout(main_layout)

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.query_button.pressed.connect(self.query_database)

    def get_current_session(self):
        return self.parent.ipdatabase_connect_widget.session

    def query_database(self):

        session = self.get_current_session()
        if session is None:
            self.errorPopup("No current active session")
            return

        # get the start date and time into a UTCDateTime
        date = self.start_date_edit.date().toPyDate()
        time = self.start_time_edit.time().toPyTime()
        trace_length = self.duration_edit.value()
        utc_string = str(date) + 'T' + str(time)
        start_time = UTCDateTime(utc_string)
        end_time = start_time + trace_length

        if self.net_edit.text() in ['*', '']:
            net = '%'
        else:
            net = self.net_edit.text()

        if self.sta_edit.text() in ['*', '']:
            sta = '%'
        else:
            sta = self.sta_edit.text()

        if self.cha_edit.text() in ['*', '']:
            cha = '%'
        else:
            cha = self.cha_edit.text()

        wfs = database.query_db(session, start_time=start_time, end_time=end_time, sta=sta, cha=cha, return_type='wfdisc_rows')
        self.parent.ipdatabase_query_results_table.setData(wfs)

    #def update_query_string(self):


    @pyqtSlot(str, str)
    def errorPopup(self, message, title="Oops..."):
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setText(message)
        msg_box.setWindowTitle(title)
        msg_box.exec_()