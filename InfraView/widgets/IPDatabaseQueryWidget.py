from PyQt5.QtWidgets import (QDateEdit, QDateTimeEdit, QDoubleSpinBox, QFormLayout, QFrame, QHBoxLayout, QLabel, QLineEdit, 
                             QPushButton, QSpinBox, QTimeEdit, QVBoxLayout,
                             QPlainTextEdit, QSizePolicy)
from PyQt5.QtCore import QDate, QTime, pyqtSlot, Qt

from obspy.core.utcdatetime import UTCDateTime

from infrapy.utils import database

from InfraView.widgets import IPUtils


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
        validator = IPUtils.CapsValidator()

        title_label = QLabel("\tEvent Query")
        title_label.setStyleSheet("QLabel {font-weight:bold; color: white; background-color: black}")

        self.evid_edit = QLineEdit()
        self.evid_edit.setMaximumWidth(150)
        self.evid_edit.setToolTip('')

        self.startDateTime_edit = QDateTimeEdit()
        self.startDateTime_edit.setDisplayFormat('yyyy-MM-ddTHH:mm:ss.zzz')
        #self.startDateTime_edit.setDateTime(self.startDateTime_edit.minimumDateTime())

        self.endDateTime_edit = QDateTimeEdit()
        self.endDateTime_edit.setDisplayFormat('yyyy-MM-ddTHH:mm:ss.zzz')
        #self.endDateTime_edit.setDateTime(self.endDateTime_edit.minimumDateTime())

        self.lat_edit = QDoubleSpinBox()
        self.lat_edit.setRange(-90.1,90.0)  # the -90.1 is used as the "unset" value 
        self.lat_edit.setDecimals(2)
        self.lat_edit.setMaximumWidth(80)
        self.lat_edit.setSpecialValueText('deg')
        self.lat_edit.setValue(self.lat_edit.minimum())

        self.lon_edit = QDoubleSpinBox()
        self.lon_edit.setRange(-180.1, 180.0)
        self.lon_edit.setDecimals(2)
        self.lon_edit.setMaximumWidth(80)
        self.lon_edit.setSpecialValueText('deg')
        self.lon_edit.setValue(self.lon_edit.minimum())

        self.radius_edit = QSpinBox()
        self.radius_edit.setMinimum(0)
        self.radius_edit.setMaximumWidth(80)
        self.radius_edit.setValue(0)
        self.radius_edit.setSpecialValueText('km')

        #latlon_layout = QHBoxLayout()
        #latlon_layout.addWidget

        form_layout1 = QFormLayout()
        form_layout1.addRow("EVID: ", self.evid_edit)
        #form_layout1.addRow("Start Day/Time: ", self.startDateTime_edit)
        #form_layout1.addRow("End Day/Time: ", self.endDateTime_edit)
        #form_layout1.addRow("Lat: ", self.lat_edit)
        #form_layout1.addRow("Lon: ", self.lon_edit)
        #form_layout1.addRow("Radius: ", self.radius_edit)

        self.query_textEdit = QPlainTextEdit()
        self.query_textEdit.setMaximumHeight(120)
        self.query_textEdit.setReadOnly(True)
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

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.clear_button.clicked.connect(self.clear_form)
        self.evid_edit.textEdited.connect(self.update_query_text)
        self.query_button.clicked.connect(self.query_database)

    def clear_form(self):
        self.evid_edit.setText("")
        self.query_textEdit.setPlainText("")

    def update_query_text(self):
        q = self.query_database(asquery=True)
        self.query_textEdit.setPlainText(str(q))

    def get_current_session(self):
        return self.parent.ipdatabase_connect_widget.session

    def get_tables(self):
        table_dictionary= self.parent.ipdatabase_connect_widget.table_dialog.get_tables_from_text()
        return table_dictionary

    def get_schema(self):
        return self.parent.ipdatabase_connect_widget.schema_type_combo.currentText()

    def query_database(self, asquery=False):
        session = self.get_current_session()
        if session is None:
            IPUtils.errorPopup("No current active session")
            return

        # first thing to do is assemble the info for the query
        tables = self.get_tables()
        db_tables = database.make_tables_from_dict(tables=tables, schema=self.get_schema())

        if asquery:
            return database.eventID_query(session, self.evid_edit.text(), db_tables, asquery=True)
        else:
            origins = database.eventID_query(session, self.evid_edit.text(), db_tables, asquery=False)
            self.parent.ipevent_query_results_table.setData(origins) 


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

        validator = IPUtils.CapsValidator()

        title_label = QLabel("\tWaveform Query")
        title_label.setStyleSheet("QLabel {font-weight:bold; color: white; background-color: black}")

        self.sta_edit = QLineEdit()
        self.sta_edit.setMaximumWidth(150)
        self.sta_edit.setToolTip('Wildcards OK \nOne or more SEED station codes. \nMultiple codes are comma-separated (e.g. "ANMO,PFO")')
        self.sta_edit.setValidator(validator)

        self.cha_edit = QLineEdit()
        self.cha_edit.setText("BDF")
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
        form_layout1.addRow("Station: ", self.sta_edit)
        form_layout1.addRow("Channel: ", self.cha_edit)
        form_layout1.addRow("Start date: ", self.start_date_edit)
        form_layout1.addRow("Start time(UTC): ", self.start_time_edit)
        form_layout1.addRow("Duration (s): ", self.duration_edit)

        self.query_textEdit = QPlainTextEdit()
        self.query_textEdit.setMaximumHeight(120)
        self.query_textEdit.setReadOnly(True)
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
        self.clear_button.clicked.connect(self.clear_form)
        self.query_button.clicked.connect(self.query_database)
        self.sta_edit.textEdited.connect(self.update_query_string)
        self.cha_edit.textEdited.connect(self.update_query_string)
        self.start_date_edit.dateChanged.connect(self.update_query_string)
        self.start_time_edit.timeChanged.connect(self.update_query_string)
        self.duration_edit.valueChanged.connect(self.update_query_string)

    def clear_form(self):
        self.sta_edit.setText("")
        self.cha_edit.setText("BDF")
        self.start_date_edit.setDate(self.start_date_edit.minimumDate())
        self.start_time_edit.setTime(QTime(00,00,00))
        self.duration_edit.setValue(600)
        self.query_textEdit.setPlainText("")

    def update_query_string(self):
        session = self.get_current_session()
        if session is None:
            IPUtils.errorPopup("No current active session")
            return

        # get the start date and time into a UTCDateTime
        date = self.start_date_edit.date().toPyDate()
        time = self.start_time_edit.time().toPyTime()
        trace_length = self.duration_edit.value()
        utc_string = str(date) + 'T' + str(time)
        start_time = UTCDateTime(utc_string)
        end_time = start_time + trace_length

        if self.sta_edit.text() in ['*', '']:
            sta = '%'
        else:
            sta = self.sta_edit.text()

        if self.cha_edit.text() in ['*', '']:
            cha = '%'
        else:
            cha = self.cha_edit.text()

        tables = self.parent.ipdatabase_connect_widget.table_dialog.get_tables_from_text()
        try:
            new_query = database.query_db(session, tables, start_time=start_time, end_time=end_time, sta=sta, cha=cha, return_type='wfdisc_rows', asquery=True)
            self.query_textEdit.setPlainText(str(new_query))
        except KeyError as e:
            IPUtils.errorPopup(str(e) + " is not defined.  Have you defined all your tables?")

    def get_current_session(self):
        return self.parent.ipdatabase_connect_widget.session

    def get_startstop_times(self):
        # returns UTCDateTime objects of the start and stop times
        starttime_str = self.start_date_edit.date().toString(Qt.ISODate) + "T" + self.start_time_edit.time().toString(Qt.ISODate)
        starttime = UTCDateTime(starttime_str)
        stoptime = UTCDateTime(starttime) + self.duration_edit.value()
        return starttime, stoptime

    def query_database(self):

        session = self.get_current_session()
        if session is None:
            IPUtils.errorPopup("No current active session")
            return

        # get the start date and time into a UTCDateTime
        date = self.start_date_edit.date().toPyDate()
        time = self.start_time_edit.time().toPyTime()
        trace_length = self.duration_edit.value()
        utc_string = str(date) + 'T' + str(time)
        start_time = UTCDateTime(utc_string)
        end_time = start_time + trace_length

        if '*' in self.sta_edit.text():
            sta = self.sta_edit.text().replace('*', '%')
        elif self.sta_edit.text().strip() == '':
            sta = '%'
        else:
            sta = self.sta_edit.text()

        if '*' in self.cha_edit.text():
            cha = self.cha_edit.text().replace('*', '%')
        elif self.cha_edit.text().strip() == '':
            cha = '%'
        else:
            cha = self.cha_edit.text()

        tables = self.parent.ipdatabase_connect_widget.table_dialog.get_tables_from_text()
        
        wfs = database.query_db(session, tables, start_time=start_time, end_time=end_time, sta=sta, cha=cha, return_type='wfdisc_rows')
        if len(wfs) > 0:
            self.parent.ipdatabase_query_results_table.setData(wfs)
        else:
            IPUtils.errorPopup("No results found")
