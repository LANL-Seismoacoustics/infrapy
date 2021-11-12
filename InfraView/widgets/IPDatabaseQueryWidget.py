from PyQt5.QtWidgets import (QDateEdit, QFormLayout, QFrame, QHBoxLayout, QLabel, QLineEdit, 
                             QMessageBox, QPushButton, QSpinBox, QTimeEdit, QVBoxLayout)
from PyQt5.QtCore import QDate, pyqtSlot

from obspy.core.utcdatetime import UTCDateTime

import pandas as pd
# from tabulate import tabulate

from infrapy.utils import database


class IPDatabaseQueryWidget(QFrame):

    def __init__(self, parent):
        super().__init__()
        self.setFrameStyle(QFrame.Box|QFrame.Plain)
        self.parent = parent  # reference to the IPBeamformingWidget to which this belongs
        self.session = None
        self.buildUI()

    def buildUI(self):

        title_label = QLabel("\tDatabase Query")
        title_label.setStyleSheet("QLabel {font-weight:bold; color: white; background-color: black}")

        self.net_edit = QLineEdit()
        self.net_edit.setMaximumWidth(150)
        self.net_edit.setToolTip('Wildcards OK \nCan be SEED network codes or data center defined codes. \nMultiple codes are comma-separated (e.g. "IU,TA").')
        self.sta_edit = QLineEdit()
        self.sta_edit.setMaximumWidth(150)
        self.sta_edit.setToolTip('Wildcards OK \nOne or more SEED station codes. \nMultiple codes are comma-separated (e.g. "ANMO,PFO")')
        self.loc_edit = QLineEdit()
        self.loc_edit.setMaximumWidth(150)
        self.loc_edit.setToolTip('One or more SEED location identifiers. \nMultiple identifiers are comma-separated (e.g. "00,01"). \nAs a special case “--“ (two dashes) will be translated to a string of two space characters to match blank location IDs.')
        self.cha_edit = QLineEdit()
        self.cha_edit.setMaximumWidth(150)
        self.cha_edit.setToolTip('Wildcards OK \nOne or more SEED channel codes. \nMultiple codes are comma-separated (e.g. "BHZ,HHZ")')

        form_layout1 = QFormLayout()
        form_layout1.addRow("Net: ", self.net_edit)
        form_layout1.addRow("Sta: ", self.sta_edit)
        form_layout1.addRow("Loc: ", self.loc_edit)
        form_layout1.addRow("Cha: ", self.cha_edit)

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
        self.duration_edit.setValue(600)

        form_layout2 = QFormLayout()
        form_layout2.addRow("Start date: ", self.start_date_edit)
        form_layout2.addRow("Start time(UTC): ", self.start_time_edit)
        form_layout2.addRow("Duration (s): ", self.duration_edit)

        horiz_layout = QHBoxLayout()
        horiz_layout.addLayout(form_layout1)
        horiz_layout.addLayout(form_layout2)
        horiz_layout.addStretch()

        self.query_button = QPushButton("Query Database")

        horiz_layout_2 = QHBoxLayout()
        horiz_layout_2.addStretch()
        horiz_layout_2.addWidget(self.query_button)
        horiz_layout_2.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addWidget(title_label)
        main_layout.addLayout(horiz_layout)
        main_layout.addLayout(horiz_layout_2)
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
            net = '*'
        else:
            net = self.net_edit.text()

        if self.sta_edit.text() in ['*', '']:
            sta = '*'
        else:
            sta = self.sta_edit.text()

        if self.loc_edit.text() in ['*', '']:
            loc = '*'
        else:
            loc = self.loc_edit.text()

        if self.cha_edit.text() in ['*', '']:
            cha = '*'
        else:
            cha = self.cha_edit.text()

        df = database.query_db(session, start_time, end_time, net="*", sta="*", loc="*", cha="*")
        print(df)
    @pyqtSlot(str, str)
    def errorPopup(self, message, title="Oops..."):
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setText(message)
        msg_box.setWindowTitle(title)
        msg_box.exec_()
