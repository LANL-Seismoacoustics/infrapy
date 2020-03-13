from PyQt5 import QtGui
from PyQt5.QtWidgets import (QDialog, QWidget, QVBoxLayout,
                             QFormLayout, QComboBox, QDialogButtonBox,
                             QLabel, QLineEdit)

from PyQt5.QtCore import Qt, QRegExp
from PyQt5.QtGui import QRegExpValidator

import string


class IPNewDetectionDialog(QDialog):

    def __init__(self, parent):
        super(IPNewDetectionDialog, self).__init__(parent)

        self.__buildUI__()

    def exec_(self, pick_time, f_stat, trace_vel, back_az, lat, lon, ele, method, fr):

        self.lineEdit_pick_name.setText('')
        self.lineEdit_event.setText('')
        self.lineEdit_note.setText('')

        data_string = ("Time (UTC): {t}\n"
                       "F Statistic: {fs:.2f}\n"
                       "Trace Velocity: {tv:.2f} m/s\n"
                       "Back Azimuth: {baz:+.2f}\n"
                       "Latitude: {la:+.2f}\n"
                       "Longitude: {lo:+.2f}\n"
                       "Elevation (m): {el:.2f} m\n"
                       "Method: {meth}\n"
                       "Frequency Range: ({f1:.2f}, {f2:.2f})").format(t=pick_time,
                                                                       fs=f_stat,
                                                                       tv=trace_vel,
                                                                       baz=back_az,
                                                                       la=lat,
                                                                       lo=lon,
                                                                       el=ele,
                                                                       meth=method,
                                                                       f1=fr[0],
                                                                       f2=fr[1])

        self.data_label.setText(data_string)

        self.lineEdit_pick_name.setFocus()

        return super().exec_()

    def __buildUI__(self):

        self.setWindowTitle(self.tr('InfraView - New Detection'))

        self.data_label = QLabel('')

        pick_name_label = QLabel("Detection Name: ")
        self.lineEdit_pick_name = QLineEdit()

        pick_note_label = QLabel("Note (optional): ")
        self.lineEdit_note = QLineEdit()

        pick_event_label = QLabel("Event: (optional)")
        self.lineEdit_event = QLineEdit()

        formBox = QFormLayout()
        formBox.addRow(pick_name_label, self.lineEdit_pick_name)
        formBox.addRow(pick_event_label, self.lineEdit_event)
        formBox.addRow(pick_note_label, self.lineEdit_note)

        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
                                   Qt.Horizontal,
                                   self)
        buttons.button(QDialogButtonBox.Ok).setText("Add Detection")
        buttons.button(QDialogButtonBox.Cancel).setText("Discard")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        v_layout = QVBoxLayout()
        v_layout.addWidget(self.data_label)
        # v_layout.addWidget(QWidget.HLine())  # This is a horizontal line
        v_layout.addLayout(formBox)
        # v_layout.addWidget(QWidget.HLine())  # This is a horizontal line
        v_layout.addWidget(buttons)

        self.setLayout(v_layout)

    def reset(self):
        pass

    def getName(self):
        return self.lineEdit_pick_name.text()

    def getNote(self):
        return self.lineEdit_note.text()

    def getEvent(self):
        return self.lineEdit_event.text()
