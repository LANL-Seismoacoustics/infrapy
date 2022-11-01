from PyQt5 import QtGui
from PyQt5.QtWidgets import (QDialog, QPushButton, QButtonGroup, QHBoxLayout, QVBoxLayout,
                             QFormLayout, QCheckBox, QDialogButtonBox, QGroupBox, QLabel, QLineEdit)

from PyQt5.QtCore import Qt

class IPNewDetectionsDialog(QDialog):

    detections = []
    lat = None
    lon = None
    ele = None
    method = ''
    fr = None

    current_idx = 0
    detection_count = 0

    names = []      # list of strings
    events = []     # list of strings
    notes = []      # list of strings

    def __init__(self, parent):
        super(IPNewDetectionsDialog, self).__init__(parent)

        self.buildUI()

    def exec_(self, detections, lat, lon, ele, method, fr):
        self.current_idx = 0

        self.detections = detections
        self.lat = lat
        self.lon = lon
        self.ele = ele
        self.method = method
        self.fr = fr

        self.keep = [True] * len(detections)
        self.names = [''] * len(detections)
        self.notes = [''] * len(detections)
        self.events = [''] * len(detections)
        

        self.lineEdit_detection_name.setText('')
        self.lineEdit_event.setText('')
        self.lineEdit_note.setText('')

        # dont call this until you've populated the variables
        self.update_widget()

        self.lineEdit_detection_name.setFocus()

        return super().exec_()

    def buildUI(self):

        self.setWindowTitle(self.tr('InfraView: New Detection(s)'))

        self.detection_group_box = QGroupBox()

        self.data_label = QLabel('')

        detection_name_label = QLabel("Detection Name: ")
        self.lineEdit_detection_name = QLineEdit()
        self.lineEdit_detection_name.textChanged.connect(self.update_names)

        detection_note_label = QLabel("Note (optional): ")
        self.lineEdit_note = QLineEdit()
        self.lineEdit_note.textChanged.connect(self.update_notes)

        detection_event_label = QLabel("Event: (optional)")
        self.lineEdit_event = QLineEdit()
        self.lineEdit_event.textChanged.connect(self.update_events)

        formBox = QFormLayout()
        formBox.addRow(detection_name_label, self.lineEdit_detection_name)
        formBox.addRow(detection_event_label, self.lineEdit_event)
        formBox.addRow(detection_note_label, self.lineEdit_note)

        self.keep_checkbox = QCheckBox('Keep')
        self.keep_checkbox.clicked.connect(self.update_keep)
        self.keep_checkbox.setChecked(True)
        self.discard_checkbox = QCheckBox('Discard')
        self.discard_checkbox.clicked.connect(self.update_keep)
        self.check_button_group = QButtonGroup()
        self.check_button_group.setExclusive(True)
        self.check_button_group.addButton(self.keep_checkbox)
        self.check_button_group.addButton(self.discard_checkbox)
        checkbox_layout = QHBoxLayout()
        checkbox_layout.addStretch(1)
        checkbox_layout.addWidget(self.keep_checkbox)
        checkbox_layout.addWidget(self.discard_checkbox)
        checkbox_layout.addStretch(1)

        self.next_button = QPushButton("Next")
        self.next_button.clicked.connect(self.increment_current_idx)
        self.prev_button = QPushButton("Prev")
        self.prev_button.clicked.connect(self.decrement_current_idx)
        next_button_layout = QHBoxLayout()
        next_button_layout.addWidget(self.prev_button)
        next_button_layout.addStretch(1)
        next_button_layout.addWidget(self.next_button)

        v_layout = QVBoxLayout()
        v_layout.addWidget(self.data_label)
        #v_layout.addWidget(QWidget.HLine())  # This is a horizontal line
        v_layout.addLayout(formBox)
        # v_layout.addWidget(QWidget.HLine())  # This is a horizontal line
        v_layout.addLayout(checkbox_layout)
        v_layout.addStretch(0)
        v_layout.addLayout(next_button_layout)

        self.detection_group_box.setLayout(v_layout)

        # Close Button
        buttons = QDialogButtonBox(QDialogButtonBox.Ok,
                                   Qt.Horizontal,
                                   self)
        buttons.button(QDialogButtonBox.Ok).setText("Finish")
        buttons.accepted.connect(self.accept)

        main_layout = QVBoxLayout()        
        main_layout.addWidget(self.detection_group_box)
        main_layout.addWidget(buttons)

        self.setLayout(main_layout)

    def update_notes(self):
        self.notes[self.current_idx] = self.lineEdit_note.text()

    def update_events(self):
        self.events[self.current_idx] = self.lineEdit_event.text()

    def update_names(self):
        self.names[self.current_idx] = self.lineEdit_detection_name.text()

    def increment_current_idx(self):
        self.current_idx = self.current_idx + 1
        self.update_widget()

    def decrement_current_idx(self):
        self.current_idx = self.current_idx - 1
        self.update_widget()

    def update_keep(self):
        self.keep[self.current_idx] = self.keep_checkbox.isChecked()

    def update_widget(self):

        d = self.detections[self.current_idx]
        self.detection_count = len(self.detections)

        self.next_button.setEnabled(self.current_idx + 1 < self.detection_count)
        self.prev_button.setEnabled(self.current_idx > 0)

        self.detection_group_box.setTitle(str(self.current_idx + 1) + ' of ' + str(self.detection_count))
        # beamforming_new returns detections in the order det_time, det_start, det_end, back_az, trc_vel, fstat
        data_string = ("Time (UTC): {t}\n"
                       "F Statistic: {fs:.2f}\n"
                       "Trace Velocity: {tv:.2f} m/s\n"
                       "Back Azimuth: {baz:+.2f}\n"
                       "Latitude: {la:+.2f}\n"
                       "Longitude: {lo:+.2f}\n"
                       "Elevation (m): {el:.2f} m\n"
                       "Method: {meth}\n"
                       "Detection Start/Stop: ({start:.2f}, {stop:.2f})\n"
                       "Frequency Range: ({f1:.2f}, {f2:.2f})\n ").format(t=d[0],
                                                                       fs=d[5],
                                                                       tv=d[4],
                                                                       baz=d[3],
                                                                       la=self.lat,
                                                                       lo=self.lon,
                                                                       el=self.ele,
                                                                       meth=self.method,
                                                                       start=d[1],
                                                                       stop=d[2],
                                                                       f1=self.fr[0],
                                                                       f2=self.fr[1])
        self.data_label.setText(data_string)
        self.lineEdit_detection_name.setText(self.names[self.current_idx])
        self.lineEdit_event.setText(self.events[self.current_idx])
        self.lineEdit_note.setText(self.notes[self.current_idx])
        self.keep_checkbox.setChecked(self.keep[self.current_idx])

    def get_detections(self):
        # This is called after the dialog is closed.  Will only return the detections that
        # we're marked to keep
        detections_to_keep = [i for idx, i in enumerate(self.detections) if self.keep[idx]]

        return detections_to_keep

    def get_names(self):
        names_to_keep = [i for idx, i in enumerate(self.names) if self.keep[idx]]
        return names_to_keep

    def get_notes(self):
        notes_to_keep = [i for idx, i in enumerate(self.notes) if self.keep[idx]]
        return notes_to_keep

    def get_events(self):
        events_to_keep = [i for idx, i in enumerate(self.events) if self.keep[idx]]
        return events_to_keep
    


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

        self.setWindowTitle(self.tr('InfraView: New Detection'))

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
