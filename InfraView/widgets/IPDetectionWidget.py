import os
import json

from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import (QFileDialog, QWidget, QPushButton,
                             QLabel, QGridLayout, QHBoxLayout,
                             QVBoxLayout, QLayout)

from PyQt5.QtCore import pyqtSignal, pyqtSlot, QSettings

import pandas as pd

from obspy import UTCDateTime

from InfraView.widgets import IPDetectionTableView, IPPickItem, IPNewDetectionDialog, IPPickLine
from InfraView.widgets import IPUtils


class IPDetectionWidget(QWidget):

    signal_detections_changed = pyqtSignal(list)
    signal_detections_cleared = pyqtSignal()

    _savefile = ("","")        # this should hold the file with it's full path

    _detections = []

    _moving_detection = None

    def __init__(self, parent):
        super().__init__(parent)

        self.parent = parent

        self.buildUI() 

        self.show()

    def buildUI(self):

        self.buildIcons()

        self.detection_view = IPDetectionTableView.IPDetectionTableView(self)

        self.fileLabel = QLabel('current file: ')

        self.clearButton = QPushButton(' Clear')
        self.clearButton.setIcon(self.clearIcon)

        self.saveButton = QPushButton(' Save')
        self.saveButton.setIcon(self.saveIcon)

        self.saveAsButton = QPushButton(' Save As...')
        self.saveAsButton.setIcon(self.saveAsIcon)

        self.loadButton = QPushButton(' Load...')
        self.loadButton.setIcon(self.openIcon)

        savebutton_group = QWidget()
        savebutton_layout = QGridLayout()

        savebutton_group.setLayout(savebutton_layout)
        savebutton_layout.addWidget(self.loadButton, 1, 0)
        savebutton_layout.addWidget(self.clearButton, 1, 1)
        savebutton_layout.addWidget(self.saveButton, 0, 0)
        savebutton_layout.addWidget(self.saveAsButton, 0, 1)
        savebutton_layout.setSizeConstraint(QLayout.SetFixedSize)

        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.fileLabel)
        vertical_layout.addWidget(savebutton_group)
        vertical_layout.addStretch(1)

        horizontal_layout = QHBoxLayout()
        horizontal_layout.addWidget(self.detection_view)
        horizontal_layout.addLayout(vertical_layout)

        self.setLayout(horizontal_layout)

        # Create the newpickdialog for later use
        self.newDetectionDialog = IPNewDetectionDialog.IPNewDetectionDialog(self)
        self.new_detections_dialog = IPNewDetectionDialog.IPNewDetectionsDialog(self)

        self.connectSignalsAndSlots()

    def buildIcons(self):
        self.clearIcon = QIcon.fromTheme("edit-clear")
        self.openIcon = QIcon.fromTheme("document-open")
        self.saveIcon = QIcon.fromTheme("document-save")
        self.saveAsIcon = QIcon.fromTheme("document-save-as")

    def connectSignalsAndSlots(self):
        # buttons
        self.loadButton.clicked.connect(self.loadDetections)
        self.clearButton.clicked.connect(self.clearDetections)
        self.saveButton.clicked.connect(self.saveDetections)
        self.saveAsButton.clicked.connect(self.saveDetectionsAs)

        self.signal_detections_changed.connect(self.set_data)

        self.detection_view.signal_delete_detections.connect(self.delete_detections)

    def get_model(self):
        return self.detection_view.get_model()

    def get_data(self):
        return self._detections

    @pyqtSlot(list)
    @pyqtSlot(pd.DataFrame)
    def set_data(self, data):
        """ SLOT:

        Data comes in as either a list of IPPickItems or dataframes, but for the tableview, we need a pandas dataframe
        so first convert the list to a menu list, then to a dataframe (there must be a more elegant way to do this?)
        """

        if isinstance(data, list):
            if len(data) < 1:
                # this sets the data to an empty pandas dataframe with the correct headers
                df = pd.DataFrame(columns=self.detection_view.columns)
                self.detection_view.set_data(df)
            else:
                pi_menu = []
                for p_i in data:
                    pi_menu.append(p_i.generateDict())
                self.detection_view.set_data(pd.DataFrame(pi_menu))
        else:  # it must be a dataframe already
            self.detection_view.set_data(data)

    def new_detections(self, 
                       detections,
                       lat,
                       lon,
                       elev=None,
                       event=None,
                       element_cnt=None,
                       method=None,
                       fr=None):

        if self.new_detections_dialog.exec_(detections, lat, lon, elev, method, fr):
            detections = self.new_detections_dialog.get_detections()
            names = self.new_detections_dialog.get_names()
            events = self.new_detections_dialog.get_events()
            notes = self.new_detections_dialog.get_notes()

            for idx, detection in enumerate(detections):



                new_detection = IPPickItem.IPPickItem(names[idx])
                new_detection.set_peakF_UTCtime(detection[0])
                new_detection.set_start(detection[1])
                new_detection.set_end(detection[2])
                new_detection.set_back_azimuth(detection[3])
                new_detection.set_trace_velocity(detection[4])
                new_detection.set_peakF_value(detection[5])
                new_detection.set_lat(lat)
                new_detection.set_lon(lon)
                new_detection.set_event_id(events[idx])
                new_detection.set_note(notes[idx])
                new_detection.set_method(method)
                new_detection.set_freq_range(fr)

                # optional info
                if elev is not None:
                    new_detection.set_ele(elev)
                if element_cnt is not None:
                    new_detection.set_array_dim(element_cnt)

                # check for duplication
                duplicate = False
                for jdx, d in enumerate(self._detections):
                    if new_detection.is_equal_to(d):
                        # we already have that detection in the list, so pop up warning and discard the new one
                        duplicate = True
                        IPUtils.errorPopup("A detection is already in the list at index {}.  The new detection will be discarded".format(jdx), title="Duplicate Detection")

                if not duplicate:
                    # append the new detection
                    self._detections.append(new_detection)
            
            self.signal_detections_changed.emit(self._detections)
            return True
        else:
            return False


    def newDetection(self,
                     pick_name,
                     time,
                     F_value,
                     trace_vel,
                     back_az,
                     lat,
                     lon,
                     elev=None,
                     event=None,
                     start_end=None,
                     element_cnt=None,
                     method=None,
                     fr=None):

        
        if self.newDetectionDialog.exec_(time,
                                         F_value,
                                         trace_vel,
                                         back_az,
                                         lat,
                                         lon,
                                         elev,
                                         method,
                                         fr):

            name = self.newDetectionDialog.getName()
            event = self.newDetectionDialog.getEvent()
            note = self.newDetectionDialog.getNote()

            new_detection = IPPickItem.IPPickItem(name)

            # required info
            new_detection.set_name(name)
            new_detection.set_peakF_UTCtime(time)
            new_detection.set_peakF_value(F_value)
            new_detection.set_lat(lat)
            new_detection.set_lon(lon)

            # optional info
            if elev is not None:
                new_detection.set_ele(elev)
            if event is not None:
                new_detection.set_event_id(event)
            if start_end is not None:
                new_detection.set_start(start_end[0])
                new_detection.set_end(start_end[1])
            if trace_vel is not None:
                new_detection.set_trace_velocity(trace_vel)
            if back_az is not None:
                new_detection.set_back_azimuth(back_az)
            if element_cnt is not None:
                new_detection.set_array_dim(element_cnt)
            if method is not None:
                new_detection.set_method(method)
            if fr is not None:
                new_detection.set_freq_range(fr)
            if note is not None:
                new_detection.set_note(note)

            # append the new detection
            self._detections.append(new_detection)

            self.signal_detections_changed.emit(self._detections)

            return True

        else:
            return False

    @pyqtSlot()
    def update_plots(self):
        # TODO
        # self.parent.plotViewer.plot_detection_lines(self.get_data())
        pass

    def saveDetections(self):

        # if there is no current filename, prompt for one...
        # TODO: if there is an open project, default to that
        if self._savefile == ("",""):
            self.saveDetectionsAs()
        else:
            data_to_save = []
            for entry in self._detections:
                data_to_save.append(entry.generateDict())
            with open(self._savefile[0], 'w') as of:
                json.dump(data_to_save, of, indent=4)
                path = os.path.dirname(self._savefile[0])
                settings = QSettings('LANL', 'InfraView')
                settings.setValue("last_detectionfile_directory", path)
                fileText = 'savefile: ' + self._savefile[0]

                # this bit is to shorten long filenames for pretty display
                if len(fileText) > 40:  # the 40 here is arbitrary, maybe there's a better way to determine that value?
                    fileText = 'savefile: ...' + fileText[-28:]
                self.fileLabel.setText(fileText)

   

    def saveDetectionsAs(self):

        if len(self._detections) == 0:
            IPUtils.errorPopup('No Detections to Save')
            return

        if self.parent.getProject() is None:
            # force a new filename...
            default_data_dir = os.path.join(os.path.dirname(__file__), '../../examples/data')
            settings = QSettings('LANL', 'InfraView')
            previous_directory = settings.value("last_detectionsfile_directory", default_data_dir)
        else:
            # There is an open project, so make the default save location correspond to what the project wants
            previous_directory = str(self.parent.getProject().get_detectionsPath())

        self._savefile = QFileDialog.getSaveFileName(self, 'Save File', previous_directory)

        if self._savefile[0]:
            data_to_save = []
            for entry in self._detections:
                data_to_save.append(entry.generateDict())
            with open(self._savefile[0], 'w') as of:
                json.dump(data_to_save, of, indent=4)
                path = os.path.dirname(self._savefile[0])
                settings = QSettings('LANL', 'InfraView')
                settings.setValue("last_detectionfile_directory", path)
                fileText = 'file: ' + self._savefile[0]

                # this bit is to shorten long filenames for pretty displayß
                if len(fileText) > 40:  # the 40 here is arbitrary, maybe there's a better way to determine that value?
                    fileText = 'file: ...' + fileText[-28:]
                self.fileLabel.setText(fileText)

    def loadDetections(self):
        # open a file and load the picks
        if self.parent.getProject() is None:
            # force a new filename...
            default_data_dir = os.path.join(os.path.dirname(__file__), '../../examples/data')
            settings = QSettings('LANL', 'InfraView')
            previous_directory = settings.value("last_detectionfile_directory", default_data_dir)
        else:
            # There is an open project, so make the default save location correspond to what the project wants
            previous_directory = str(self.parent.getProject().get_detectionsPath())

        self._openfile, _ = QFileDialog.getOpenFileName(self, 'Open File', previous_directory)

        if self._openfile == '':
            return

        with open(self._openfile, 'r') as infile:
            settings = QSettings('LANL', 'InfraView')
            settings.setValue("last_detectionfile_directory", os.path.abspath(self._openfile))
            try:
                newdata = json.load(infile)
            except json.decoder.JSONDecodeError:
                IPUtils.errorPopup("Error loading file. Perhaps this isn't a json file?")
                return

            # this is a bit of a hack to make sure opened files have the correct data headers
            for entry in newdata:
                newDetection = IPPickItem.IPPickItem()
                newDetection.fillFromDict(entry)
                self._detections.append(newDetection)

        self.signal_detections_changed.emit(self._detections)

    @pyqtSlot(list)
    def delete_detections(self, detections_to_delete):
        """SLOT:

        This method will delete a set of detections based on a list of their indexes.
        """
        detections_to_delete.sort(reverse=True)
        for idx in detections_to_delete:
            del self._detections[idx]

        self.signal_detections_changed.emit(self._detections)

    @pyqtSlot(IPPickLine.IPPickLine)
    def delete_detection(self, detectionLine):
        """SLOT:

        This method will delete a single detection based on a referece to its detectionLine, this is
        generally connected to a delete_me signal from a pickline
        """
        for detection in self._detections:
            if detection.getAssociatedPickLine() is detectionLine:
                self._detections.remove(detection)
                self.signal_detections_changed.emit(self._detections)

    @pyqtSlot(IPPickLine.IPPickLine, float)
    def detectionLineMoving(self, detectionLine, pos):
        """SLOT:
        slot called when someone is dragging a detection line.  The main purpose of this is to
        1. update the detection table with the new time
        2. assign a detection to the _moving_detection variable so that detectionLineMoved will know what to update (this might be superfluous)
        """
        est = UTCDateTime(self.parent.get_earliest_start_time())
        for detection in self._detections:
            if detection.getAssociatedPickLine() is detectionLine:
                detection.set_peakF_UTCtime(est + pos)
                self.set_data(self._detections)
                self._moving_detection = detection

    @pyqtSlot(float)
    def detectionLineMoved(self, EndingPos, fstat, backaz, tracev):
        """ SLOT:

            This slot will update the table when it receives a signal that the detection
            line was moved.
        """
        newUTCtime = UTCDateTime(self.parent.get_earliest_start_time()) + EndingPos
        self._moving_detection.set_peakF_UTCtime(newUTCtime)
        self._moving_detection.set_peakF_value(fstat)
        self._moving_detection.set_back_azimuth(backaz)
        self._moving_detection.set_trace_velocity(tracev)
        self.signal_detections_changed.emit(self._detections)

    def clearDetections(self):
        """
        Clears the accumulated list of detections and updates the table
        """
        self._detections.clear()
        self.set_data(self._detections)

        self.signal_detections_changed.emit(self._detections)
        self.signal_detections_cleared.emit()

    @pyqtSlot(IPPickLine.IPPickLine)
    def createNewDetectionStartEndRegion(self, detectionline):
        """ SLOT:
            This will add start/end information for a single detection in the detectionwidget,
            This should emit a detections_changed signal to trigger redrawing of the detection lines,
            the map, and the detection table.
        """
        if detectionline is None:
            return

        for detection in self._detections:
            if detection.getAssociatedPickLine() is detectionline:
                # we found the correct pickline, lets figure out the 'default' start/end based on
                # the width of the plots
                viewRange = self.parent.fstatPlot.viewRange()
                t_span = viewRange[0][1] - viewRange[0][0]
                start_end = [-t_span / 50., t_span / 50.]   # +/- 1 percent of the visable range

                detection.set_start(start_end[0])
                detection.set_end(start_end[1])
                detectionline.addStartEndBars(start_end)

                self.signal_detections_changed.emit(self._detections)   # should trigger a redraw

                return  # and we're done here

    @pyqtSlot(IPPickLine.IPPickLine)
    def removeDetectionStartEndRegion(self, detectionline):
        """ SLOT:

            This will remove the start/end information for a single detection in the detectionwidget,
            This should emit a detections_changed signal to trigger redrawing of the detection lines,
            the map, and the detection table.
        """
        for detection in self._detections:
            if detection.getAssociatedPickLine() is detectionline:
                detection.set_start(None)
                detection.set_end(None)
                self.signal_detections_changed.emit(self._detections)
                return  # and we're done here

    @pyqtSlot(IPPickLine.IPPickLine, list)
    def updateDetectionStartEnd(self, detectionline, start_end):
        for detection in self._detections:
            if detection.getAssociatedPickLine() is detectionline:
                detection.set_start(start_end[0])
                detection.set_end(start_end[1])
                self.set_data(self._detections)
                return  # and we're done here
