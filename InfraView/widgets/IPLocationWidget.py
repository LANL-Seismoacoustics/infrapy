import sys
import matplotlib

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import (QCheckBox, QLabel, QWidget, QBoxLayout, QHBoxLayout,
                             QVBoxLayout, QDoubleSpinBox, QSpinBox,
                             QFormLayout, QFrame, QPushButton,
                             QSplitter, QTextEdit, QComboBox)

from PyQt5.QtCore import Qt, QObject, QThread, pyqtSignal, pyqtSlot, QSettings

import numpy as np

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
from scipy.spatial.distance import pdist, squareform

import cartopy.feature as cfeature

from infrapy.location import bisl
from infrapy.association import hjl

from InfraView.widgets import IPMapWidget
from InfraView.widgets import IPEventWidget
from InfraView.widgets import IPUtils

import pyqtgraph as pg
from pyqtgraph.GraphicsScene import exportDialog

# Make sure that we are using QT5
matplotlib.use('Qt5Agg')


class IPLocationWidget(QWidget):

    bisl_result = None

    detections = []
    trimmed_detections = []

    mp_pool = None  # multiprocessing pool

    signal_start_dist_calc = pyqtSignal()
    signal_start_BISL_calc = pyqtSignal()
    signal_start_cluster_calc = pyqtSignal()

    def __init__(self, parent, pool):
        super().__init__()
        self.parent = parent

        self.mp_pool = pool

        self.buildUI()

    def buildUI(self):

        # BottomTab widgets go here...
        self.consoleBox = QTextEdit()
        self.consoleBox.setReadOnly(True)

        # set up the map widget
        self.mapWidget = IPMapWidget.IPMapWidget(self)

        # set up the distance matrix viewer
        self.dm_view = IPDistanceMatrixWidget(self)

        # set up dendrogram widget
        self.dendrogram = IPDendrogramWidget(self)

        # set up the bisl settings widget
        self.bislSettings = BISLSettings(self)

        # set up showgroundtruth widget
        self.showgroundtruth = ShowGroundTruth(self)
        self.showgroundtruth.event_widget.sigEventWidgetChanged.connect(self.parent.waveformWidget.plotViewer.pl_widget.updateEventLines)
        self.showgroundtruth.event_widget.sigEventWidgetChanged.connect(self.parent.waveformWidget.plotViewer.pl_widget.plotEventLines)
        self.showgroundtruth.event_widget.sigEventWidgetChanged.connect(self.mapWidget.plot_ground_truth)
        self.showgroundtruth.event_widget.showGT_cb.stateChanged.connect(self.mapWidget.show_hide_ground_truth)


        # set up association settings widget
        self.assocSettings = AssociationSettings(self)

        # right hand widgets layout holds the settings widgets
        rh_widget = QWidget()
        rh_layout = QVBoxLayout()
        rh_layout.addWidget(self.bislSettings)
        rh_layout.addWidget(self.assocSettings)
        rh_layout.addWidget(self.showgroundtruth)
        rh_layout.addStretch()
        rh_widget.setLayout(rh_layout)

        # splitter holding the association plots
        self.assoc_splitter = QSplitter(Qt.Vertical)
        self.assoc_splitter.addWidget(self.dm_view)
        self.assoc_splitter.addWidget(self.dendrogram)
        self.assoc_splitter.setSizes([100000, 100000])

        # splitter holding the map canvas and the association plots
        self.loc_splitter = QSplitter(Qt.Horizontal)
        self.loc_splitter.addWidget(self.mapWidget)
        self.loc_splitter.addWidget(self.assoc_splitter)

        # large splitter holding the map, association plots, and the console
        self.mapSplitter = QSplitter(Qt.Vertical)
        self.mapSplitter.addWidget(self.loc_splitter)
        self.mapSplitter.addWidget(self.consoleBox)

        self.mainSplitter = QSplitter(Qt.Horizontal)
        self.mainSplitter.addWidget(self.mapSplitter)
        self.mainSplitter.addWidget(rh_widget)

        main_layout = QBoxLayout(QBoxLayout.TopToBottom)
        main_layout.addWidget(self.mainSplitter)
        self.setLayout(main_layout)

        self.connectSignalsAndSlots()

        # Create threads for the distancematrix calculation, BISL, and clustering
        self.dmThread = QThread()
        self.bislThread = QThread()
        self.clusterThread = QThread()

    def connectSignalsAndSlots(self):

        self.bislSettings.run_bisl_button.clicked.connect(self.run_bisl)
        self.bislSettings.update_dm_button.clicked.connect(self.calc_distance_matrix)
        self.bislSettings.rng_max_edit.valueChanged.connect(self.mapWidget.update_range_max)

        self.bislSettings.confidence_edit.valueChanged.connect(self.bislSettings.enable_update_dm_button)
        self.bislSettings.confidence_edit.valueChanged.connect(self.calc_conf_ellipse)

        self.assocSettings.dist_max_edit.valueChanged.connect(self.dm_adjust_max_distance)
        self.assocSettings.threshold_edit.valueChanged.connect(self.cluster_adjust_threshold)
        self.assocSettings.update_assoc_button.clicked.connect(self.calc_associations)

        self.dm_view.signal_trim_detections.connect(self.trim_detections)
        self.dendrogram.signal_new_colors.connect(self.dm_view.set_colors)

    @pyqtSlot()
    def detections_cleared(self):
        self.mapWidget.clear_plot()
        self.detections = []
        self.dm_view.clear()
        self.dendrogram.clear_plot()
        self.consoleBox.clear()
    
    def get_detections(self):
        return self.detections

    def get_trimmed_detections(self):
        return self.trimmed_detections

    @pyqtSlot(list)
    def update_detections(self, new_detections, detection_type='ip_detections', recalc_assoc=True):
        if new_detections is None:
            return  # Nothing to do

        if len(new_detections) < 1:
            self.detections_cleared()
            return

        self.detections = []
        if detection_type == "ip_detections":   # we need to covert to InfrasoundDetections
            for detection in new_detections:
                self.detections.append(detection.to_InfrasoundDetection())

        else:
            for detection in new_detections:
                self.detections.append(detection)

        self.trimmed_detections = self.detections

        self.mapWidget.update_detections()

        if recalc_assoc:
            self.calc_distance_matrix()

    @pyqtSlot(list, str)
    def trim_detections(self, indicies, linecolor='gray'):

        # the detections to show has been changed, which means we probably don't 
        # want the bisl results showing anymore.  So lets remove those first.
        self.mapWidget.remove_conf_ellipse()
        self.mapWidget.remove_bisl_result()

        self.trimmed_detections = []
        if len(self.detections) < 1:
            return  # nothing to do

        # lets pick out the detections that we want to show
        for index in indicies:
            self.trimmed_detections.append(self.detections[index])
            self.trimmed_detections[-1].index = index

        self.mapWidget.update_detections(line_color=linecolor)

    def run_bisl(self):
        if not self.dm_view.is_group_selected():
            IPUtils.errorPopup("You need to select a cluster in the Distance Matrix to run BISL on.")
            return  # nothing to do
        if self.trimmed_detections is None:
            IPUtils.errorPopup("no detections loaded. \n You need at least two detections to run BISL.")
            return  # nothing to do

        if len(self.trimmed_detections) < 2:
            IPUtils.errorPopup("not enough detections loaded. \n You need two or more detections to run BISL.")
            return  # you need at least 2 detections to calculate the dist matrix

        # if there are previous bisl results on the map, remove them now
        self.mapWidget.remove_bisl_result()
        self.mapWidget.remove_conf_ellipse()

        self.bisl_workerObject = BISLWorkerObject(self.trimmed_detections,
                                                  beam_width=self.bislSettings.bm_width_edit.value(),
                                                  rad_min=self.bislSettings.rad_min_edit.value(),
                                                  rad_max=self.bislSettings.rad_max_edit.value(),
                                                  rng_max=self.bislSettings.rng_max_edit.value(),
                                                  resol=self.bislSettings.resolution_edit.value())

        self.bisl_workerObject.moveToThread(self.bislThread)

        self.signal_start_BISL_calc.connect(self.bisl_workerObject.run)
        self.bisl_workerObject.signal_runFinished.connect(self.bisl_run_finished)

        # start the thread
        self.consoleBox.setText("...Calculating...")
        self.bislThread.start()
        self.signal_start_BISL_calc.emit()

    @pyqtSlot(dict)
    def bisl_run_finished(self, result):
        self.bisl_result = result

        self.consoleBox.setText(bisl.summarize(result, self.bislSettings.confidence_edit.value()))

        self.calc_conf_ellipse(self.bislSettings.confidence_edit.value())

    @pyqtSlot(int)
    def calc_conf_ellipse(self, confidence):

        if self.bisl_result is None:
            return  # nothing to plot

        conf_dx, conf_dy = bisl.calc_conf_ellipse([0.0, 0.0],
                                                  [self.bisl_result['EW_stdev'],
                                                  self.bisl_result['NS_stdev'],
                                                  self.bisl_result['covar']],
                                                  confidence)
        # tell the mapWidget to plot the results
        self.mapWidget.plot_bisl_result(self.bisl_result['lon_mean'],
                                        self.bisl_result['lat_mean'])

        self.mapWidget.plot_conf_ellipse(self.bisl_result['lon_mean'],
                                         self.bisl_result['lat_mean'],
                                         conf_dx,
                                         conf_dy)

    @pyqtSlot()
    def calc_distance_matrix(self):

        if len(self.detections) < 1:
            IPUtils.errorPopup("No detections loaded.\n You need two or more detections to calculate a distance matrix.")
            return  # nothing to do

        if len(self.detections) < 2:
            # IPUtils.errorPopup("not enough detections loaded. \n You need 2 or more detections to calculate a distance matrix.")
            return  # you need at least 2 detections to calculate the dist matrix

        self.dist_matrix = None

        self.dm_workerObject = DistanceMatrixWorkerObject(self.detections,
                                                          beam_width=self.bislSettings.bm_width_edit.value(),
                                                          rng_max=self.bislSettings.rng_max_edit.value(),
                                                          rad_min=self.bislSettings.rad_min_edit.value(),
                                                          rad_max=self.bislSettings.rad_max_edit.value(),
                                                          resol=self.bislSettings.resolution_edit.value(),
                                                          pool=self.mp_pool)

        self.dm_workerObject.moveToThread(self.dmThread)

        self.signal_start_dist_calc.connect(self.dm_workerObject.run)
        self.dm_workerObject.signal_runFinished.connect(self.dm_run_finished)

        # start the thread
        self.dmThread.start()
        self.signal_start_dist_calc.emit()
        self.dm_view.showCalculatingText()

        self.bislSettings.update_dm_button.setEnabled(False)

    @pyqtSlot(np.ndarray)
    def dm_run_finished(self, data):
        self.dist_matrix_orig = data    # keep this around incase someone twiddles with the max_distance setting
        self.dm_adjust_max_distance()

        self.dm_view.hideCalculatingText()

        if self.dist_matrix is not None:
            self.dm_view.set_data(self.dist_matrix)

        # Now that the distance matrix is set, calculate the association dendrogram
        self.calc_associations()

    def dm_adjust_max_distance(self):
        self.dist_matrix = self.dist_matrix_orig.copy()
        self.dist_matrix[self.dist_matrix_orig > self.assocSettings.dist_max_edit.value()] = self.assocSettings.dist_max_edit.value()
        self.assocSettings.update_assoc_button.setEnabled(True)

    def cluster_adjust_threshold(self):
        self.assocSettings.update_assoc_button.setEnabled(True)

    @pyqtSlot()
    def calc_associations(self):

        if self.dist_matrix is None:
            IPUtils.errorPopup("No distance matrix...I need a distance matrix")
            return  # Nothing to do

        self.cluster_workerObject = ClusterWorkerObject(self.dist_matrix,
                                                        threshold=self.assocSettings.threshold_edit.value())

        self.cluster_workerObject.moveToThread(self.clusterThread)

        self.signal_start_cluster_calc.connect(self.cluster_workerObject.run)
        self.cluster_workerObject.signal_runFinished.connect(self.cluster_run_finished)

        # start the thread
        self.clusterThread.start()
        self.signal_start_cluster_calc.emit()

        self.assocSettings.update_assoc_button.setEnabled(False)

    @pyqtSlot(np.ndarray, np.ndarray)
    def cluster_run_finished(self, links, labels):

        self.dendrogram.set_data(links, self.assocSettings.threshold_edit.value())

        # Sort the distance matrix using the labels
        det_cnt = len(self.dist_matrix)
        sorting = np.array([])
        for n in range(max(labels + 1)):
            sorting = np.concatenate((sorting, np.arange(det_cnt)[labels == n]))
        sorting = sorting.astype(int)

        distance_matrix_sorted = np.empty_like(self.dist_matrix)
        for n1 in range(det_cnt):
            for n2 in range(det_cnt):
                distance_matrix_sorted[n1][n2] = self.dist_matrix[sorting[n1], sorting[n2]]

        self.dm_view.set_data(distance_matrix_sorted, labels)
        self.update_detections(self.detections, detection_type='detections', recalc_assoc=False)

    def saveWindowGeometrySettings(self):
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('LocationWidget')
        settings.setValue("windowSize", self.size())
        settings.setValue("windowPos", self.pos())
        settings.setValue("mapSplitterSettings", self.mapSplitter.saveState())
        settings.setValue("mainSplitterSettings", self.mainSplitter.saveState())
        settings.setValue("assocSplitterSettings", self.assoc_splitter.saveState())
        settings.setValue("loc_splitterSettings", self.loc_splitter.saveState())
        settings.endGroup()

    def restoreWindowGeometrySettings(self):
        # Restore settings
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('LocationWidget')

        mapSplitterSettings = settings.value("mapSplitterSettings")
        if mapSplitterSettings:
            self.mapSplitter.restoreState(mapSplitterSettings)

        mainSplitterSettings = settings.value("mainSplitterSettings")
        if mainSplitterSettings:
            self.mainSplitter.restoreState(mainSplitterSettings)

        assocSplitterSettings = settings.value("assocSplitterSettings")
        if assocSplitterSettings:
            self.assoc_splitter.restoreState(assocSplitterSettings)

        locSplitterSettings = settings.value("loc_splitterSettings")
        if locSplitterSettings:
            self.loc_splitter.restoreState(locSplitterSettings)

        settings.endGroup()


class BISLSettings(QFrame):

    earth_radius = 6378.1   # km

    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.buildUI()

    def buildUI(self):

        self.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)

        title_label = QLabel('BISL')
        title_label.setToolTip('Bayesian Infrasound Source Localization')
        title_label.setStyleSheet("font-weight: bold;")
        title_label.setAlignment(Qt.AlignCenter)

        self.bm_width_edit = QDoubleSpinBox()
        self.bm_width_edit.setMinimum(2.5)
        self.bm_width_edit.setMaximum(45.0)
        self.bm_width_edit.setValue(10)
        self.bm_width_edit.setSuffix(' deg')
        self.bm_width_edit.valueChanged.connect(self.enable_update_dm_button)

        self.rad_min_edit = QDoubleSpinBox()
        self.rad_min_edit.setMinimum(50)
        self.rad_min_edit.setMaximum(np.pi * self.earth_radius)
        self.rad_min_edit.setValue(100.0)
        self.rad_min_edit.setSuffix(' km')
        self.rad_min_edit.valueChanged.connect(self.enable_update_dm_button)

        self.rad_max_edit = QDoubleSpinBox()
        self.rad_max_edit.setMinimum(50)
        self.rad_max_edit.setMaximum(np.pi * self.earth_radius)
        self.rad_max_edit.setValue(1000.0)
        self.rad_max_edit.setSuffix(' km')
        self.rad_max_edit.valueChanged.connect(self.enable_update_dm_button)

        self.rng_max_edit = QDoubleSpinBox()
        self.rng_max_edit.setMinimum(10)
        self.rng_max_edit.setMaximum(np.pi * self.earth_radius)
        self.rng_max_edit.setValue(3000.0)
        self.rng_max_edit.setSuffix(' km')
        self.rng_max_edit.valueChanged.connect(self.enable_update_dm_button)

        self.resolution_edit = QSpinBox()
        self.resolution_edit.setMinimum(10)
        self.resolution_edit.setMaximum(10000)
        self.resolution_edit.setValue(180)
        self.resolution_edit.valueChanged.connect(self.enable_update_dm_button)

        self.confidence_edit = QSpinBox()
        self.confidence_edit.setMinimum(1)
        self.confidence_edit.setMaximum(99)
        self.confidence_edit.setValue(95)
        self.confidence_edit.setSuffix(' %')

        layout = QFormLayout()
        layout.addRow(self.tr('Beam Width: '), self.bm_width_edit)
        layout.addRow(self.tr('Radius Min.: '), self.rad_min_edit)
        layout.addRow(self.tr('Radius Max.: '), self.rad_max_edit)
        layout.addRow(self.tr('Range Max.: '), self.rng_max_edit)
        layout.addRow(self.tr('Resolution'), self.resolution_edit)
        layout.addRow(self.tr('Confidence'), self.confidence_edit)

        self.run_bisl_button = QPushButton('Run BISL')
        self.run_bisl_button.setMaximumWidth(110)

        self.update_dm_button = QPushButton('Update Dist. Matrix')

        mainlayout = QVBoxLayout()
        mainlayout.addWidget(title_label)
        mainlayout.addLayout(layout)

        buttonLayout = QHBoxLayout()
        buttonLayout.addWidget(self.run_bisl_button)
        buttonLayout.addWidget(self.update_dm_button)

        mainlayout.addLayout(buttonLayout)

        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        self.setLayout(mainlayout)

    @pyqtSlot(float)
    @pyqtSlot(int)
    def enable_update_dm_button(self, _):
        self.update_dm_button.setEnabled(True)


class ShowGroundTruth(QFrame):

    sig_event_changed = pyqtSignal(float, float)

    def __init__(self, parent):
        super().__init__()

        self.parent = parent
        self.buildUI()

    def buildUI(self):
        self.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        
        title_label = QLabel('Event/Ground Truth')
        title_label.setStyleSheet("font-weight: bold;")
        title_label.setAlignment(Qt.AlignCenter)

        self.event_widget = IPEventWidget.IPEventWidget(self)

        layout = QVBoxLayout()
        layout.addWidget(title_label)
        layout.addWidget(self.event_widget)

        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        layout.setAlignment(Qt.AlignCenter)
        self.setLayout(layout)
    
    @pyqtSlot(dict)
    def eventChanged(self, event_dict):
        if event_dict['Latitude']:
            self.lat_label.setText("{:.5f}".format(event_dict['Latitude']))
        else:
            self.lat_label.setText("0.0")
        if event_dict['Longitude']:
            self.lon_label.setText("{:.5f}".format(event_dict['Longitude']))
        else:
            self.lon_label.setText("0.0")
        
        self.sig_event_changed.emit(float(self.lon_label.text()), float(self.lat_label.text()))

    def show_gt(self):
        return self.showGT_cb.isChecked()


class IPDistanceMatrixWidget(QWidget):

    N = 5
    calc_text = None

    labels = None
    sorted_labels = None
    current_group = None

    greenPen = pg.mkPen(color='g')
    whitePen = pg.mkPen(color='w')
    bluePen = pg.mkPen(color='b')
    redPen = pg.mkPen(color='r')

    color_palette_pens = []
    color_palette_str = ['b', 'r', 'g', 'c', 'm', 'y']

    # used to keep track of wether a group has been clicked on.
    group_selected = False

    signal_trim_detections = pyqtSignal(list, str)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.buildUI()

        self.s1 = IPScatterPlotItem(x=[], y=[], symbol='s', pxMode=False)
        # self.s1.signal_point_hovered.connect(self.handle_point_hovered)
        # self.s1.signal_hover_leave.connect(self.handle_hover_leave)
        self.s1.sigClicked.connect(self.handle_mouse_click)

        for c in self.color_palette_str:
            self.color_palette_pens.append(pg.mkPen(color=c))

    def buildUI(self):
        self.dm_plotitem = IPDistanceMatrixPlot()

        gl_layout = pg.GraphicsLayoutWidget()
        gl_layout.addItem(self.dm_plotitem)

        instruct_label = QLabel("Click on a cluster to choose detections to run BISL on.")

        layout = QVBoxLayout()
        layout.addWidget(gl_layout)
        layout.addWidget(instruct_label)
        self.setLayout(layout)

    def showCalculatingText(self):
        self.calc_text = pg.TextItem('...Calculating...', color=(20, 20, 20), fill=(255, 255, 255), anchor=(0.5, 0.5), border={'color': 'k', 'width': 1})
        self.dm_plotitem.addItem(self.calc_text)
        self.calc_text.setPos(self.N / 2., self.N / 2.)

    def hideCalculatingText(self):
        self.dm_plotitem.removeItem(self.calc_text)

    def set_data(self, dist_data, labels=None):
        self.s1.clear()
        self.dm_plotitem.clear()

        self.current_group = None   # on changing the distance matrix, this fixes bug where clicking on a group wouldn't do anything if it had been previously highlighted

        self.labels = labels

        self.N = dist_data.shape[0]
        squares = []
        max_dist = np.amax(dist_data, axis=(0, 1))
        self.dm_plotitem.setXRange(0, self.N, padding=0)
        self.dm_plotitem.setYRange(0, self.N, padding=0)

        # pos = [0.0, 0.5*max_dist, 0.25*max_dist, 0.75*max_dist, max_dist]
        # colors = np.array([[50,50,50,255], [255,50,50,255], [50,255,50,255], [50,50,255,255], (255,255,255,255)], dtype=np.ubyte)
        # cmap = pg.ColorMap(pos, colors)
        # map_colors = cmap.map(dist_data)

        # testcm = cm.get_cmap("nipy_spectral")

        for i in range(self.N):
            for j in range(self.N):
                if max_dist == 0:
                    # this means the data has been cleared
                    color = (255, 255, 255)
                else:
                    color = (255.0 * dist_data[i][j] / max_dist, 255.0 * dist_data[i][j] / max_dist, 255.0 * dist_data[i][j] / max_dist)

                squares.append({'pos': (i, j), 'pen': {'color': 'w', 'width': 1}, 'brush': color, 'data': dist_data[i][j]})

        self.s1.addPoints(squares)
        self.s1.setSize(1)
        self.dm_plotitem.addItem(self.s1)

        # draw labels
        # first clear out all the previous labels
        for item in reversed(self.dm_plotitem.items):
            if type(item) is pg.TextItem:
                self.dm_plotitem.removeItem(item)
                del item

        # since cluster() doesn't return the sorted labels, we need to do it here
        # Sort the distance matrix using the labels
        if self.labels is not None:
            self.sorted_labels = np.array([])
            for n in range(max(self.labels + 1)):
                self.sorted_labels = np.concatenate((self.sorted_labels, np.arange(self.N)[labels == n]))
            self.sorted_labels = self.sorted_labels.astype(int)
        else:
            self.sorted_labels = None

        # now add the sorted labels to the plot
        for i in range(self.N):
            # x-axis
            if self.sorted_labels is not None:
                tx = pg.TextItem(str(self.sorted_labels[i]), anchor=(0.5, 0), color=(0, 0, 0))
                ty = pg.TextItem(str(self.sorted_labels[i]), anchor=(0, 0.5), color=(0, 0, 0))
            else:
                tx = pg.TextItem(str(i), anchor=(0.5, 0))
                ty = pg.TextItem(str(i), anchor=(0, 0.5))
            tx.setPos(i, -0.5)
            ty.setPos(-1, i)

            self.dm_plotitem.addItem(tx)
            self.dm_plotitem.addItem(ty)

        self.xlabel = pg.TextItem('Detection Number', anchor=(0.5, 0), color=(0, 0, 0))
        self.xlabel.setPos((self.N - 1) / 2., -1.5)

        self.ylabel = pg.TextItem('Detection Number', anchor=(0.5, 0), color=(0, 0, 0), angle=90)
        self.ylabel.setPos(-2, (self.N - 1) / 2.)

        self.dm_plotitem.addItem(self.xlabel)
        self.dm_plotitem.addItem(self.ylabel)

        self.dm_plotitem.enableAutoRange()

    def clear(self):
        self.s1.clear()
        self.N = 5
        initial_data = np.zeros((self.N, self.N))
        self.set_data(initial_data)
        self.sorted_labels = None
        self.labels = None
        self.xlabel = None
        self.ylabel = None

    @pyqtSlot(pg.SpotItem)
    def handle_point_hovered(self, point):

        if self.sorted_labels is None or self.labels is None:
            return

        # for convenience
        pos_x = int(point.pos().x())
        pos_y = int(point.pos().y())

        if self.labels[self.sorted_labels[pos_x]] == self.labels[self.sorted_labels[pos_y]]:
            # the point is in a group, so we want to highlight the group
            group_num = self.labels[self.sorted_labels[pos_x]]
            # find indicies in labels array that have that grouping
            indicies = [i for i, value in enumerate(self.labels) if value == group_num]

            for pnt in self.s1.points():
                pnt.setPen(self.whitePen)

            for index_i in indicies:
                for index_j in indicies:
                    a = np.where(self.sorted_labels == index_i)[0]
                    b = np.where(self.sorted_labels == index_j)[0]
                    ps = self.s1.pointsAt(pg.Point(a, b))
                    for p in ps:
                        p.setPen(self.bluePen)

        else:
            for pnt in self.s1.points():
                pnt.setPen(self.whitePen)

    @pyqtSlot()
    def handle_hover_leave(self):
        for pnt in self.s1.points():
            pnt.setPen(self.whitePen)

    def is_group_selected(self):
        if self.current_group is None:
            return False
        else:
            return True

    @pyqtSlot(object, object)
    def handle_mouse_click(self, scatterPlot, points):

        # mouse
        for pnt in points:
            pos_x = int(pnt.pos().x())
            pos_y = int(pnt.pos().y())

            if self.labels[self.sorted_labels[pos_x]] == self.labels[self.sorted_labels[pos_y]]:
                # the point is in a group, so we want to highlight that group with a colored outline
                group_num = self.labels[self.sorted_labels[pos_x]]

                # now check to see if the clicked group is a new one, or the one currently clicked
                if group_num != self.current_group:
                    # we have a new group, so first reset all points
                    for pnt in self.s1.points():
                        pnt.setPen(self.whitePen)

                    # update current group
                    self.current_group = group_num

                    # find indicies of detections not in the group
                    indicies = [i for i, value in enumerate(self.labels) if value == group_num]

                    if group_num < len(self.color_palette_str):
                        self.signal_trim_detections.emit(indicies, self.color_palette_str[group_num])
                    else:
                        self.signal_trim_detections.emit(indicies, 'gray')

                    for index_i in indicies:
                        for index_j in indicies:
                            a = np.where(self.sorted_labels == index_i)[0]
                            b = np.where(self.sorted_labels == index_j)[0]
                            ps = self.s1.pointsAt(pg.Point(a, b))
                            for p in ps:
                                if group_num < len(self.color_palette_pens):
                                    p.setPen(self.color_palette_pens[group_num])
                                else:
                                    p.setPen(self.whitePen)
            else:
                # the clicked point is not in a group, so don't highlight anything
                self.current_group = None

                for pnt in self.s1.points():
                    pnt.setPen(self.whitePen)

                self.signal_trim_detections.emit(self.sorted_labels.tolist(), 'gray')

    @pyqtSlot(list)
    def set_colors(self, new_colors):
        self.color_palette_str = new_colors
        self.color_palette_pens.clear()
        for c in self.color_palette_str:
            self.color_palette_pens.append(pg.mkPen(color=c))


class IPScatterPlotItem(pg.ScatterPlotItem):

    signal_point_hovered = pyqtSignal(pg.SpotItem)
    signal_hover_leave = pyqtSignal()

    last_val = -1
    last_point = None

    def __init__(self, *args, **kargs):
        super().__init__(*args, **kargs)

        self.setAcceptHoverEvents(True)

    def hoverMoveEvent(self, evt):
        pts = self.pointsAt(evt.pos())
        if len(pts) > 0:
            for point in pts:
                if point is not self.last_point:
                    self.last_point = point
                    self.signal_point_hovered.emit(point)

    def hoverLeaveEvent(self, evt):
        self.signal_hover_leave.emit()


class IPDistanceMatrixPlot(pg.PlotItem):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setAspectLocked(lock=True, ratio=1)
        self.hideAxis('bottom')
        self.hideAxis('left')
        self.hideAxis('right')
        self.hideAxis('top')
        self.setTitle('Distance Matrix')

    def mouseClickEvent(self, evt):
        if evt.button() == Qt.RightButton:
            self.export_dialog = exportDialog.ExportDialog(self.scene())
            self.export_dialog.show()
            evt.accept()


class DistanceMatrixWorkerObject(QObject):

    signal_runFinished = pyqtSignal(np.ndarray)

    def __init__(self, detections,
                 beam_width=10,
                 rng_max=np.pi / 2.0 * 6370.0,
                 rad_min=100.,
                 rad_max=1000.,
                 resol=180,
                 pool=None):

        super().__init__()
        self.detections = detections
        self.beam_width = beam_width
        self.rng_max = rng_max
        self.rad_min = rad_min
        self.rad_max = rad_max
        self.resol = resol
        self.pool = pool

        self.thread_stopped = True

    @pyqtSlot()
    def run(self):

        if len(self.detections) == 0:
            return  # nothing to do

        self.thread_stopped = False

        try:
            self.dist_matrix = hjl.build_distance_matrix(self.detections,
                                                         bm_width=self.beam_width,
                                                         rng_max=self.rng_max,
                                                         rad_min=self.rad_min,
                                                         rad_max=self.rad_max,
                                                         resol=self.resol,
                                                         pool=self.pool)
        except Exception:
            IPUtils.errorPopup("Error while calculating the distance matrix: {}".format(sys.exc_info()[0]))
            self.thread_stopped = True
            return

        self.signal_runFinished.emit(self.dist_matrix)

    @pyqtSlot()
    def stop(self):
        self.thread_stopped = True


class BISLWorkerObject(QObject):

    signal_runFinished = pyqtSignal(dict)

    def __init__(self, detections,
                 beam_width=10,
                 rad_min=100.,
                 rad_max=1000.,
                 rng_max=np.pi / 2.0 * 6370.0,
                 resol=180):

        super().__init__()
        self.detections = detections
        self.beam_width = beam_width
        self.rng_max = rng_max
        self.rad_min = rad_min
        self.rad_max = rad_max
        self.resol = resol

        self.thread_stopped = True

    @pyqtSlot()
    def run(self):
        if len(self.detections) == 0:
            return  # nothing to do

        self.thread_stopped = False

        # run bisl
        try:
            self.bisl_result = bisl.run(self.detections,
                                         bm_width=self.beam_width,
                                         rad_min=self.rad_min,
                                         rad_max=self.rad_max,
                                         rng_max=self.rng_max,
                                         resol=self.resol)
        except Exception:
            IPUtils.errorPopup("Error while running BISL: {}".format(sys.exc_info()[0]))
            self.thread_stopped = True
            return

        self.signal_runFinished.emit(self.bisl_result)

    @pyqtSlot()
    def stop(self):
        self.threadStopped = True


class ClusterWorkerObject(QObject):

    signal_runFinished = pyqtSignal(np.ndarray, np.ndarray)

    def __init__(self, dm,
                 threshold,
                 linkage_method='weighted'):

        super().__init__()
        self.dist_matrix = dm
        self.threshold = threshold
        self.linkage_method = linkage_method

        self.thread_stopped = True

    @pyqtSlot()
    def run(self):

        det_cnt = len(self.dist_matrix)
        if det_cnt == 0:
            return  # nothing to do

        self.thread_stopped = False

        # run clustering
        try:
            links = linkage(squareform(self.dist_matrix), self.linkage_method)
        except Exception:
            IPUtils.errorPopup("Error while calculating the linkage: {}".format(sys.exc_info()))
            self.thread_stopped = True
            return

        try:
            labels = fcluster(links, self.threshold, criterion='distance') - 1
        except Exception:
            IPUtils.errorPopup("Error while calculating the labels: {}".format(sys.exc_info()))
            self.thread_stopped = True
            return

        self.signal_runFinished.emit(links, labels)

    @pyqtSlot()
    def stop(self):
        self.thread_stopped = True


class IPDendrogramWidget(QWidget):

    signal_new_colors = pyqtSignal(list)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.fig = Figure()
        self.axes = self.fig.add_subplot(111)
        self.axes.set_title('Associations')
        self.axes.title.set_size(10)

        self.axes.tick_params(axis='both', labelsize=8)

        self.axes.set_xlabel('Detection Number')
        self.axes.set_ylabel('Distance')
        self.axes.xaxis.label.set_size(8)
        self.axes.yaxis.label.set_size(8)

        self.canvas = FigureCanvas(self.fig)

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)

        self.setLayout(layout)

    def set_data(self, links, threshold):
        self.axes.clear()
        self.axes.set_title('Associations', fontsize=10)

        # The link color palette needs to match the color palette in the distance matrix widget!!!
        set_link_color_palette(['#006ba6', '#ce1126', '#428a17', '#ffcc33', '#008080', 'm', '#ff4570', '#ff9000', 'b', 'g', 'c'])

        den = dendrogram(links, ax=self.axes, leaf_rotation=0., leaf_font_size=8, color_threshold=threshold, above_threshold_color='0.5')

        den_colors = []
        for c in den['color_list']:
            if not self.is_number(c):
                den_colors.append(c)

        # remove duplicates
        den_colors = list(dict.fromkeys(den_colors))

        self.signal_new_colors.emit(den_colors)

        self.axes.axhline(y=threshold)

        self.axes.set_xlabel('Detection Number')
        self.axes.set_ylabel('Distance')
        self.axes.xaxis.label.set_size(8)
        self.axes.yaxis.label.set_size(8)

        self.fig.canvas.draw()  # update matlabplot
        self.repaint()          # update widget

    def clear_plot(self):
        self.axes.clear()
        self.axes.set_title('Associations')
        self.fig.canvas.draw()
        self.repaint()

    def is_number(self, str):
        try:
            float(str)
            return True
        except ValueError:
            return False


class AssociationSettings(QFrame):

    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.buildUI()

    def buildUI(self):

        self.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)

        title_label = QLabel('Association Settings')
        title_label.setStyleSheet("font-weight: bold;")
        title_label.setAlignment(Qt.AlignCenter)

        self.threshold_edit = QDoubleSpinBox()
        self.threshold_edit.setMinimum(0.0)
        self.threshold_edit.setMaximum(1000.0)
        self.threshold_edit.setValue(5.0)

        self.dist_max_edit = QDoubleSpinBox()
        self.dist_max_edit.setMinimum(0.0)
        self.dist_max_edit.setMaximum(1000.0)
        self.dist_max_edit.setValue(10.0)

        layout = QFormLayout()
        layout.addRow(self.tr('Threshold: '), self.threshold_edit)
        layout.addRow(self.tr('Maximum Distance.: '), self.dist_max_edit)

        self.update_assoc_button = QPushButton('Update Associations')

        mainlayout = QVBoxLayout()
        mainlayout.addWidget(title_label)
        mainlayout.addLayout(layout)

        buttonLayout = QHBoxLayout()
        buttonLayout.addWidget(self.update_assoc_button)

        mainlayout.addLayout(buttonLayout)
        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        self.setLayout(mainlayout)


class MapSettings(QFrame):
    # NOT CURRENTLY USING THIS!!!!
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.buildUI()

    def buildUI(self):

        self.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)

        label_title = QLabel("Map Settings")
        label_title.setStyleSheet("font-weight: bold;")

        label_central_lon = QLabel(self.tr('Central Longitude (deg): '))
        self.central_lon_cb = QComboBox()
        self.central_lon_cb.addItem('0')
        self.central_lon_cb.addItem('180')

        label_resolution = QLabel(self.tr('Resolution'))
        self.resolution_cb = QComboBox()
        self.resolution_cb.addItem('110m')
        self.resolution_cb.addItem('50m')
        self.resolution_cb.addItem('10m')

        label_features = QLabel(self.tr('Features'))
        label_features.setStyleSheet("font-weight: bold")
        self.draw_states_check = QCheckBox('States and Provences')
        self.draw_states_check.setChecked(True)
        self.draw_lakes_check = QCheckBox('Lakes')
        self.draw_lakes_check.setChecked(True)
        self.draw_rivers_check = QCheckBox('Rivers')
        self.draw_rivers_check.setChecked(True)
        self.draw_borders_check = QCheckBox('Borders')
        self.draw_borders_check.setChecked(True)

        layout = QFormLayout()
        # layout.addRow(label_central_lon, self.central_lon_cb)
        layout.addRow(label_resolution, self.resolution_cb)
        layout.addRow(label_features)
        layout.addRow(self.draw_states_check)
        layout.addRow(self.draw_lakes_check)
        layout.addRow(self.draw_rivers_check)
        layout.addRow(self.draw_borders_check)

        mainlayout = QVBoxLayout()
        mainlayout.setAlignment(Qt.AlignCenter)
        mainlayout.addWidget(label_title)
        mainlayout.addLayout(layout)

        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        self.setLayout(mainlayout)


class Draw_Map_Worker_Object(QObject):

    thread_stopped = True
    signal_mapFinished = pyqtSignal()

    def __init__(self, figure, axes, settingsWidget):
        super().__init__()
        self.mapSettings = settingsWidget
        self.fig = figure
        self.axes = axes

    @pyqtSlot()
    def run(self):
        self.thread_stopped = False

        # self.fig.clf()
        self.axes.clear()

        resolution = self.mapSettings.resolution_cb.currentText()
        cent_lon = int(self.mapSettings.central_lon_cb.currentText())

        land = cfeature.NaturalEarthFeature('physical', 'land', resolution,
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'],
                                            linewidth=0.5)

        states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                        name='admin_1_states_provinces_lines',
                                                        scale=resolution,
                                                        facecolor='none')

        self.axes.add_feature(land)
        self.axes.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)

        self.axes.add_feature(cfeature.OCEAN.with_scale(resolution), facecolor=(22. / 255., 43. / 255., 72. / 255., 0.5))
        self.axes.add_feature(cfeature.LAKES.with_scale(resolution))
        self.axes.add_feature(cfeature.BORDERS.with_scale(resolution), linewidth=0.5)
        self.axes.add_feature(cfeature.COASTLINE.with_scale(resolution))

        self.fig.canvas.draw()  # update matlabplot

        self.signal_mapFinished.emit()

    @pyqtSlot()
    def stop(self):
        self.thread_stopped = True
