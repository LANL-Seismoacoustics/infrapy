
from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QVBoxLayout, 
                             QSplitter, QTabWidget, QAction,
                             QScrollArea, QToolBar, QToolButton)

from PyQt5.QtCore import pyqtSignal, pyqtSlot, Qt, QThread, QCoreApplication, QSettings
from PyQt5 import QtCore, QtGui
from PyQt5.QtGui import QIcon, QPainterPath, QColor, QCursor

import pyqtgraph as pg
from pyqtgraph import ViewBox

import warnings

import numpy as np
from pathlib import Path

# import infraview widgets here
from InfraView.widgets import IPDetectionWidget
from InfraView.widgets import IPDetectorSettingsWidget
from InfraView.widgets import IPNewDetectionDialog
from InfraView.widgets import IPPickLine
from InfraView.widgets import IPPlotItem
from InfraView.widgets import IPBeamformingSettingsWidget
from InfraView.widgets import IPPolarPlot
from InfraView.widgets import IPSaveBeamformingResultsDialog
from InfraView.widgets import IPUtils

# import infrapy modules here
from infrapy.detection import beamforming_new
from infrapy.utils import data_io

# import obspy modules here
from obspy.core import UTCDateTime


class IPBeamformingWidget(QWidget):

    signal_startBeamforming = pyqtSignal()
    signal_stopBeamforming = pyqtSignal()

    streams = None

    hline = None  # the horizontal crosshair line in the waveform window
    vline = None  # the vertical crosshair line in the waveform window
    position_label = None  # list to hold the labels that show the position of the crosshairs
    value_label = None  # list to hold the labels that show the y-value of the crosshairs

    _plot_list = []     # list to hold references to the four main plots

    _slowness_collection = []   # This will hold the slowness plots for the current run
    _max_projection_data = None

    _t = []
    _trace_vel = []
    _back_az = []
    _f_stats = []

    waveform_data_item = None

    _mp_pool = None

    lanl_blue = QColor(10, 44, 71)
    lanl_light_blue = QColor(34, 77, 122)
    lanl_green = QColor(67, 137, 23)
    lanl_orange = QColor(255, 144, 0)

    def __init__(self, parent, pool):
        super().__init__()

        self.parent = parent

        self._mp_pool = pool

        self.buildUI()

    def make_crosshair(self):
        crosshair = QPainterPath()
        crosshair.moveTo(0, -0.5)
        crosshair.lineTo(0,  0.5)
        crosshair.moveTo(-0.5, 0)
        crosshair.lineTo(0.5, 0)
        return crosshair

    def set_textitem_fontsize(self, item, int):
        font = item.textItem.font()
        font.setPointSize(int)
        item.textItem.setFont(font)

    def buildUI(self):

        self.make_toolbar()
        # crosshair_symbol = self.make_crosshair()

        self.lhWidget = pg.GraphicsLayoutWidget()
        self.lhWidget.setMouseTracking(True)

        self.waveformPlot = IPPlotItem.IPPlotItem(mode='waveform', est=None)
        self.waveformPlot.setLabel('left', 'Waveform')
        self.waveformPlot.hideButtons()

        self.fstatPlot = IPPlotItem.IPPlotItem(mode='waveform', est=None)
        self.fstatPlot.hideButtons()
        self.fstatPlot.setYRange(0, 1, padding=0)
        self.fstatPlot.disableAutoRange(axis=ViewBox.XAxis)
        self.fstatPlot.showGrid(x=True, y=True, alpha=0.3)
        self.fstatPlot.setLabel('left', 'F-Statistic')
        self.fstat_marker = pg.PlotDataItem([],[], symbol='+', symbolSize='25')
        self.fstatPlot.addItem(self.fstat_marker)
        self.fstat_marker_label = pg.TextItem('', color=(150,150,150), anchor=(0,1))
        self.fstat_marker_label.setZValue(15)
        self.set_textitem_fontsize(self.fstat_marker_label, 10)
        self.fstat_slowness_marker = pg.PlotDataItem([], [], symbol = 'o', symbolSize='10', color=self.lanl_blue)

        self.threshold_line = pg.InfiniteLine(pos=0.0, angle=0.0, pen=pg.mkPen('b', width=2, moveable=True, style=QtCore.Qt.DotLine))
        self.threshold_label = pg.InfLineLabel(line=self.threshold_line, text='', movable=True, position=0.04, anchors=[(0.5,1), (0.5,1)])
        self.threshold_label.setColor((0,0,255))
        t_font = self.threshold_label.textItem.font()
        t_font.setPointSize(10)
        self.threshold_label.textItem.setFont(t_font)
        self.fstatPlot.addItem(self.threshold_line)
        # this is the label that pops up to alert someone that the program is calculating the threshold
        self.threshold_calculating_label = pg.TextItem('Calculating Threshold...', color=(0,0,0))

        self.traceVPlot = IPPlotItem.IPPlotItem(mode='waveform', est=None)
        self.traceVPlot.hideButtons()
        self.traceVPlot.showGrid(x=True, y=True, alpha=0.3)
        self.traceVPlot.setYRange(0, 500, padding=0)
        self.traceVPlot.disableAutoRange(axis=ViewBox.XAxis)
        self.traceVPlot.setLabel('left', 'Trace Velocity (m/s)')
        self.traceV_marker = pg.PlotDataItem([],[], symbol='+', symbolSize='25')
        self.traceVPlot.addItem(self.traceV_marker)
        self.traceV_marker_label = pg.TextItem('', color=(150,150,150), anchor=(0,1))
        self.traceV_marker_label.setZValue(15)
        self.set_textitem_fontsize(self.traceV_marker_label, 10)
        self.traceV_slowness_marker = pg.PlotDataItem([], [], symbol = 'o', symbolSize='10', color=self.lanl_green)

        self.backAzPlot = IPPlotItem.IPPlotItem(mode='waveform', est=None)
        self.backAzPlot.hideButtons()
        self.backAzPlot.showGrid(x=True, y=True, alpha=0.3)
        self.backAzPlot.setYRange(-180, 180, padding=0)
        # I want to make sure this plot has meaningful ticks
        la = self.backAzPlot.getAxis('left')
        ba_ticks = [-180.0, -90.0, 0, 90.0, 180.0]
        la.setTicks([[(tic, str(tic)) for tic in ba_ticks]])
        self.backAzPlot.disableAutoRange(ViewBox.XAxis)
        self.backAzPlot.setLabel('left', 'Back Azimuth (deg)')
        self.backAz_marker = pg.PlotDataItem([], [], symbol='+', symbolSize='25')
        self.backAzPlot.addItem(self.backAz_marker)
        self.backAz_marker_label = pg.TextItem('', color=(150,150,150), anchor=(0,1))
        self.backAz_marker_label.setZValue(15)
        self.set_textitem_fontsize(self.backAz_marker_label, 10)
        self.backAz_slowness_marker = pg.PlotDataItem([], [], symbol = 'o', symbolSize='10', color=self.lanl_orange)

        self.resultPlots = {'fplot': self.fstatPlot, 'tracePlot': self.traceVPlot, 'backPlot': self.backAzPlot}

        self.fstatPlot.setXLink(self.waveformPlot)
        self.traceVPlot.setXLink(self.waveformPlot)
        self.backAzPlot.setXLink(self.waveformPlot)

        self.lhWidget.addItem(self.waveformPlot)
        self.lhWidget.nextRow()
        self.lhWidget.addItem(self.fstatPlot)
        self.lhWidget.nextRow()
        self.lhWidget.addItem(self.traceVPlot)
        self.lhWidget.nextRow()
        self.lhWidget.addItem(self.backAzPlot)

        self._plot_list.append(self.waveformPlot)
        self._plot_list.append(self.fstatPlot)
        self._plot_list.append(self.traceVPlot)
        self._plot_list.append(self.backAzPlot)

        self.addCrosshairs()

        # --------------------------------------------
        # this is where I create the linear region item that specifies the current portion of waveform being evaluated

        self.timeRangeLRI = pg.LinearRegionItem()
        self.timeRangeLRI.setMovable(False)
        brush = QtGui.QBrush(QtGui.QColor(50, 50, 50, 50))
        self.timeRangeLRI.setBrush(brush)

        # --------------------------------------------
        # the slownessWidget will hold the slowness plot and the projection plot

        slownessWidget = pg.GraphicsLayoutWidget()

        # Create the slowness plot and its dataitem
        self.slownessPlot = IPPolarPlot.IPPolarPlot()
        self.spi = pg.ScatterPlotItem(pxMode=False, pen=pg.mkPen(None))
        slowness_pen = pg.mkPen(color=(60, 60, 60), width=2)
        self.max_line = pg.PlotDataItem(x=[],
                                        y=[],
                                        pen=slowness_pen,
                                        symbol=None)
        self.max_line.setZValue(20)
        self.slownessPlot.addItem(self.max_line)

        # Create the slowness widget and its dataitem
        self.projectionPlot = IPPlotItem.IPPlotItem()
        
        self.projectionCurve = pg.PlotDataItem(x=[],
                                               y=[],
                                               pen=(60, 60, 60),
                                               symbol=None)
        self.max_projectionCurve = pg.PlotDataItem(x=[],
                                                   y=[],
                                                   pen=(100, 100, 100),
                                                   symbol=None)

        self.projectionPlot.showGrid(x=True, y=True, alpha=0.3)
        self.projectionPlot.addItem(self.max_projectionCurve)
        self.projectionPlot.addItem(self.projectionCurve)

        self.projectionPlot.setLabel('left', 'Avg. Beam Power')
        self.projectionPlot.setLabel('bottom', 'Azimuth')
        self.projectionPlot.setXRange(-180, 180)
        self.projectionPlot.getAxis('bottom').setTicks([[(-180, '-180'), (-90, '-90'), (0, '0'), (90, '90'), (180, '180')]])

        self.slowness_time_label = pg.LabelItem('t = ', color=QColor(44, 44, 44))
        self.slowness_backAz_label = pg.LabelItem('Back Azimuth (deg) = ', color=QColor(44, 44, 44))
        self.slowness_traceV_label = pg.LabelItem('Trace Velocity (m/s) = ', color=QColor(44, 44, 44))

        slownessWidget.addItem(self.slownessPlot)
        slownessWidget.nextRow()
        slownessWidget.addItem(self.slowness_time_label)
        slownessWidget.nextRow()
        slownessWidget.addItem(self.slowness_backAz_label)
        slownessWidget.nextRow()
        slownessWidget.addItem(self.slowness_traceV_label)
        slownessWidget.nextRow()
        slownessWidget.addItem(self.projectionPlot)


        # ---------------------------------------------
        # the bottomWidget will hold the beamforming settings widget, and the detection widget...
        bottomWidget = QWidget()
        self.bottomTabWidget = QTabWidget()

        self.bottomSettings = IPBeamformingSettingsWidget.IPBeamformingSettingsWidget(self)
        self.bottomSettings_scrollarea = QScrollArea()
        self.bottomSettings_scrollarea.setWidget(self.bottomSettings)
        
        self.detector_settings = IPDetectorSettingsWidget.IPDetectorSettingsWidget(self)
        self.detector_settings_scrollarea = QScrollArea()
        self.detector_settings_scrollarea.setWidget(self.detector_settings)
        
        self.detectionWidget = IPDetectionWidget.IPDetectionWidget(self)

        self.detectiontab_idx = self.bottomTabWidget.addTab(self.detectionWidget, 'Detections')
        self.settingstab_idx = self.bottomTabWidget.addTab(self.bottomSettings_scrollarea, 'Beamformer Settings')
        self.det_settings_tab_idx = self.bottomTabWidget.addTab(self.detector_settings, 'Detector Settings')

        bottomLayout = QHBoxLayout()
        bottomLayout.addWidget(self.bottomTabWidget)

        bottomWidget.setLayout(bottomLayout)

        # ---------------------------------------------

        self.splitterTop = QSplitter(Qt.Horizontal)
        self.splitterTop.addWidget(self.lhWidget)
        self.splitterTop.addWidget(slownessWidget)
        self.splitterBottom = QSplitter(Qt.Horizontal)
        self.splitterBottom.addWidget(bottomWidget)

        # ---------------------------------------------

        self.main_splitter = QSplitter(Qt.Vertical)
        self.main_splitter.addWidget(self.splitterTop)
        self.main_splitter.addWidget(self.splitterBottom)

        self.main_layout = QVBoxLayout()
        self.main_layout.setMenuBar(self.toolbar)
        self.main_layout.addWidget(self.main_splitter)

        self.setLayout(self.main_layout)

        self.addCrosshairs()

        self.connectSignalsAndSlots()

        # Create a thread for the beamforming to run in
        self.bfThread = QThread()

        #Temporary
        self.new_detections_dialog = IPNewDetectionDialog.IPNewDetectionsDialog(self)
        self.save_results_dialog = IPSaveBeamformingResultsDialog.IPSaveBeamformingResultsDialog(self)

    def make_toolbar(self):
        self.toolbar = QToolBar()

        # self.toolbar.setStyleSheet("QToolButton:!hover { padding-left:5px; padding-right:5px; padding-top:2px; padding-bottom:2px} QToolBar {background-color: rgb(0,107,166)}")
        # self.toolbar.setStyleSheet("QToolButton:!hover {background-color:blue} QToolButton:hover { background-color: lightgray }")

        toolButton_start = QToolButton()
        toolButton_stop = QToolButton()
        toolButton_clear = QToolButton()
        toolButton_export = QToolButton()

        self.runAct = QAction(QIcon.fromTheme("media-playback-start"), "Run Beamforming", self)
        self.runAct.triggered.connect(self.runBeamforming)
        toolButton_start.setToolButtonStyle(Qt.ToolButtonTextOnly)
        toolButton_start.setDefaultAction(self.runAct)

        self.stopAct = QAction(QIcon.fromTheme("media-playback-stop"), 'Stop', self)
        self.stopAct.setEnabled(False)
        toolButton_stop.setToolButtonStyle(Qt.ToolButtonTextOnly)
        toolButton_stop.setDefaultAction(self.stopAct)

        self.clearAct = QAction(QIcon.fromTheme("edit-clear"), 'Clear', self)
        self.clearAct.triggered.connect(self.clearResultPlots)
        toolButton_clear.setToolButtonStyle(Qt.ToolButtonTextOnly)
        toolButton_clear.setDefaultAction(self.clearAct)

        self.exportAct = QAction(QIcon.fromTheme("document-save-as"), 'Export Results', self)
        self.exportAct.triggered.connect(self.exportResults)
        toolButton_export.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        toolButton_export.setDefaultAction(self.exportAct)

        self.toolbar.addWidget(toolButton_start)
        self.toolbar.addWidget(toolButton_stop)
        self.toolbar.addWidget(toolButton_clear)
        self.toolbar.addSeparator()
        self.toolbar.addWidget(toolButton_export)


    def addCrosshairs(self):
        # This adds the crosshairs that follow the mouse around, as well as the position labels which display the
        # UTC time in the top right corner of the plots
        
        self.vline = pg.InfiniteLine(angle=90, movable=False, pen='k')
        self.hline = pg.InfiniteLine(angle=0, movable=False, pen='k')
        self.position_label = pg.TextItem(color=(0, 0, 0), html=None, anchor=(1, 0))
        self.value_label = pg.TextItem(color=(0, 0, 0), html=None, anchor=(1, 0))

        self.vline.setZValue(10)
        self.hline.setZValue(11)

        self.waveformPlot.addItem(self.vline, ignoreBounds=True)
        self.waveformPlot.addItem(self.hline, ignoreBounds=True)
        self.waveformPlot.addItem(self.position_label, ignoreBounds=True)
        self.waveformPlot.addItem(self.value_label, ignoreBounds=True)

    def connectSignalsAndSlots(self):
        # keep as many signal and slot connections as possible together in one place
        self.lhWidget.scene().sigMouseMoved.connect(self.myMouseMoved)
        self.lhWidget.scene().sigMouseClicked.connect(self.myMouseClicked)
        self.detectionWidget.signal_detections_changed.connect(self.plotDetectionLines)

    def setStreams(self, streams):
        # keep a local reference for the streams that will be analyzied
        self.streams = streams

    def get_earliest_start_time(self):
        return self.parent.waveformWidget.plotViewer.pl_widget.earliest_start_time

    def plotDetectionLines(self):
        """
        This is the routine that draws the detection lines on the fstat plot

        Plotting detection lines makes no sense if there are no waveforms loaded to set the date and
        time for the plots.  One way to check this is to see if earliest_start_time is None,
        and if it is, bail until plots are loaded.

        If detections are determined to exist, cycle through them, create a new line, connect
        it to the appropriate slots, and add it to the fstat plot.
        """
        e_s_t = self.get_earliest_start_time()
        if e_s_t is None:
            return

        # the detectioningWidget is where the detection data lives
        detection_data = self.detectionWidget.get_data()

        # Data may have changed, so first clear out old detection lines, and then
        # we'll repopulate
        self.clearDetectionLines()

        # if no detections to plot, return
        if len(detection_data) == 0:
            return

        for detection in detection_data:

            starting_position = detection.get_peakF_UTCtime(type='obspy') - UTCDateTime(e_s_t)

            newDetectionLine = IPPickLine.IPPickLine(detection, starting_pos=starting_position)

            # These connections need to be made for each new detection line
            newDetectionLine.sigPickLineMoving.connect(self.detectionWidget.detectionLineMoving)
            newDetectionLine.sigPickLineMoved.connect(self.detectionLineMoved)
            newDetectionLine.sigDeleteMe.connect(self.detectionWidget.delete_detection)

            newDetectionLine.sigCreateStartEndBars.connect(self.detectionWidget.createNewDetectionStartEndRegion)
            newDetectionLine.sigRemoveStartEndBars.connect(self.detectionWidget.removeDetectionStartEndRegion)
            newDetectionLine.sigStartEndBarsChanged.connect(self.detectionWidget.updateDetectionStartEnd)

            detection.setAssociatedPickLine(newDetectionLine)

            # add the detectionline to the fstat plot (or others eventually?) and
            # if it has one, a start/end linear region item
            self.fstatPlot.addItem(newDetectionLine)

            if newDetectionLine.startEndBars() is not None:
                self.fstatPlot.addItem(newDetectionLine.startEndBars())

    @pyqtSlot(float)
    def detectionLineMoved(self, pos):
        # somebody moved a detection line, so we need to first find the nearest _t value to the new position, and set the detection line
        # to that _t value
        nearest_idx = self.nearest_in_t(pos)
        self.detectionWidget.detectionLineMoved(self._t[nearest_idx], 
                                                self._f_stats[nearest_idx], 
                                                self._back_az[nearest_idx], 
                                                self._trace_vel[nearest_idx])
        self.plot_slowness_at_idx(nearest_idx)

        t_nearest = self._t[nearest_idx]
        f_nearest = self._f_stats[nearest_idx]
        ba_nearest = self._back_az[nearest_idx]
        tv_nearest = self._trace_vel[nearest_idx]

        self.fstat_slowness_marker.setData([t_nearest], [f_nearest])
        self.backAz_slowness_marker.setData([t_nearest], [ba_nearest])
        self.traceV_slowness_marker.setData([t_nearest], [tv_nearest])

        self.fstatPlot.addItem(self.fstat_slowness_marker)
        self.backAzPlot.addItem(self.backAz_slowness_marker)
        self.traceVPlot.addItem(self.traceV_slowness_marker)
        

    def clearDetectionLines(self):
        """
        Remove all detection lines from all plots, note that this does not remove detections
        from the list of detections in the detectionwidget
        """
        for plot in self._plot_list:
            for item in reversed(plot.items):
                if type(item) is IPPickLine.IPPickLine:
                    plot.removeItem(item)
                    del item
        self.clearDetectionStartEndRegions()

    def clearDetectionStartEndRegions(self):
        """
        Remove all start/end regions from all plots, note that this does not remove any info
        from the list of detections in the detectionwidget
        """
        for plot in self._plot_list:
            for item in reversed(plot.items):
                if type(item) is IPPickLine.IPStartEndRegionItem:
                    plot.removeItem(item)
                    del item

    @pyqtSlot(pg.PlotDataItem, tuple, str)
    def setWaveform(self, plotLine, region, plot_label=None):
        initial = False
        if self.waveform_data_item is not None:
            self.waveform_data_item.clear()
        else:
            self.waveform_data_item = pg.PlotDataItem()
            initial = True

        # bringing in a new waveform, we might have a new earliest_start_time, so update that in the 
        # plots so that the x-axes will be correct
        est = self.get_earliest_start_time()
        self.waveformPlot.setEarliestStartTime(est)
        self.fstatPlot.setEarliestStartTime(est)
        self.backAzPlot.setEarliestStartTime(est)
        self.traceVPlot.setEarliestStartTime(est)


        # need to make a copy of the currently active plot and give it to the beamformingwidget for display
        self.waveform_data_item.setData(plotLine.xData, plotLine.yData)
        self.waveform_data_item.setPen(pg.mkPen(color=(100, 100, 100), width=1))
        self.waveformPlot.enableAutoRange(axis=ViewBox.YAxis)
        if initial:
            # only need to add the item if it wasn't already added
            self.waveformPlot.addItem(self.waveform_data_item)
        if plot_label is not None:
            self.waveformPlot.setPlotLabel(plot_label)
        self.waveformPlot.setXRange(region[0], region[1], padding=0)

    @pyqtSlot(tuple)
    def updateWaveformRange(self, new_range):
        self.waveformPlot.setXRange(new_range[0], new_range[1], padding=0)
        # we want to set the title of the plot to reflect the current start time of the view
        self.start_time = self.get_earliest_start_time() + new_range[0]
        self.waveformPlot.setTitle(str(self.start_time))

    def myMouseMoved(self, evt):
        # This takes care of the crosshairs
        if len(self._plot_list) == 0:
            return

        if len(self._t) == 0:
            return

        e_s_t = self.get_earliest_start_time()
        if e_s_t is None:
            return

        # so save some cpu time, and since the plots have the same x axis, lets get the nearest time data
        # outside of the plot loop
        mouse_point_x = (self._plot_list[0].vb.mapSceneToView(evt)).x()

        nearest_idx = self.nearest_in_t(mouse_point_x)

        # get the xy values of the point nearest to the cursor
        try:
            t_nearest = self._t[nearest_idx]
            f_nearest = self._f_stats[nearest_idx]
            ba_nearest = self._back_az[nearest_idx]
            tv_nearest = self._trace_vel[nearest_idx]
        except IndexError:
            return


        self.fstat_marker_label.setText(' [{:.2f}, {:.2f}]'.format(t_nearest, f_nearest))
        self.fstat_marker_label.setPos(t_nearest, f_nearest)
        self.backAz_marker_label.setText(' [{:.2f}, {:.2f}]'.format(t_nearest, ba_nearest))
        self.backAz_marker_label.setPos(t_nearest, ba_nearest)
        self.traceV_marker_label.setText(' [{:.2f}, {:.2f}]'.format(t_nearest, tv_nearest))
        self.traceV_marker_label.setPos(t_nearest, tv_nearest)

        #set the data of all of the cursors
        self.fstat_marker.setData([t_nearest], [f_nearest])
        self.backAz_marker.setData([t_nearest], [ba_nearest])
        self.traceV_marker.setData([t_nearest], [tv_nearest])

        mouse_point_x = (self.waveformPlot.vb.mapSceneToView(evt)).x()
        mouse_point_y = (self.waveformPlot.vb.mapSceneToView(evt)).y()

        myRange = self.waveformPlot.viewRange()
        self.position_label.setPos(myRange[0][1], myRange[1][1])
        
        self.position_label.setText("UTC = {0}".format(e_s_t + mouse_point_x))
        self.vline.setPos(mouse_point_x)
        self.hline.setPos(mouse_point_y)

        #for idx, my_plot in enumerate(self._plot_list):

        #    mouse_point_y = (my_plot.vb.mapSceneToView(evt)).y()

        #    if my_plot.sceneBoundingRect().contains(evt):
        #        mouse_in_plot = True

        #        if idx == 0:
        #            self.position_label.setVisible(True)
        #            self.value_label.setVisible(True)
        #            self.position_label.setText("UTC = {0}".format(e_s_t + mouse_point_y))

                # myRange = my_plot.viewRange()
                # vb = my_plot.getViewBox()
                # _, sy = vb.viewPixelSize()  # this is to help position the valueLabels below the positionLabels

                # self.position_labels[idx].setVisible(True)
                # self.position_labels[idx].setPos(myRange[0][1], myRange[1][1])
                # self.position_labels[idx].setText("UTC = {0}".format(e_s_t + mouse_point_y))

                # self._value_labels[idx].setVisible(True)
                # self._value_labels[idx].setPos(myRange[0][1], myRange[1][1] - sy * self.position_labels[idx].boundingRect().height())
                # self._value_labels[idx].setText("{}".format(round(mouse_point_y, 4)))

            #else:
            #    self.position_labels[idx].setVisible(False)
            #    self._value_labels[idx].setVisible(False)

        #if not mouse_in_plot:
            # clear markers
        #    self.fstat_marker.setData([], [])
        #    self.backAz_marker.setData([], [])
        #    self.traceV_marker.setData([], [])

        #    self.fstat_marker_label.setText('')
        #    self.backAz_marker_label.setText('')
        #    self.traceV_marker_label.setText('')

    def myMouseClicked(self, evt):

        # if there's no data loaded, return immediately
        if len(self._f_stats) == 0:
            return

        # Gather up any keyboard modifiers to check for Ctrl, or Shift, or
        # other keypresses
        modifiers = self.parent.ipApp.keyboardModifiers()
        if modifiers == QtCore.Qt.ShiftModifier:
            # Shift+click
            pass

        elif modifiers == QtCore.Qt.ControlModifier:
            # if control click on the plot, then draw a Linear Region Item on
            # the plot
            if evt.button() == QtCore.Qt.LeftButton:
                self.mouseClick_ControlLeft(evt)
            
        elif modifiers == (QtCore.Qt.ControlModifier | QtCore.Qt.ShiftModifier):
            # Cntrl+Shift+Click
            pass

        else:
            # Handle regular Left Button Click
            if evt.button() == QtCore.Qt.LeftButton:
                # This is the primary way of adding a pick to a plot
                self.mouseClick_Left(evt)
            elif evt.button() ==QtCore.Qt.RightButton:
                self.mouseClick_right(evt)

    def mouseClick_right(self, evt):
        pass

    def mouseClick_Left(self, evt):
        # Go ahead and grab the position of the mouse click and also generate a
        # QPoint out of it for some uses
        p = QCursor.pos()  # this is the global coordinate of the mouse
        scenePos = evt.scenePos()

        for my_plot in self._plot_list:

            # screenGeometry is the global rectangle of the viewbox:
            if my_plot.vb.screenGeometry().contains(p):
                # get the index of the point nearest to the click
                mouse_point_x = (my_plot.vb.mapSceneToView(scenePos)).x()
                nearest_idx = self.nearest_in_t(mouse_point_x)

                self.plot_slowness_at_idx(nearest_idx)

                t_nearest = self._t[nearest_idx]
                f_nearest = self._f_stats[nearest_idx]
                ba_nearest = self._back_az[nearest_idx]
                tv_nearest = self._trace_vel[nearest_idx]

                self.fstat_slowness_marker.setData([t_nearest], [f_nearest])
                self.backAz_slowness_marker.setData([t_nearest], [ba_nearest])
                self.traceV_slowness_marker.setData([t_nearest], [tv_nearest])

                self.fstatPlot.addItem(self.fstat_slowness_marker)
                self.backAzPlot.addItem(self.backAz_slowness_marker)
                self.traceVPlot.addItem(self.traceV_slowness_marker)

                # move the waveform time region to reflect the location of the current selected point
                t_range = self.timeRangeLRI.getRegion()
                t_half_width = (t_range[1] - t_range[0]) / 2.
                t_region = [t_nearest - t_half_width, t_nearest + t_half_width]
                self.timeRangeLRI.setRegion(t_region)

    def mouseClick_ControlLeft(self, evt):
        # TODO: a lot of this is redundant with mouseClick_left, can they be combined in any way?

        # Go ahead and grab the position of the mouse click and also generate a
        # QPoint out of it for some uses
        p = QCursor.pos()  # this is the global coordinate of the mouse
        scenePos = evt.scenePos()

        for my_plot in self._plot_list:
            # screenGeometry is the global rectangle of the viewbox
            if my_plot.vb.screenGeometry().contains(p):
                # get the index of the point nearest to the click
                mouse_point_x = (my_plot.vb.mapSceneToView(scenePos)).x()
                nearest_idx = self.nearest_in_t(mouse_point_x)

                # plot the slowness plot for that index
                self.plot_slowness_at_idx(nearest_idx)

                t_nearest = self._t[nearest_idx]
                f_nearest = self._f_stats[nearest_idx]
                ba_nearest = self._back_az[nearest_idx]
                tv_nearest = self._trace_vel[nearest_idx]

                self.fstat_slowness_marker.setData([t_nearest], [f_nearest])
                self.backAz_slowness_marker.setData([t_nearest], [ba_nearest])
                self.traceV_slowness_marker.setData([t_nearest], [tv_nearest])

                center = self.parent.waveformWidget.stationViewer.get_current_center()
                # since we are manually adding a detection, the start and end need to be estimated...
                # lets make them +/- 5% of the window width
                window_range = self.fstatPlot.getViewBox().viewRange()
                window_width = window_range[0][1] - window_range[0][0]
                det_start = -window_width/20.0
                det_end = window_width/20.0

                det_time = self.get_earliest_start_time() + t_nearest
                dets = [[det_time, det_start, det_end, ba_nearest, tv_nearest, f_nearest]]
                self.detectionWidget.new_detections(dets,
                                    center[0],
                                    center[1],
                                    elev=center[2],
                                    event='',
                                    element_cnt=len(self.streams),
                                    method='manual',
                                    fr=self.bottomSettings.getFreqRange())

                # move the waveform time region to reflect the location of the f_max
                t_range = self.timeRangeLRI.getRegion()
                t_half_width = (t_range[1] - t_range[0]) / 2.
                t_region = [t_nearest - t_half_width, t_nearest + t_half_width]
                self.timeRangeLRI.setRegion(t_region)
                
                self.bottomTabWidget.setCurrentIndex(self.detectiontab_idx)

                

    def nearest_in_t(self, value):
        if len(self._t) < 1:
            return

        ## Return the index of the time array that is closest to value
        a = np.asarray(self._t)
        return (np.abs(a-value)).argmin()

    def getProject(self):
        return self.parent.getProject()

    def saveWindowGeometrySettings(self):
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('BeamFormingWidget')
        settings.setValue("windowSize", self.size())
        settings.setValue("windowPos", self.pos())
        settings.setValue("bfmainSplitterSettings", self.main_splitter.saveState())
        settings.setValue("splitterTopSettings", self.splitterTop.saveState())
        settings.setValue("splitterBottomSettings", self.splitterBottom.saveState())
        settings.endGroup()

    def restoreWindowGeometrySettings(self):
        # Restore settings
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('BeamFormingWidget')

        splitterTopSettings = settings.value("splitterTopSettings")
        if splitterTopSettings:
            self.splitterTop.restoreState(splitterTopSettings)

        splitterBottomSettings = settings.value("splitterBottomSettings")
        if splitterBottomSettings:
            self.splitterBottom.restoreState(splitterBottomSettings)

        splitterMainSettings = settings.value("bfmainSplitterSettings")
        if splitterMainSettings:
            self.main_splitter.restoreState(splitterMainSettings)

        settings.endGroup()

    @pyqtSlot(bool)
    def show_calculating_threshold_label(self, show):
        if show:
            xRange = self.fstatPlot.viewRange()[0]
            yRange = self.fstatPlot.viewRange()[1]
            self.fstatPlot.addItem(self.threshold_calculating_label)
            self.threshold_calculating_label.setPos(xRange[0], yRange[1])
        else:
            self.fstatPlot.removeItem(self.threshold_calculating_label)

    def runBeamforming(self):
        if self.streams is None:
            IPUtils.errorPopup('You should have at least 3 streams loaded to run beamfinder')
            return

        if len(self.streams) < 3:
            IPUtils.errorPopup('You should have at least 3 waveforms loaded to run beamfinder')
            return

        if self.parent.waveformWidget.get_inventory() is None:
            IPUtils.errorPopup('There are no stations loaded.  Station Lat and Lon information is required to do beamforming.')
            return

        if self.parent.waveformWidget.stationViewer.getStationCount() != self.streams.count():
            IPUtils.errorPopup('The number of stations is not equal to the number of waveforms. Each waveform must have a matching station with Lat./Lon. information in it.')
            return

        # we only want the slowness plot to show data at end of run
        self.spi.clear()
        self.spi.addPoints([])

        # clear out previous run
        self.clearResultPlots()
        self.fstatPlot.addItem(self.fstat_marker)
        self.fstatPlot.addItem(self.fstat_marker_label, ignoreBounds=True)

        self.traceVPlot.addItem(self.traceV_marker)
        self.traceVPlot.addItem(self.traceV_marker_label, ignoreBounds=True)

        self.backAzPlot.addItem(self.backAz_marker)
        self.backAzPlot.addItem(self.backAz_marker_label, ignoreBounds=True)

        # First lets create some new curves, and add them to the pertinent plots
        self._t = []
        self._trace_vel = []
        self._back_az = []
        self._f_stats = []
        self.slownessX = np.array([])
        self.slownessY = np.array([])
        self.beam_power = np.array([])


        self.max_projection = None
        self.max_projection_curve = None
        self.max_projection_index = None

        # add the slownessPlotItem to the slowness plot
        self.slownessPlot.addItem(self.spi)

        self.resultData = {'t': self._t,
                           'tracev': self._trace_vel,
                           'backaz': self._back_az,
                           'fstats': self._f_stats,
                           'slownessX': self.slownessX,
                           'slownessY': self.slownessY,
                           'beampower': self.beam_power}

        method = self.bottomSettings.getMethod()
        if method == 'bartlett':
            symb = 'o'
            fcolor = self.lanl_blue
            tcolor = self.lanl_green
            bcolor = self.lanl_orange
        elif method == 'gls':
            symb = '+'
            fcolor = (220, 0, 0)
            tcolor = (0, 220, 0)
            bcolor = (0, 0, 220)
        elif method == 'bartlett_covar':
            symb = 't'
            fcolor = (190, 0, 0)
            tcolor = (0, 190, 0)
            bcolor = (0, 0, 190)
        elif method == 'capon':
            symb = 's'
            fcolor = (160, 0, 0)
            tcolor = (0, 160, 0)
            bcolor = (0, 0, 160)
        elif method == 'music':
            symb = 'd'
            fcolor = (130, 0, 0)
            tcolor = (0, 130, 0)
            bcolor = (0, 0, 130)

        symbol_size = '5'

        self.fval_curve = pg.PlotDataItem(x=self._t,
                                          y=self._f_stats,
                                          pen=None,
                                          brush=fcolor,
                                          symbol=symb,
                                          symbolPen=fcolor,
                                          symbolBrush=fcolor,
                                          symbolSize=symbol_size)

        self.fval_curve.sigPointsClicked.connect(self.pointsClicked)
        self.fstatPlot.addItem(self.fval_curve)

        self.trace_curve = pg.PlotDataItem(x=self._t,
                                           y=self._trace_vel,
                                           pen=None,
                                           brush=tcolor,
                                           symbol=symb,
                                           symbolPen=tcolor,
                                           symbolBrush=tcolor,
                                           symbolSize=symbol_size)

        # self.trace_curve.sigPointsClicked.connect(self.pointsClicked)
        self.traceVPlot.addItem(self.trace_curve)
        self.backaz_curve = pg.ScatterPlotItem(x=self._t,
                                               y=self._back_az,
                                               pen=None,
                                               brush=bcolor,
                                               symbol=symb,
                                               symbolPen=bcolor,
                                               symbolBrush=bcolor,
                                               symbolSize=symbol_size)

        # self.backaz_curve.sigClicked.connect(self.pointsClicked)
        self.backAzPlot.addItem(self.backaz_curve)

        self._slowness_collection = []  # Clear this array for the new run

        # do any checks of the input here before you create the worker object.
        # The first check is to make sure the back azimuth start angle is less than the back azimuth end angle 
        baz_start, baz_end = self.bottomSettings.getBackAzRange()
        if baz_start >= baz_end:
            IPUtils.errorPopup('The back azimuth start angle must be less than the end angle. Please correct this in the Beamformer Settings tab.')
            self.bottomTabWidget.setCurrentIndex(self.settingstab_idx)
            # and bail out before going farther
            return

        # Ditto for the trace velocity range
        tv_min, tv_max = self.bottomSettings.getTraceVRange()
        if tv_min >= tv_max:
            IPUtils.errorPopup('The minimum trace velocity must be less than the max.  Please correct this in the Beamformer Settings tab.')
            self.bottomTabWidget.setCurrentIndex(self.settingstab_idx)
            return

        self.bfWorker = BeamformingWorkerObject(self.streams,
                                                self.resultData,
                                                self.bottomSettings.getNoiseRange(),
                                                self.bottomSettings.getSignalRange(),
                                                self.bottomSettings.getFreqRange(),
                                                self.bottomSettings.getWinLength(),
                                                self.bottomSettings.getWinStep(),
                                                self.bottomSettings.getMethod(),
                                                self.bottomSettings.getNumSigs(),
                                                self.bottomSettings.getSubWinLength(),
                                                self.parent.waveformWidget.get_inventory(),
                                                self._mp_pool,
                                                self.bottomSettings.getBackAzResolution(),
                                                self.bottomSettings.getTraceVelResolution(),
                                                self.bottomSettings.getTraceVRange(),
                                                self.bottomSettings.getBackAzRange(),
                                                self.detector_settings.is_auto_threshold())

        self.bfWorker.moveToThread(self.bfThread)

        self.signal_startBeamforming.connect(self.bfWorker.run)
        # self.stopButton.clicked.connect(self.bfWorker.stop)
        self.stopAct.triggered.connect(self.bfWorker.stop)

        self.bfWorker.signal_dataUpdated.connect(self.updateCurves)
        self.bfWorker.signal_slownessUpdated.connect(self.updateSlowness)
        self.bfWorker.signal_projectionUpdated.connect(self.updateProjection)
        self.bfWorker.signal_timeWindowChanged.connect(self.updateWaveformTimeWindow)
        self.bfWorker.signal_runFinished.connect(self.runFinished)
        self.bfWorker.signal_threshold_calc_is_running.connect(self.show_calculating_threshold_label)
        self.bfWorker.signal_threshold_calculated.connect(self.detector_settings.set_auto_threshold_level)
        self.bfWorker.signal_error_popup.connect(IPUtils.errorPopup)
        self.bfWorker.signal_reset_beamformer.connect(self.reset_run_buttons)

        # show the time range
        if self.timeRangeLRI not in self.waveformPlot.items:
            self.waveformPlot.addItem(self.timeRangeLRI)
        self.timeRangeLRI.setRegion((self.bottomSettings.getSignalRange()[0], self.bottomSettings.getSignalRange()[0] + self.bottomSettings.getWinLength()))

        # disable some buttons
        self.runAct.setEnabled(False)
        self.clearAct.setEnabled(False)
        self.exportAct.setEnabled(False)
        self.stopAct.setEnabled(True)

        # reset the run_step
        self.run_step = 0

        # start the thread
        self.bfThread.start()

        self.signal_startBeamforming.emit()

    def pointsClicked(self, pdi, points_clicked):
        # print('type(pdi) = {}'.format(type(pdi)))
        # print('type(points_clicked) = {}'.format(type(points_clicked)))
        # for idx, point in enumerate(points_clicked):
        #     print('{}: x={}, y={}'.format(idx, point.x(), point.y()))
        pass

    def updateCurves(self):

        self.fval_curve.setData(self._t, self._f_stats)
        self.trace_curve.setData(self._t, self._trace_vel)
        self.backaz_curve.setData(self._t, self._back_az)

        f_yrange = self.fstatPlot.vb.viewRange()[1]
        f_max = max(self._f_stats)
        if f_max > f_yrange[1]:
            self.fstatPlot.setYRange(0, f_max * 1.1, padding=0)

        t_yrange = self.traceVPlot.vb.viewRange()[1]
        t_max = max(self._trace_vel)
        if t_max > t_yrange[1]:
            self.traceVPlot.setYRange(0, t_max * 1.1, padding=0)

    @pyqtSlot(np.ndarray)
    def updateSlowness(self, slowness):
        # adds slowness to the slowness_collection
        save_slowness = True
        if save_slowness:
            self._slowness_collection.append(slowness[:])
            slowness = []


    def plot_slowness_at_idx(self, idx):
        pg.setConfigOptions(antialias=True)

        self.dots = []

        slowness = self._slowness_collection[idx]

        max_slowness = np.max(slowness[:, -1])
        max_slowness_idx = np.argmax(slowness[:,-1])
        min_slowness = np.min(slowness[:, -1])

        self.max_line.setData([0, slowness[max_slowness_idx, 0]], [0, slowness[max_slowness_idx, 1]])
        if self.max_line not in self.slownessPlot.items:
            self.slownessPlot.addItem(self.max_line)

        method = self.bottomSettings.getMethod()
        if method == "music" or method == "capon":
            max_slowness = np.max(slowness[:, -1])
            min_slowness = np.min(slowness[:, -1])
            scaled_slowness = 100 * (1.0 - (slowness[:, -1] - min_slowness) / (max_slowness - min_slowness))
        elif method == "gls":
            scaled_slowness = 100 * (1.0 - (slowness[:, -1]) / np.max(slowness[:, -1]))

        # need to auto adjust the size of the dots so they fill up the slowness nicely
        # the largest value on the plot will be 1/minimum_trace velocity
        lv = 1.0/self.bottomSettings.tracev_min_spin.value()
        dot_size = lv/15.

        for i in range(slowness.shape[0]):
            if method == "bartlett_covar" or method == "bartlett":
                brush = pg.intColor(100 * (slowness[i, -1]), hues=100, values=1)
            else:
                brush = pg.intColor(scaled_slowness[i], hues=100, values=1)

            self.dots.append({'pos': (slowness[i, 0], slowness[i, 1]), 'brush': brush, 'size': dot_size})
        
        self.spi.clear()
        self.spi.addPoints(self.dots)

        self.slowness_time_label.setText('t = {:.2f}'.format(self._t[idx]))
        self.slowness_backAz_label.setText('Back Azimuth (deg) =  {:.2f}'.format(self._back_az[idx]))
        self.slowness_traceV_label.setText('Trace Velocity (m/s) = {:.2f}'.format(self._trace_vel[idx]))
        self.slownessPlot.setRange(self.bottomSettings.tracev_min_spin.value())

        self.fstat_slowness_marker.setData([self._t[idx]], [self._f_stats[idx]])
        self.backAz_slowness_marker.setData([self._t[idx]], [self._back_az[idx]])
        self.traceV_slowness_marker.setData([self._t[idx]], [self._trace_vel[idx]])

        if self.fstat_slowness_marker not in self.fstatPlot.items:
            self.fstatPlot.addItem(self.fstat_slowness_marker)
        if self.backAz_slowness_marker not in self.backAzPlot.items:
            self.backAzPlot.addItem(self.backAz_slowness_marker)
        if self.traceV_slowness_marker not in self.traceVPlot.items:
            self.traceVPlot.addItem(self.traceV_slowness_marker)

        pg.setConfigOptions(antialias=False)

    @pyqtSlot(np.ndarray, np.ndarray)
    def updateProjection(self, projection, avg_beam_power):

        self.projectionCurve.setData(projection)
        if self.max_projection is None:
            self.max_projection = np.amax(projection[:, 1])
            self._max_projection_data = projection.copy()
            self.max_projectionCurve.setData(projection)
            self.max_projection_index = self.run_step
        else:
            _max = np.amax(projection[:, 1])
            if _max > self.max_projection:
                self.max_projection = _max
                self._max_projection_data = projection.copy()
                self.max_projectionCurve.setData(self._max_projection_data)
                self.max_projection_index = self.run_step

        method = self.bottomSettings.getMethod()
        if method == "bartlett_covar" or method == "bartlett" or method == "gls":
            self.projectionPlot.setYRange(0, 1, padding=0)
        else:
            pass

        self.projectionPlot.setXRange(-180, 180)

    @pyqtSlot(tuple)
    def updateWaveformTimeWindow(self, window):
        # this will update the linearregionitem that displays the timewindow currently evaluated
        self.timeRangeLRI.setRegion(window)
        self.run_step += 1

    @pyqtSlot()
    def reset_run_buttons(self):
        self.runAct.setEnabled(True)
        self.clearAct.setEnabled(True)
        self.exportAct.setEnabled(True)
        self.stopAct.setEnabled(False)

    @pyqtSlot()
    def runFinished(self):
        if len(self._f_stats) < 1:
            # we haven't finished a single step, so bail out
            return

        self.reset_run_buttons()

        # add a detection at the place were fstat was a maximum
        center = self.parent.waveformWidget.stationViewer.get_current_center()

        # Gather data for input into the beamforming peak detector

        # can this be done in one line?
        beam_results = np.array([self._back_az, self._trace_vel, self._f_stats])
        beam_results = beam_results.T 
        # beamforming_new uses numpy datetime64 to hold the times, so we have to 
        # convert our times to that format
        num_times = []
        for t in self._t:
            num_times.append(np.datetime64(self.get_earliest_start_time() + t))
        num_times = np.asarray(num_times)

        channel_count = len(self.streams)
        det_window_length = 300
        det_threshold = 0.99
        tb_prod = 400
        back_az_lim = 10
        
        if self.detector_settings.is_auto_threshold():
            fixed_threshold = self.detector_settings.get_auto_threshold_level()
        else:
            fixed_threshold = self.detector_settings.get_manual_threshold_level()

        self.threshold_line.setPos(fixed_threshold)
        self.threshold_label.setText('Threshold = {:.1f}'.format(fixed_threshold))
        self.fstatPlot.addItem(self.threshold_line)

        with warnings.catch_warnings(record=True) as w_array:
            dets = beamforming_new.detect_signals(num_times, 
                                                  beam_results, 
                                                  det_window_length, 
                                                  tb_prod, 
                                                  channel_count, 
                                                  det_p_val=det_threshold, 
                                                  min_seq=self.detector_settings.min_peak_width.value(), 
                                                  back_az_lim=self.detector_settings.back_az_limit.value(),
                                                  fixed_thresh=fixed_threshold)

            
            for w in w_array:
                IPUtils.errorPopup(str(w.message), "Warning")
        

        if len(dets) == 0:
            IPUtils.errorPopup("No Detections Found", "Results")
            return

        self.detectionWidget.new_detections(dets,
                                            center[0],
                                            center[1],
                                            elev=center[2],
                                            event='',
                                            element_cnt=len(self.streams),
                                            method=self.bottomSettings.getMethod(),
                                            fr=self.bottomSettings.getFreqRange())
        
        # find peak F-value location and the corresponding back azimuth and trace velocity
        f_max = max(self._f_stats)
        f_max_idx = self._f_stats.index(f_max)
        f_max_time = self._t[f_max_idx]

        # make the slowness plot show the data at the time of fstat max
        self.plot_slowness_at_idx(f_max_idx)

        # make the projection plot show the data at the time of fstat max
        if self._max_projection_data is not None:
            self.projectionCurve.setData(self._max_projection_data)

        # move the waveform time region to reflect the location of the f_max
        t_range = self.timeRangeLRI.getRegion()
        t_half_width = (t_range[1] - t_range[0]) / 2.
        t_region = [f_max_time - t_half_width, f_max_time + t_half_width]
        self.timeRangeLRI.setRegion(t_region)

        self.bottomTabWidget.setCurrentIndex(self.detectiontab_idx)

    def exportResults(self):
        if len(self._t) == 0:
            IPUtils.errorPopup("There is no data to export", "Warning")
            return  # nothing to do

        project = self.getProject()
        if project is not None:
            results_path = project.get_beamformResultsPath()
        else:
            results_path = Path.home()

        earliest_start_time = self.get_earliest_start_time()

        if self.save_results_dialog.exec_(results_path):
            filename = self.save_results_dialog.getFilename()
            t_utc = []
            for t in self._t:
                t_utc.append(earliest_start_time + t)
            
            data_io.export_beam_results_to_csv(filename, t_utc, self._f_stats, self._back_az, self._trace_vel)

            if self.save_results_dialog.wavefileIsChecked():
                # here we want to save the data that is in the visible portion of the waveform chart at the top of the beamfinder window
                wavefilename = self.save_results_dialog.getWaveFilename()
                
                xdata, ydata = self.waveform_data_item.getData()

                t_utc =[]
                for t in xdata:
                    t_utc.append(earliest_start_time + t)

                data_io.export_waveform_to_csv(wavefilename, t, ydata)
        else:
            pass

    def clearResultPlots(self):
        self.fstatPlot.clear()
        self.fstatPlot.setYRange(0, 1, padding=0)
        self._f_stats.clear()

        self.backAzPlot.clear()
        self._back_az.clear()

        self.traceVPlot.clear()
        self.traceVPlot.setYRange(0, 500, padding=0)
        self._trace_vel.clear()

        self.projectionCurve.clear()
        self.max_projectionCurve.clear()

        self.slownessPlot.clear()
        self.slownessPlot.drawPlot(self.bottomSettings.tracev_min_spin.value())
        self.spi.clear()

        self._t.clear()

        # clearing removes the crosshairs, so lets put them back
        #self.addCrosshairs()

        self.slowness_time_label.setText('t = ')
        self.slowness_backAz_label.setText('Back Azimuth (deg) = ')
        self.slowness_traceV_label.setText('Trace Velocity (m/s) = ')

    def clearWaveformPlot(self):
        self.waveform_data_item = None
        self.waveformPlot.clear()
        self.waveformPlot.setTitle(None)
        self.waveformPlot.clearPlotLabel()
        self.waveformPlot.setYRange(0, 1, padding=0)
        self.clearResultPlots()     # it doesn't make sense to have results and no waveform


class BeamformingWorkerObject(QtCore.QObject):

    signal_runFinished = pyqtSignal()
    signal_dataUpdated = pyqtSignal()
    signal_slownessUpdated = pyqtSignal(np.ndarray)
    signal_projectionUpdated = pyqtSignal(np.ndarray, np.ndarray)
    signal_timeWindowChanged = pyqtSignal(tuple)
    signal_threshold_calc_is_running = pyqtSignal(bool)
    signal_threshold_calculated = pyqtSignal(float)
    signal_error_popup = pyqtSignal(str, str)
    signal_reset_beamformer = pyqtSignal()

    def __init__(self, streams, resultData, noiseRange, sigRange, freqRange,
                 win_length, win_step, method, signal_cnt, sub_window_len,
                 inventory, pool, back_az_resol, tracev_resol, tracev_range,
                 back_az_range, auto_thresh):
        super().__init__()
        self.resultData = resultData
        self.streams = streams
        self.noiseRange = noiseRange
        self.sigRange = sigRange
        self.freqRange = freqRange
        self.win_length = win_length
        self.win_step = win_step
        self.method = method
        self.signal_cnt = signal_cnt
        self._inv = inventory
        self._pool = pool
        self._back_az_resolution = back_az_resol
        self._back_az_start = back_az_range[0]
        self._back_az_end = back_az_range[1]
        self._trace_v_resolution = tracev_resol
        self._trace_v_range = tracev_range
        self.is_auto_threshold = auto_thresh

        self.threadStopped = True

        if sub_window_len is None:
            self.sub_window_len = self.win_length
        else:
            self.sub_window_len = sub_window_len

        # Items below could be added to the beamformer settings
        self.sub_window_overlap = 0.5
        self.fft_window = 'hanning'
        self.normalize_windowing = False
        self.normalize_beam = True

    @pyqtSlot()
    def stop(self):
        self.threadStopped = True
        self.signal_reset_beamformer.emit()

    @staticmethod
    def window_beamforming_map_wrapper(args):
        return window_beamforming_map(*args)
    

    # function and wrapper to beamform different windows using pool
    # def window_beamforming(self, x, t, window, geom, delays, ns_covar_inv):
    #     X, S, f = beamforming_new.fft_array_data(x, t, window, sub_window_len=sub_window_len, sub_window_overlap=sub_window_overlap, fft_window=fft_window, normalize_windowing=normalize_windowing)
    #     beam_power = beamforming_new.run(X, S, f, geom, delays, [freq_min, freq_max], method=beam_method, ns_covar_inv=ns_covar_inv, signal_cnt=sig_cnt, normalize_beam=normalize_beam)
    #     return beamforming_new.find_peaks(beam_power, back_az_vals, trc_vel_vals, signal_cnt=sig_cnt)

    @pyqtSlot()
    def run(self):

        self.threadStopped = False

        back_az_vals = np.arange(self._back_az_start, self._back_az_end, self._back_az_resolution)
        # note if you make the 300 and 750 into a control, then you need to do that when you calculate the slowness size as well
        # TODO: the trace velocity range should be added to the settings widget
        #trc_vel_vals = np.arange(300.0, 750.0, self._trace_v_resolution)
        trc_vel_vals = np.arange(self._trace_v_range[0], self._trace_v_range[1], self._trace_v_resolution)

        det_p_val = 0.01

        latlon = []

        # we want to build the latlon array so that it has the same order as the streams
        location_count = 0
        for trace in self.streams:

            id_bits = trace.id.split('.')
            # TODO... this is a bit of a hack to help deal with horrible people who make sac files with absent network/station codes
            #         see for instance, sac_to_inventory for the other half of this
            if id_bits[0] == '':
                id_bits[0] = '###'
            if id_bits[1] == '':
                id_bits[1] = '###'
            stream_station_id = id_bits[0] + '.' + id_bits[1]

            if len(self._inv.networks) > 0:
                for network in self._inv.networks:
                    for station in network.stations:
                        station_id = network.code + '.' + station.code
                        if station_id == stream_station_id:
                            latlon.append([station.latitude, station.longitude])
                            location_count += 1
            

        #if location_count != len(self.streams):
        #    self.signal_error_popup.emit("Trace IDs don't seem to match with the inventory station list. Please check each carefully and make sure you have a matching inventory entry for each stream \n Aborting", "Inventory and Stream mismatch")
        #    self.signal_reset_beamformer.emit() # currently this will just reset the buttons
        #    return

        x, t, _, geom = beamforming_new.stream_to_array_data(self.streams, latlon)
        M, _ = x.shape

        # define slowness_grid... these are the x,y values that correspond to the beam_power values
        slowness = beamforming_new.build_slowness(back_az_vals, trc_vel_vals)
        delays = beamforming_new.compute_delays(geom, slowness)

        # Compute the noise covariance if using GLS and the detection threshold
        if self.method == "gls":
            _, S, _ = beamforming_new.fft_array_data(x, t, window=[self.noiseRange[0], self.noiseRange[1]], sub_window_len=self.sub_window_len)
            ns_covar_inv = np.empty_like(S)
            for n in range(S.shape[2]):
                S[:, :, n] += 1.0e-3 * np.mean(np.diag(S[:, :, n])) * np.eye(S.shape[0])
                ns_covar_inv[:, :, n] = np.linalg.inv(S[:, :, n])
        else:
            ns_covar_inv = None

        # #################################################### 
        # Compute detection threshold here...
        if self.is_auto_threshold:
            self.signal_threshold_calc_is_running.emit(True)
            if self._pool:
                args = []
                for window_start in np.arange(self.noiseRange[0], self.noiseRange[1], self.win_step):
                    if window_start + self.win_length > self.noiseRange[1]:
                        break

                    args = args + [[x, t, [window_start, window_start + self.win_length], geom, delays, ns_covar_inv, 
                                        self.sub_window_len, self.sub_window_overlap, self.fft_window, self.normalize_windowing, self.freqRange, 
                                        self.method, self.signal_cnt, self.normalize_beam, back_az_vals, trc_vel_vals]]

                try:
                    beam_results = np.array(self._pool.map(self.window_beamforming_map_wrapper, args))[:, 0, :]
                except IndexError:
                    IPUtils.errorPopup('Index Error...This usually occurs because the width \n'
                                    'of your noise window is less than the length of your \n'
                                    'beamforming window.  Correct that and try rerunning.')
                    self.stop()
                    return

            else:
                beam_results = []
                for window_start in np.arange(self.noiseRange[0], self.noiseRange[1], self.win_step):
                    if window_start + self.win_length > self.noiseRange[1]:
                        break
                    peaks = self.window_beamforming(x, t, [window_start, window_start + self.win_length], geom, delays, ns_covar_inv)
                    for j in range(self.signal_cnt):
                        beam_results = beam_results + [[peaks[j][0], peaks[j][1], peaks[j][2]]]
                beam_results = np.array(beam_results)

            f_vals = beam_results[:, 2] / (1.0 - beam_results[:, 2]) * (x.shape[0] - 1)
            det_thresh = beamforming_new.calc_det_thresh(f_vals, det_p_val, self.win_length * (self.freqRange[1] - self.freqRange[0]), M)

            self.signal_threshold_calc_is_running.emit(False)
            self.signal_threshold_calculated.emit(det_thresh)

        # ########## Finished calculating threshold ######################

        # Run beamforming in windowed data and write to file

        for self.window_start in np.arange(self.sigRange[0], self.sigRange[1], self.win_step):

            # In order to catch the stop button clicks, need to force process the events
            QCoreApplication.processEvents()

            if self.threadStopped:
                self.signal_runFinished.emit()
                return

            if self.window_start + self.win_length > self.sigRange[1]:
                self.signal_runFinished.emit()
                return

            self.signal_timeWindowChanged.emit((self.window_start, self.window_start + self.win_length))

            X, S, f = beamforming_new.fft_array_data(x, t, window=[self.window_start, self.window_start + self.win_length], sub_window_len=self.sub_window_len)

            beam_power = beamforming_new.run(X,
                                             S,
                                             f,
                                             geom,
                                             delays,
                                             self.freqRange,
                                             method=self.method,
                                             normalize_beam=True,
                                             signal_cnt=self.signal_cnt,
                                             pool=self._pool,
                                             ns_covar_inv=ns_covar_inv)

            # Compute relative beam power and average over frequencies
            avg_beam_power = np.average(beam_power, axis=0)

            # Analyze distribution to find peaks and compute the f-value of the peak
            peaks = beamforming_new.find_peaks(beam_power, back_az_vals, trc_vel_vals, signal_cnt=self.signal_cnt)

            self.resultData['t'].append(self.window_start + self.win_length / 2.0)
            self.resultData['backaz'].append(peaks[0][0])
            self.resultData['tracev'].append(peaks[0][1])

            sig_est, residual = beamforming_new.extract_signal(X, f, np.array([peaks[0][0], peaks[0][1]]), geom)
            signal_wvfrm = np.fft.irfft(sig_est)/(t[1]-t[0])

            if self.method == "bartlett_covar" or self.method == "bartlett" or self.method == "gls":
                fisher_val = peaks[0][2] / (1.0 - peaks[0][2]) * (M - 1)
                self.resultData['fstats'].append(fisher_val)

            else:
                self.resultData['fstats'].append(peaks[0][2])

            # Data is updated, signal plots to change
            self.signal_dataUpdated.emit()

            # Compute back azimuth projection of distribution
            az_proj, _ = beamforming_new.project_beam(beam_power, back_az_vals, trc_vel_vals)
            projection = np.c_[back_az_vals, az_proj]

            # signal projection plot to update
            self.signal_projectionUpdated.emit(projection, avg_beam_power)

            # signal slowness plot to update
            self.signal_slownessUpdated.emit(np.c_[slowness, avg_beam_power])

    # the pathos multiprocessing pool map can't pickle a QObject.  So we can't pass self to this method.
    # as a result we must make this a static method and pass all variables through the call instead of using the
    # standard "self.x" way.

def window_beamforming_map(x, 
                           t, 
                           window, 
                           geom, 
                           delays, 
                           ns_covar_inv, 
                           sub_win_len, 
                           sub_win_over, 
                           fft_win, 
                           norm_win,
                           f_range,
                           method,
                           sig_count,
                           norm_beam,
                           back_az_vals,
                           trace_vel_vals):

    X, S, f = beamforming_new.fft_array_data(x,
                                             t,
                                             window,
                                             sub_window_len=sub_win_len, 
                                             sub_window_overlap=sub_win_over, 
                                             fft_window=fft_win, 
                                             normalize_windowing=norm_win)
        
    beam_power = beamforming_new.run(X, 
                                     S, 
                                     f, 
                                     geom,
                                     delays, 
                                     f_range, 
                                     method=method, 
                                     ns_covar_inv=ns_covar_inv, 
                                     signal_cnt=sig_count, 
                                     normalize_beam=norm_beam)
       
    return beamforming_new.find_peaks(beam_power, back_az_vals, trace_vel_vals, signal_cnt=sig_count)

    