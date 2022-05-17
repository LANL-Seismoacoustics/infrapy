import pyqtgraph as pg
import numpy as np

from PyQt5 import QtCore
from PyQt5.QtCore import Qt, pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import (QWidget, QDoubleSpinBox,
                             QLabel, QMessageBox,
                             QHBoxLayout, QVBoxLayout,
                             QScrollArea, QSplitter)

from InfraView.widgets import IPPlotWidget
from InfraView.widgets import IPWaveformSelectorWidget

from obspy.core import UTCDateTime
from obspy.core.stream import Stream


class IPPlotViewer(QSplitter):

    def __init__(self, parent, fs_widget):
        super().__init__(parent)

        self._parent = parent
        self.settings = parent.settings
        self.fs_widget = fs_widget
        self.buildUI()

    def buildUI(self):

        self.pl_widget = IPPlotLayoutWidget(self)
        self.waveform_selector = IPWaveformSelectorWidget.IPWaveformSelectorWidget(self)
        self.title = QLabel("")
        waveform_selector_scrollarea = QScrollArea()
        waveform_selector_scrollarea.setWidget(self.waveform_selector)
        self.lr_settings_widget = IPLinearRegionSettingsWidget(self)
        
        rhs_widget = QWidget()
        rhs_layout = QVBoxLayout()

        rhs_layout.addWidget(self.title)
        rhs_layout.setAlignment(self.title, Qt.AlignHCenter)
        rhs_layout.addWidget(self.pl_widget)
        rhs_layout.addWidget(self.lr_settings_widget)
        rhs_widget.setLayout(rhs_layout)

        self.addWidget(self.waveform_selector)
        self.addWidget(rhs_widget)

    def set_streams(self, st, st_filtered, c_f_d_s):
        # c_f_d_s is the current filter display settings
        self.waveform_selector.update_selections(st)
        self.pl_widget.plot_traces(st, st_filtered, c_f_d_s)

    def clear(self):
        self.pl_widget.sts = None
        self.pl_widget.plot_list.clear()
        self.pl_widget.filtered_plot_lines.clear()
        self.pl_widget.clear()
        self.waveform_selector.clear_form()
        self.title.setText("")

    @pyqtSlot(Stream, Stream)
    def update(self, sts, filtered_sts):
        pass

    @pyqtSlot(dict)
    def show_hide_lines(self, current_filter_display_settings):
        self.pl_widget.draw_plot_lines(current_filter_display_settings)

    def get_plot_lines(self):
        return self.pl_widget.plot_lines

    def get_filtered_plot_lines(self):
        return self.pl_widget.filtered_plot_lines


class IPPlotLayoutWidget(pg.GraphicsLayoutWidget):

    sig_active_plot_changed = pyqtSignal(int, list, list, tuple)

    plot_list = []              # this list will hold references to the plots
    plot_lines = []             # this list will hold references to the unfiltered plot lines
    filtered_plot_lines = []    # this list will hold references to the filtered plot lines

    t = []  # this will hold the list of time series for all the plots
    freqs = [] # this will hold the frequency list for the spectrograms

    earliest_start_time = None
    latest_end_time = None

    v_lines = []
    h_lines = []
    position_labels = []

    pick_line_list = []  # array to hold the references to pick lines

    active_plot = 0
    last_range = []

    def __init__(self, parent):
        super().__init__(parent=parent)
        self._parent = parent
        self.settings = parent.settings     # for convenience

        self.connect_signals_and_slots()
        self.setMouseTracking(True)

    def errorPopup(self, message):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setText(message)
        msgBox.setWindowTitle("Oops...")
        msgBox.exec_()

    def connect_signals_and_slots(self):
        self.scene().sigMouseMoved.connect(self.myMouseMoved)
        self.scene().sigMouseClicked.connect(self.myMouseClicked)

    def get_active_plot(self):
        return self.active_plot

    def calc_spectrograms(self, streams):
        self.freqs.clear()
        self.times.clear()
        self.hist.clear()
        self.spectrograms.clear()

        self.sts = streams

        if streams is None:
            return

        for trace in streams:
            freqs, times, spectrogram = signal.spectrogram(trace.data, fs=trace.stats['sampling_rate'])
            self.freqs.append(freqs)
            self.times.append(times)
            self.spectrograms.append(np.log10(spectrogram))

        #self.update_images()

    def plot_traces(self,
                    sts,
                    filtered_sts,
                    current_filter_display_settings):

        pg.setConfigOptions(antialias=False)

        self.plot_list.clear()
        self.plot_lines.clear()
        self.filtered_plot_lines.clear()
        self.position_labels.clear()
        self.t.clear()
        self.clearAll()

        self.active_plot = 0

        if sts is None:
            # nothing to do
            return

        # populate self.t
        # the filtered lines will have the same range as the unfiltered, so just pass sts
        self.getAllTraceTimeSeries(sts)

        # title will always be at 0,0 if it exists
        self._parent.title.setText(self.earliest_start_time.isoformat())

        for idx, trace in enumerate(sts):

            #####################################
            # First let's assemble the trace data, raw and, if needed, filtered

            # trace.data is the unfiltered data for the current trace.  Lets create a plot line for it
            self.plot_lines.append(pg.PlotDataItem(self.t[idx], trace.data))
            self.plot_lines[idx].setPen(pg.mkPen(color=(100, 100, 100), width=1))
            self.plot_lines[idx].setZValue(5)

            # generate the filtered lines from the filtered streams
            self.filtered_plot_lines.append(pg.PlotDataItem(self.t[idx], filtered_sts[idx].data))

            # create a new plot and add it to the layout
            new_plot = IPPlotWidget.IPPlotWidget(mode='waveform')
            new_plot.addItem(self.plot_lines[idx], name=trace.id)
            new_plot.addItem(self.filtered_plot_lines[idx])

            ###########################################################################
            # Now lets set up the signal and noise linear region items
            #
            # connect signals and slots so that if the lris are changed, the spins change, and vice versa
            self._parent.lr_settings_widget.noiseSpinsChanged.connect(new_plot.getNoiseRegion().setRegion)
            self._parent.lr_settings_widget.signalSpinsChanged.connect(new_plot.getSignalRegion().setRegion)
            new_plot.getSignalRegion().sigRegionChanged.connect(self._parent.lr_settings_widget.updateSpinValues)
            new_plot.getNoiseRegion().sigRegionChanged.connect(self._parent.lr_settings_widget.updateSpinValues)

            # itialize the positions of the linear region items
            tstart = self.t[idx][0]
            tend = self.t[idx][-1]
            li_range = (tend - tstart) / 10.0

            new_plot.setSignalRegionRange([tend - 2 * li_range, tend - li_range])
            new_plot.setNoiseRegionRange([tstart + li_range, tstart + 2 * li_range])

            # This is the bit that links the signal and noise regions so that when you move one, they all move
            if idx > 0:
                new_plot.getSignalRegion().sigRegionChanged.connect(self.plot_list[0].copySignalRange)
                new_plot.getNoiseRegion().sigRegionChanged.connect(self.plot_list[0].copyNoiseRange)
                self.plot_list[0].getSignalRegion().sigRegionChanged.connect(new_plot.copySignalRange)
                self.plot_list[0].getNoiseRegion().sigRegionChanged.connect(new_plot.copyNoiseRange)

            new_plot.getNoiseRegion().sigRegionChanged.connect(self._parent._parent.update_noise_PSD)
            new_plot.getSignalRegion().sigRegionChanged.connect(self._parent._parent.update_signal_PSD)

            ##################################################################################
            #  Now we can initialize the background colors of the plots
            #
            # we want to highlight the active plot, so set the background color here
            # if idx == self.active_plot:
            #     new_plot.setBackgroundColor(255, 255, 255)
            # else:
            #     new_plot.setBackgroundColor(200, 200, 200)

            # cluge because setting background color covers axis for some reason
            new_plot.getAxis("top").setZValue(0)
            new_plot.getAxis("bottom").setZValue(0)
            new_plot.getAxis("left").setZValue(0)
            new_plot.getAxis("right").setZValue(0)
            #
            ####################################################################################
            # create the legend/label for the new plot
            li = pg.LegendItem()
            li.setParentItem(new_plot.vb)
            li.anchor(itemPos=(0, 1), parentPos=(0, 1))
            li.addItem(self.plot_lines[idx], trace.id)

            # set the plot lines' color and width (if you chance something here, it needs to also
            # be changed in updateTraces)

            # keep the new_plot reference in a list
            self.plot_list.append(new_plot)

        # Now set up the crosshairs and position labels for each plot. 
        self.v_lines.clear()
        self.h_lines.clear()
        for idx, my_plot in enumerate(self.plot_list):
            self.v_lines.append(pg.InfiniteLine(angle=90, movable=False, pen='k'))
            self.h_lines.append(pg.InfiniteLine(angle=0, movable=False, pen='k'))
            self.v_lines[idx].setZValue(10)
            self.h_lines[idx].setZValue(11)

            self.position_labels.append(pg.TextItem(color=(0, 0, 0), html=None, anchor=(1, 0)))

            my_plot.addItem(self.v_lines[idx], ignoreBounds=True)
            my_plot.addItem(self.h_lines[idx], ignoreBounds=True)
            my_plot.addItem(self.position_labels[idx], ignoreBounds=True)

        # now draw the plots...
        self.draw_plots()
        # and finally draw the data lines on the plots
        self.draw_plot_lines(current_filter_display_settings)
        # updateAxes currently only links the x-axes
        self.updateAxes()
        # signal to the psd viewer that the plots are loaded and to draw the initial psd
        self.sig_active_plot_changed.emit(self.active_plot,
                                              self.plot_lines,
                                              self.filtered_plot_lines,
                                              self.plot_list[self.active_plot].getSignalRegion().getRegion())

    def draw_plots(self):
        # The idea here is to draw the plots that are currently selected in the waveformselectorwidget

        # first step is to clear out the current layout
        self.clear()

        # now, iterate through the current checked items and add the plots that are at each index
        values = self._parent.waveform_selector.get_value_list()

        # if there are no checked plots...
        if not any(values):
            self.sig_active_plot_changed.emit(-1,
                                              self.plot_lines,
                                              self.filtered_plot_lines,
                                              self.plot_list[0].getSignalRegion().getRegion())
            return

        # if the active_plot is not currently checked, make the first checked plot the active plot
        if not values[self.active_plot]:
            self.active_plot = values.index(True)
            self.sig_active_plot_changed.emit(self.active_plot,
                                              self.plot_lines,
                                              self.filtered_plot_lines,
                                              self.plot_list[self.active_plot].getSignalRegion().getRegion())

        for idx, plot in enumerate(self.plot_list):
            # set the color of the plots
            if idx == self.active_plot:
                plot.setBackgroundColor(255, 255, 255)
            else:
                plot.setBackgroundColor(230, 230, 230)
            # add the checked plots in the waveformselector to the layout
            if values[idx]:
                self.nextRow()
                self.addItem(plot)

    def draw_plot_lines(self, current_filter_display_settings):

        # Now we need to make sure the correct plot lines are turned on/off
        if current_filter_display_settings['apply']:
            for idx, _ in enumerate(self.plot_lines):
                self.filtered_plot_lines[idx].setVisible(True)
                if current_filter_display_settings['showUnfiltered']:
                    self.plot_lines[idx].setVisible(True)
                else:
                    self.plot_lines[idx].setVisible(False)
        else:
            for idx, _ in enumerate(self.plot_lines):
                self.filtered_plot_lines[idx].setVisible(False)
                self.plot_lines[idx].setVisible(True)

        # Make sure they have correct colors/linestyles etc
        # this could be contained in the setVisible bit above, but I separate
        # them for clarity
        if current_filter_display_settings['apply']:
            for idx, _ in enumerate(self.plot_lines):
                self.filtered_plot_lines[idx].setPen(pg.mkPen(color=(100, 100, 100), width=1))
                self.plot_lines[idx].setPen(pg.mkPen(color=(220, 220, 220), width=1))

                self.filtered_plot_lines[idx].setZValue(5)
                self.plot_lines[idx].setZValue(4)
        else:
            for idx, _ in enumerate(self.plot_lines):
                self.plot_lines[idx].setPen(pg.mkPen(color=(100, 100, 100), width=1))
                self.plot_lines[idx].setZValue(5)

    def update_filtered_line_data(self, filtered_streams):
        for idx, trace in enumerate(filtered_streams):
            self.filtered_plot_lines[idx].setData(self.t[idx], trace.data)

    def updateAxes(self):
        # This little bit links the axes of the newly created plot to those of the first one
        # so they scale the same (this could be optional)
        self.setXaxisCoupling(True)
        # self.setYaxisCoupling(True)

        #if len(self.plot_list) > 0 and self.latest_end_time is not None and self.earliest_start_time is not None:
        #    self.plot_list[-1].setRange(xRange=(0, UTCDateTime(self.latest_end_time) - UTCDateTime(self.earliest_start_time)), padding=0)

    def normalize_traces_to_largest(self, stream):
        my_max = 0
        max_i = 0
        for idx, trace in enumerate(stream):
            m = abs(max(trace.data - np.mean(trace.data)))
            if m > my_max:
                my_max = m
                max_i = idx
        self.plot_list[max_i].getViewBox().autoRange(padding=0)
        

    @pyqtSlot(object, object)
    def adjustSignalRegionRange(self, dummy, new_range):
        if len(self.plot_list) > 0:
            self.plot_list[0].getSignalRegion().setRegion(new_range)

    def getAllTraceTimeSeries(self, sts):
        self.earliest_start_time = None
        self.latest_end_time = None
        for trace in sts:
            b = 0.0
            if trace.stats['_format'] == 'SAC':
                b = trace.stats.get('sac').get('b')

            if self.earliest_start_time is None:
                self.earliest_start_time = trace.stats.starttime
            else:
                if UTCDateTime(trace.stats.starttime) < UTCDateTime(self.earliest_start_time):
                    self.earliest_start_time = trace.stats.starttime

            if self.latest_end_time is None:
                self.latest_end_time = trace.stats.endtime
            else:
                if UTCDateTime(self.latest_end_time) < UTCDateTime(trace.stats.endtime):
                    self.latest_end_time = trace.stats.endtime

        # Now find the offset start for each trace
        offsets = []
        for trace in sts:
            # if the trace is from a SAC file, find out if there is a b value
            b = 0.0
            if trace.stats['_format'] == 'SAC':
                b = trace.stats.get('sac').get('b')
                
            offsets.append(UTCDateTime(trace.stats.starttime) -
                           UTCDateTime(self.earliest_start_time))

        self.t.clear()
        for idx, trace in enumerate(sts):
            N = trace.stats.npts
            dt = trace.stats.delta
            self.t.append(offsets[idx] + np.arange(0, N) * dt)

    # ---------------------------------------------------------------------------

    def setXaxisCoupling(self, coupled, index=0):
        if coupled:
            for idx, my_plot in enumerate(self.plot_list):
                if idx != index:
                    my_plot.setXLink(self.plot_list[index])
        else:
            for my_plot in self.plot_list:
                my_plot.setXLink(None)

    def setYaxisCoupling(self, coupled, index=0):
        if coupled:
            for idx, my_plot in enumerate(self.plot_list):
                if idx != index:
                    my_plot.setYLink(self.plot_list[index])
        else:
            for my_plot in self.plot_list:
                my_plot.setYLink(None) 
                my_plot.enableAutoRange(axis=my_plot.vb.YAxis, enable=True)

    def clearAll(self):

        self.plot_lines.clear()
        self.filtered_plot_lines.clear()
        self.plot_list.clear()
        self.t.clear()
        self.v_lines.clear()
        self.h_lines.clear()
        self.active_plot = 0

    def clearLayout(self):
        for plot in reversed(self.plot_list):
            self.removeItem(plot)

    # -----------------------------------------------------------------------------
    # Mouse Handlers

    # Mouse move handler
    def myMouseMoved(self, evt):
        # This takes care of the crosshairs
        pos = evt
        mirrorCrosshairs = True

        if self.plot_list is not None and len(self.plot_list) > 0:
            for idx, my_plot in enumerate(self.plot_list):
                try:
                    mousePoint = my_plot.vb.mapSceneToView(pos)
                except Exception:
                    return
                self.v_lines[idx].setPos(mousePoint.x())
                self.h_lines[idx].setPos(mousePoint.y())

                if my_plot.sceneBoundingRect().contains(pos):
                    self.v_lines[idx].setVisible(True)
                    self.h_lines[idx].setVisible(True)
                    self.position_labels[idx].setVisible(True)

                    myRange = my_plot.viewRange()
                    self.position_labels[idx].setPos(myRange[0][1], myRange[1][1])
                    try:
                        self.position_labels[idx].setText("UTC = {0}".format(self.earliest_start_time + mousePoint.x()))
                    except Exception:
                        # it's possible to zoom out so far the the mousepoint will be so large when you add it t the start time
                        # obspy will throw a Value Error.  This is fine, just return and don't update the labels
                        return
                else:
                    self.position_labels[idx].setVisible(False)

                    if mirrorCrosshairs:
                        self.v_lines[idx].setVisible(True)
                        self.h_lines[idx].setVisible(True)

                    else:
                        self.v_lines[idx].setVisible(False)
                        self.h_lines[idx].setVisible(False)

    ########################################################
    # Mouse click handlers

    def myMouseClicked(self, evt):

        # if there's no data loaded, return immediately
        if len(self.plot_list) < 1:
            return

        if evt.button() == QtCore.Qt.LeftButton:
            self.mouseClick_Left(evt)

    def mouseClick_Left(self, evt):

        scenePos = evt.scenePos()

        for idx, my_plot in enumerate(self.plot_list):
            # we want to ignore clicks if they aren't actually in the plot area (ie ignore axis labels etc)
            mySceneBoundingRect = my_plot.sceneBoundingRect()

            if mySceneBoundingRect.contains(scenePos):
                self.active_plot = idx
                my_plot.setBackgroundColor(255, 255, 255)
                self.sig_active_plot_changed.emit(idx,
                                                  self.plot_lines,
                                                  self.filtered_plot_lines,
                                                  self.plot_list[0].getSignalRegion().getRegion())

            else:
                my_plot.setBackgroundColor(230, 230, 230)

            # cluge because setting background color covers axis for some reason
            my_plot.getAxis("top").setZValue(0)
            my_plot.getAxis("bottom").setZValue(0)
            my_plot.getAxis("left").setZValue(0)
            my_plot.getAxis("right").setZValue(0)
        
        # TODO

        # self.updateWidgets()

    # ------------------------------------------------------------------------------
    # Key press events...

    def keyPressEvent(self, evt):
        pass

    def myNoiseSpinsChanged(self):
        start = self.noiseStartSpin.value()
        stop = start + self.noiseDurationSpin.value()
        self.noiseSpinsChanged.emit((start, stop))

    def mySignalSpinsChanged(self):
        start = self.signalStartSpin.value()
        stop = start + self.signalDurationSpin.value()
        self.signalSpinsChanged.emit((start, stop))


class IPLinearRegionSettingsWidget(QWidget):

    noiseSpinsChanged = pyqtSignal(tuple)
    signalSpinsChanged = pyqtSignal(tuple)

    def __init__(self, parent):
        super().__init__()

        self._parent = parent

        layout = QHBoxLayout(self)

        self.noiseStartSpin = QDoubleSpinBox()
        self.noiseStartSpin.setMinimum(0)
        self.noiseStartSpin.setMaximum(1000000)
        self.noiseStartSpin.setSuffix(' s')
        self.noiseDurationSpin = QDoubleSpinBox()
        self.noiseDurationSpin.setMinimum(0)
        self.noiseDurationSpin.setMaximum(1000000)
        self.noiseDurationSpin.setSuffix(' s')

        self.noiseStartSpin.valueChanged.connect(self.myNoiseSpinsChanged)
        self.noiseDurationSpin.valueChanged.connect(self.myNoiseSpinsChanged)

        self.signalStartSpin = QDoubleSpinBox()
        self.signalStartSpin.setMinimum(0)
        self.signalStartSpin.setMaximum(1000000)
        self.signalStartSpin.setSuffix(' s')

        self.signalDurationSpin = QDoubleSpinBox()
        self.signalDurationSpin.setMinimum(0)
        self.signalDurationSpin.setMaximum(1000000)
        self.signalDurationSpin.setSuffix(' s')

        self.signalStartSpin.valueChanged.connect(self.mySignalSpinsChanged)
        self.signalDurationSpin.valueChanged.connect(self.mySignalSpinsChanged)

        layout.addWidget(QLabel('Noise Region Start: '))
        layout.addWidget(self.noiseStartSpin)
        layout.addWidget(QLabel('Duration: '))
        layout.addWidget(self.noiseDurationSpin)
        layout.addStretch(1)
        layout.addWidget(QLabel('Signal Region Start: '))
        layout.addWidget(self.signalStartSpin)
        layout.addWidget(QLabel('Duration: '))
        layout.addWidget(self.signalDurationSpin)

    def updateSpinValues(self, regionItem):
        nrange = regionItem.getRegion()

        # All of the regions are linked, so pull off the values of the first one
        if type(regionItem) is IPPlotWidget.IPLinearRegionItem_Noise:
            self.noiseStartSpin.setValue(nrange[0])
            self.noiseDurationSpin.setValue(nrange[1] - nrange[0])

        elif type(regionItem) is IPPlotWidget.IPLinearRegionItem_Signal:
            self.signalStartSpin.setValue(nrange[0])
            self.signalDurationSpin.setValue(nrange[1] - nrange[0])

    def myNoiseSpinsChanged(self):
        start = self.noiseStartSpin.value()
        stop = start + self.noiseDurationSpin.value()
        self.noiseSpinsChanged.emit((start, stop))

    def mySignalSpinsChanged(self):
        start = self.signalStartSpin.value()
        stop = start + self.signalDurationSpin.value()
        self.signalSpinsChanged.emit((start, stop))