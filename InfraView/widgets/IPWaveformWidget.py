import pyqtgraph as pg
import numpy as np

from PyQt5 import QtCore
from PyQt5.QtCore import Qt, pyqtSlot, QSettings
from PyQt5.QtWidgets import (QWidget, QGridLayout, QSplitter, QTabWidget)

from InfraView.widgets import (IPFilterSettingsWidget,
                               IPPlotViewer,
                               IPPSDWidget,
                               IPStationView,
                               IPStatsView)

import copy

import obspy
from obspy.core.stream import Stream
from obspy.core.inventory import Inventory, Network, Station, Channel, Site

from InfraView.widgets import IPUtils


class IPWaveformWidget(QWidget):

    """ The IPWaveformWidget holds the waveform and inventory data.  The IPPlotViewer plots the data,
    The IPFilterSettingsWidget holds the filter settings and tells the WaveformWidget when to update that
    data.  The IPStatsView displays the trace data, and the IPStationView displays the station data.
    """

    _sts = None             # streams
    _sts_filtered = None    # filtered streams
    _inv = None             # inventory

    def __init__(self, parent=None, pool=None, project=None):
        super().__init__(parent)

        self.parent = parent
        self._mp_pool = pool

        self.buildUI()

    def buildUI(self):

        self.stationViewer = IPStationView.IPStationView(self)
        self.statsViewer = IPStatsView.IPStatsView(self)
        self.info_tabs = QTabWidget()
        self.info_tabs.addTab(self.stationViewer, 'Station Info')
        self.info_tabs.addTab(self.statsViewer, 'Trace Info')

        self.filterSettingsWidget = IPFilterSettingsWidget.IPFilterSettingsWidget(self)
        self.spectraWidget = IPPSDWidget.IPPSDWidget(self)

        self.plotViewer = IPPlotViewer.IPPlotViewer(self)

        self.lh_splitter = QSplitter(Qt.Vertical)
        self.lh_splitter.setStyleSheet("QSplitter::handle{ background-color: #DDD}")
        self.lh_splitter.addWidget(self.plotViewer)
        self.lh_splitter.addWidget(self.info_tabs)

        self.rh_splitter = QSplitter(Qt.Vertical)
        self.rh_splitter.setStyleSheet("QSplitter::handle{ background-color: #DDD}")
        self.rh_splitter.addWidget(self.spectraWidget)
        self.rh_splitter.addWidget(self.filterSettingsWidget)

        self.main_splitter = QSplitter(Qt.Horizontal)
        self.main_splitter.setStyleSheet("QSplitter::handle{ background-color: #DDD}")
        self.main_splitter.addWidget(self.lh_splitter)
        self.main_splitter.addWidget(self.rh_splitter)

        main_layout = QGridLayout()
        main_layout.addWidget(self.main_splitter)

        self.setLayout(main_layout)

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.filterSettingsWidget.sig_filter_changed.connect(self.update_filtered_data)
        self.filterSettingsWidget.sig_filter_display_changed.connect(self.plotViewer.show_hide_lines)

        self.statsViewer.removeTrace.connect(self.remove_trace)

        self.plotViewer.lr_settings_widget.noiseSpinsChanged.connect(self.parent.beamformingWidget.bottomSettings.setNoiseValues)
        self.plotViewer.lr_settings_widget.signalSpinsChanged.connect(self.parent.beamformingWidget.bottomSettings.setSignalValues)
        self.plotViewer.lr_settings_widget.signalSpinsChanged.connect(self.parent.beamformingWidget.updateWaveformRange)
        self.plotViewer.lr_settings_widget.signalSpinsChanged.connect(self.parent.singleSensorWidget.updateSignalRange)
        self.plotViewer.lr_settings_widget.noiseSpinsChanged.connect(self.parent.singleSensorWidget.updateNoiseRange)
        self.plotViewer.pl_widget.sig_active_plot_changed.connect(self.update_widgets)

        self.spectraWidget.f1_Spin.valueChanged.connect(self.parent.beamformingWidget.bottomSettings.setFmin)
        self.spectraWidget.f2_Spin.valueChanged.connect(self.parent.beamformingWidget.bottomSettings.setFmax)
        self.spectraWidget.psdPlot.getFreqRegion().sigRegionChanged.connect(self.parent.beamformingWidget.bottomSettings.setFreqValues)

    def get_project(self):
        return self.parent.getProject()

    @pyqtSlot(obspy.core.stream.Stream, obspy.core.inventory.inventory.Inventory)
    def appendTraces(self, newTraces, newInventory):
        print("appending")
        if newTraces is None:
            return

        if self._sts is None:
            self._sts = newTraces
        else:
            self._sts += newTraces

        self.update_inventory(newInventory)

        for trace in self._sts:
            trace.data = trace.data - np.mean(trace.data)
            self._sts.merge(fill_value=0)

        print("Merged traces...")

        # it's possible, if the open failed, that self.waveformWidget._sts is still None, so if it is, bail out
        # if not populate the trace stats viewer and plot the traces
        if self._sts is not None:
            # TODO...is there a better way of doing this?
            self.parent.beamformingWidget.setStreams(self._sts)
            self.stationViewer.setInventory(self._inv)
            self.statsViewer.setStats(self._sts)

            self.update_streams(self._sts)

            self.parent.setStatus("Ready", 5000)

            print('finished appending')
        else:
            return

    @pyqtSlot(obspy.core.stream.Stream, obspy.core.inventory.inventory.Inventory)
    def replaceTraces(self, newTraces, newInventory):
        # same as append, just clear out the old traces and inventory first
        self._sts = None
        self._inv = None
        self.stationViewer.setInventory(self._inv)

        self.appendTraces(newTraces, newInventory)

    @pyqtSlot(Inventory)
    def update_inventory(self, new_inventory):
        if self._inv is None:
            self._inv = new_inventory
        else:
            self._inv += new_inventory
        self.stationViewer.setInventory(self._inv)

    def remove_from_inventory(self, net, sta, loc, cha, keep_empty=False):
        new_inventory = self.inv_remove(self._inv, network=net, station=sta, location=loc, channel=cha, keep_empty=keep_empty)
        self.set_inventory(new_inventory)

    def get_inventory(self):
        return self._inv

    def set_inventory(self, new_inv):
        self._inv = None
        self.update_inventory(new_inv)

    def clear_inventory(self):
        self._inv = None
        self.stationViewer.clear()

    def get_streams(self):
        return self._sts

    def get_filtered_streams(self):
        return self._sts_filtered

    def getTraceName(self, trace):
        traceName = trace.stats['network'] + '.' + trace.stats['station'] + \
            '.' + trace.stats['location'] + '.' + trace.stats['channel']
        return traceName

    def get_earliest_start_time(self):
        return self.plotViewer.pl_widget.earliest_start_time

    @pyqtSlot(Stream)
    def update_streams(self, new_stream):
        # this should be called when you load new streams, or remove traces
        self._sts = new_stream
        self._sts_filtered = self.filter_stream(self._sts,
                                                self.filterSettingsWidget.get_filter_settings())
        self.plotViewer.set_streams(self._sts,
                                    self._sts_filtered,
                                    self.filterSettingsWidget.get_filter_display_settings())

        self.statsViewer.setStats(new_stream)

    def debug_trace(self):  # for debugging, you have to call pyqtRemoveInputHook before set_trace()
        from PyQt5.QtCore import pyqtRemoveInputHook
        from pdb import set_trace
        pyqtRemoveInputHook()
        set_trace()

    @pyqtSlot(dict)
    def update_filtered_data(self, filter_settings):

        # this should be called when settings in the filter widget are changed
        if self._sts is None:
            # Nothing to filter, clear out the filtered_streams and return
            self._sts_filtered = None
            return

        self._sts_filtered = self.filter_stream(self._sts,
                                                filter_settings)

        self.plotViewer.pl_widget.update_filtered_line_data(self._sts_filtered)
        index = self.plotViewer.pl_widget.get_active_plot()
        self.update_widgets(index, 
                            self.plotViewer.get_plot_lines(), 
                            self.plotViewer.get_filtered_plot_lines(), 
                            self.plotViewer.pl_widget.plot_list[index].getSignalRegionRange(),
                            self.plotViewer.pl_widget.plot_list[index].getNoiseRegionRange())

    def filter_stream(self, stream, cfs):

        # cfs: Current Filter Settings
        if stream is None:
            # nothing to do
            return None

        filtered_stream = Stream()

        for trace in stream:

            filtered_trace = trace.copy()

            filtType = cfs['type']

            if filtType == 'High Pass':
                try:
                    filtered_trace.filter('highpass',
                                        freq=cfs['F_high'],
                                        corners=cfs['order'],
                                        zerophase=cfs['zphase'])
                except ValueError as e:
                    IPUtils.errorPopup(str(e))

            elif filtType == 'Low Pass':
                try:
                    filtered_trace.filter('lowpass',
                                        freq=cfs['F_low'],
                                        corners=cfs['order'],
                                        zerophase=cfs['zphase'])
                except ValueError as e:
                    IPUtils.errorPopup(str(e))

            elif filtType == 'Band Pass':
                try:
                    filtered_trace.filter('bandpass',
                                        freqmin=cfs['F_high'],
                                        freqmax=cfs['F_low'],
                                        corners=cfs['order'],
                                        zerophase=cfs['zphase'])
                except ValueError as e:
                    IPUtils.errorPopup(str(e))

            else:
                IPUtils.errorPopup(filtType + ' filter not implemented yet')
                return
            filtered_stream += filtered_trace

        return filtered_stream

    def saveWindowGeometrySettings(self):
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('WaveformWidget')
        settings.setValue("main_splitterSettings", self.main_splitter.saveState())
        settings.setValue("rh_splitterSettings", self.rh_splitter.saveState())
        settings.setValue("lh_splitterSettings", self.lh_splitter.saveState())
        settings.setValue("plotviewer_splitterSettings", self.plotViewer.saveState())
        settings.endGroup()

    def restoreWindowGeometrySettings(self):
        # Restore settings
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('WaveformWidget')

        main_splitterSettings = settings.value("main_splitterSettings")
        if main_splitterSettings:
            self.main_splitter.restoreState(main_splitterSettings)

        rh_splitterSettings = settings.value("rh_splitterSettings")
        if rh_splitterSettings:
            self.rh_splitter.restoreState(rh_splitterSettings)

        lh_splitterSettings = settings.value("lh_splitterSettings")
        if lh_splitterSettings:
            self.lh_splitter.restoreState(lh_splitterSettings)

        pv_splitterSettings = settings.value("plotviewer_splitterSettings")
        if pv_splitterSettings:
            self.plotViewer.restoreState(pv_splitterSettings)
        else:
            pv_width = self.plotViewer.width()
            wsw = pv_width//6
            pww = pv_width - wsw
            self.plotViewer.setSizes([wsw, pww])

        settings.endGroup()

    @QtCore.pyqtSlot(str)
    def remove_trace(self, trace_id):

        for trace in self._sts.select(id=trace_id):
            self._sts.remove(trace)
            self.removeStation(trace.stats['network'], trace.stats['station'])

        self.statsViewer.setStats(self._sts)

        if len(self._sts) == 0:
            self._sts = None

        self.update_streams(self._sts)

    def removeStation(self, net_id, station_id):

        if self._inv is not None:
            try:
                self._inv = self.inv_remove(self._inv, network=net_id, station=station_id)

            except AttributeError as e:
                IPUtils.errorPopup(str(e))

            self.stationViewer.setInventory(self._inv)

    def inv_remove(self,
                   _inventory,
                   network='*',
                   station='*',
                   location='*',
                   channel='*',
                   keep_empty=False):

        selected = _inventory.select(network=network,
                                     station=station,
                                     location=location,
                                     channel=channel)

        selected_networks = [net for net in selected]
        selected_stations = [sta for net in selected_networks for sta in net]
        selected_channels = [cha for net in selected_networks
                             for sta in net for cha in sta]

        networks = []
        for net in _inventory:
            if net in selected_networks and station == '*' and \
                    location == '*' and channel == '*':
                continue
            stations = []
            for sta in net:
                if sta in selected_stations and location == '*' and channel == '*':
                    continue
                channels = []
                for cha in sta:
                    if cha in selected_channels:
                        continue
                    channels.append(cha)
                if not channels and not keep_empty:
                    continue
                sta = copy.copy(sta)
                sta.channels = channels
                stations.append(sta)

            if not stations and not keep_empty:
                continue
            net = copy.copy(net)
            net.stations = stations
            networks.append(net)

        return obspy.core.inventory.inventory.Inventory(networks, 'source')

    def clearWaveforms(self):
        # empty out the streams
        self._sts = None
        self._sts_filtered = None
        self.clear_inventory()

        # empty out the child widgets
        self.statsViewer.clear()
        self.plotViewer.clear()
        self.spectraWidget.clearPlot()

    @pyqtSlot(object)
    def update_signal_PSD(self, signal_region_item):

        if len(self._sts) == 0:
            self.spectraWidget.clearPlot()
            return

        signal_region = signal_region_item.getRegion()

        active_plot = self.plotViewer.pl_widget.get_active_plot()

        # calculate the PSD of the ---------------------------tart and finish
        dt = self._sts[active_plot].stats.delta
        start = int(signal_region[0] / dt)
        stop = int(signal_region[1] / dt)

        self.spectraWidget.updateSignalPSD(self._sts[active_plot][start:stop])

    @pyqtSlot(object)
    def update_noise_PSD(self, noise_region_item):

        if len(self._sts) == 0:
            self.spectraWidget.clearPlot()
            return

        noise_region = noise_region_item.getRegion()

        active_plot = self.plotViewer.pl_widget.get_active_plot()

        # calculate the PSD of the data in the current noise region
        dt = self._sts[active_plot].stats.delta
        start = int(noise_region[0] / dt)
        stop = int(noise_region[1] / dt)

        self.spectraWidget.updateNoisePSD(self._sts[active_plot][start:stop])

    @pyqtSlot(int, list, list, tuple, tuple)
    def update_widgets(self, index, lines, filtered_lines, signal_region, noise_region):
        # the -1 is sent if none of the plots are visible
        if len(self._sts) < 1 or index == -1:
            self.spectraWidget.set_title('...')
            self.spectraWidget.clearPlot()

        else:
            self.spectraWidget.set_title(self._sts[index].id)
            self.spectraWidget.set_fs(self._sts[index].stats.sampling_rate)

            noise_region_item = self.plotViewer.pl_widget.plot_list[index].getNoiseRegion()
            noise_region_item.sigRegionChanged.emit(noise_region_item)
            signal_region_item = self.plotViewer.pl_widget.plot_list[index].getSignalRegion()
            signal_region_item.sigRegionChanged.emit(signal_region_item)

            current_filter_display_settings = self.filterSettingsWidget.get_filter_display_settings()
            plot_title = self._sts[index].id
            if current_filter_display_settings['apply']:
                self.parent.beamformingWidget.setWaveform(filtered_lines[index], signal_region, plot_label=plot_title)
                self.parent.singleSensorWidget.setSignalWaveform(filtered_lines[index], signal_region, plot_label=plot_title)
                self.parent.singleSensorWidget.setNoiseWaveform(filtered_lines[index], noise_region, plot_label=plot_title)
            else:
                self.parent.beamformingWidget.setWaveform(lines[index], signal_region, plot_label=plot_title)
                self.parent.singleSensorWidget.setSignalWaveform(lines[index], signal_region, plot_label=plot_title)
                self.parent.singleSensorWidget.setNoiseWaveform(lines[index], noise_region, plot_label=plot_title)

            
