from PyQt5 import QtWidgets
from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QVBoxLayout, QAction, QCheckBox, QComboBox, QLabel, QSplitter, QSpinBox)
from PyQt5.QtCore import Qt, pyqtSlot, pyqtSignal

import matplotlib
from matplotlib import cm

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import urllib
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl as cmpl
import cartopy.mpl.geoaxes as cgeoaxes

from pyproj import Geod, transform

from obspy.core.inventory import Inventory, Network, Station, Channel, Site

import numpy as np

# Make sure that we are using QT5
matplotlib.use('Qt5Agg')


class IPMapWidget(QWidget):

    fig = None
    axes = None
    transform = None
    projection = None
    detections = []
    resolution = ''
    extent = None

    gt_marker = None

    toolbar = None

    sta_lats = []
    sta_lons = []
    evt_lat = None
    evt_lon = None

    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.buildUI()

    def buildUI(self):

        self.fig = Figure()
        self.zoom = 1
        self.mapCanvas = FigureCanvas(self.fig)
        self.mapCanvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                     QtWidgets.QSizePolicy.Expanding)

        main_splitter = QSplitter(Qt.Vertical)

        self.map_settings_widget = IPMapSettingsWidget()

        main_splitter.addWidget(self.map_settings_widget)
        main_splitter.addWidget(self.mapCanvas)

        main_splitter.setStretchFactor(0, 0)
        main_splitter.setStretchFactor(1, 1)

        # I want to make sure the splitter handle is easily seen...
        main_splitter.setStyleSheet("QSplitter::handle{ background-color: #444444}")
        main_splitter.setHandleWidth(2)

        main_layout = QVBoxLayout()
        main_layout.addWidget(main_splitter)

        self.setLayout(main_layout)

        self.compute_figure()

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):

        self.map_settings_widget.borders_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_widget.states_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_widget.lakes_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_widget.rivers_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_widget.coast_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_widget.resolution_cb.currentTextChanged.connect(self.update_resolution)

        # these technically aren't qt signals and slots, these are matplotlib callback connections
        self.fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.fig.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.fig.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        # self.fig.canvas.mpl_connect('scroll_event', self.scroll_event_callback)

    def compute_figure(self):

        self.sph_proj = Geod(ellps='WGS84')
        self.projection = ccrs.PlateCarree()
        self.transform = ccrs.PlateCarree()

        self.axes = self.fig.add_subplot(1, 1, 1, projection=self.projection)

        self.draw_map()

    def draw_map(self):
        self.axes.clear()

        resolution = self.map_settings_widget.resolution_cb.currentText()

        land = cfeature.NaturalEarthFeature('physical',
                                            'land',
                                            scale=self.map_settings_widget.resolution_cb.currentText(),
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'],
                                            linewidth=0.5)

        states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                        name='admin_1_states_provinces_lines',
                                                        scale=self.map_settings_widget.resolution_cb.currentText(),
                                                        facecolor='none')

        self.land = self.axes.add_feature(land)
        self.oceans = self.axes.add_feature(cfeature.OCEAN.with_scale(resolution), facecolor=(22. / 255., 43. / 255., 72. / 255., 0.5))

        self.states = self.axes.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
        self.lakes = self.axes.add_feature(cfeature.LAKES.with_scale(resolution))
        self.rivers = self.axes.add_feature(cfeature.RIVERS.with_scale(resolution))
        self.borders = self.axes.add_feature(cfeature.BORDERS.with_scale(resolution), linewidth=0.5)
        # cartopy stopped sending coasts?
        #self.coast = self.axes.add_feature(cfeature.COASTLINE.with_scale(resolution))   

        self.update_feature_visibilities()

    def update_feature_visibilities(self):
        # This shows/hides the various features shown on the map
        self.states.set_visible(self.map_settings_widget.states_checkbox.isChecked())
        self.lakes.set_visible(self.map_settings_widget.lakes_checkbox.isChecked())
        self.rivers.set_visible(self.map_settings_widget.rivers_checkbox.isChecked())
        self.borders.set_visible(self.map_settings_widget.borders_checkbox.isChecked())
        #self.coast.set_visible(self.map_settings_widget.coast_checkbox.isChecked())
        try:
            self.fig.canvas.draw()  # update matlabplot
        except urllib.error.URLError:
            print('problem with download...')



    def update_detections(self, ip_detections=None, linecolor='gray', autoscale=True):
        self.clear_plot(reset_zoom=False)

        if ip_detections is not None:
            self.detections = ip_detections
        
        if self.detections is None:
            return

        current_extent = self.axes.get_extent()
        rng_max = self.parent.bislSettings.rng_max_edit.value() * 1000

        lons = []
        lats = []

        # for scaling purposes, lets keep a copy of the lons and lats in a seperate array
        for detection in self.detections:
            lons.append(detection.longitude)
            lats.append(detection.latitude)

        # for scaling, lets keep track of the backaz line end points
        self.end_lats = []
        self.end_lons = []

        # this for loop draws the back azimuth lines. They will be length d (in degrees)
        for idx, detection in enumerate(self.detections):

            p_lons = [detection.longitude]
            p_lats = [detection.latitude]
            count = 0

            if hasattr(detection, 'index'):
                name = detection.index
            else:
                name = str(idx)

            # draw the back azimuth lines    
            N = 20
            # a hack... annotate doesnt correctly ingest the transform, so you have to do this...which i don't entirely understand
            # https://stackoverflow.com/questions/25416600/why-the-annotate-worked-unexpected-here-in-cartopy#_=_
            mpl_transform = ccrs.PlateCarree()._as_mpl_transform(self.axes)

            for d in np.arange(0, rng_max, rng_max / N):  # N points
                new_lon, new_lat, _ = self.sph_proj.fwd(detection.longitude, detection.latitude, detection.back_azimuth, d)
                if count == int(N / 2):
                    self.axes.annotate(name,
                                       (new_lon, new_lat),
                                       textcoords='offset points',
                                       xytext=(0, 10),
                                       xycoords=mpl_transform,
                                       ha='center',
                                       gid='detection_label')
                p_lons.append(new_lon)
                p_lats.append(new_lat)
                count += 1

            self.end_lats.append(p_lats[-1])
            self.end_lons.append(p_lons[-1])

            self.axes.plot(p_lons,
                           p_lats,
                           color=linecolor,
                           transform=self.transform,
                           gid='detection_line')
            self.axes.set_extent(current_extent)

        for detection in self.detections:
            if detection.array_dim == 3:
                symbol = '^'                    # triangle
            elif detection.array_dim == 4:
                symbol = 's'                    # square
            elif detection.array_dim == 5:
                symbol = 'p'                    # pentagon
            elif detection.array_dim == 6:
                symbol = 'H'                    # hexagon
            else:
                symbol = 'o'                    # circle

            self.axes.plot(detection.longitude,
                           detection.latitude,
                           marker=symbol,
                           markersize=7,
                           color='black',
                           transform=self.transform,
                           gid='detection_marker')
            self.axes.set_extent(current_extent)

        # the plot function will automatically zoom into the view, but we want more control than that, so return the extent
        # to what is was before plotting
        self.axes.set_extent(current_extent, crs=self.transform)

        if autoscale:
            self.autoscale_plot()

        # draw it
        try:
            self.fig.canvas.draw()
            self.repaint()
        except:
            return
        

    @pyqtSlot(str)
    def update_resolution(self, new_resolution):
        # we want to preserve the line color when we replot, so lets try to figure out what that is first
        line_color = None
        for c in self.axes.get_children():
            if c.get_gid() == 'detection_line':
                line_color = c.get_color()

        self.draw_map()
        self.update_detections(self.detections, linecolor=line_color)

    @pyqtSlot()
    def clear_detections(self):
        self.detections = []
        self.clear_plot()

    @pyqtSlot(float, float)
    def plot_ground_truth(self, lon, lat):

        if self.gt_marker is not None:
            self.gt_marker.remove()
        current_extent = self.axes.get_extent()
        self.gt_marker, = self.axes.plot(lon, lat, 'X', color='red', transform=self.transform, markersize=16, gid='ground_truth_marker')
        self.axes.set_extent(current_extent)
        self.fig.canvas.draw()
        self.show_hide_ground_truth(self.parent.showgroundtruth.show_gt())
        self.repaint()

    @pyqtSlot(bool)
    def show_hide_ground_truth(self, show):

        if self.gt_marker is not None:
            self.gt_marker.set_visible(show)
        self.fig.canvas.draw()
        self.repaint()

    def plot_bisl_result(self, result_lon, result_lat):
        for c in self.axes.get_children():
            if c.get_gid() == 'bisl_result_marker':
                c.remove()

        current_extent = self.axes.get_extent()
        self.axes.plot(result_lon, result_lat, 'o', markersize=7, color='blue', transform=self.transform, gid='bisl_result_marker')
        self.axes.set_extent(current_extent)

    def plot_conf_ellipse(self, result_lons, result_lats, conf_dx, conf_dy):
        for c in self.axes.get_children():
            if c.get_gid() == 'conf_ellipse':
                c.remove()

        conf_lons, conf_lats = self.sph_proj.fwd(np.array([result_lons] * len(conf_dx)),
                                                 np.array([result_lats] * len(conf_dy)),
                                                 np.degrees(np.arctan2(conf_dx, conf_dy)),
                                                 np.sqrt(conf_dx**2 + conf_dy**2) * 1e3)[:2]

        current_extent = self.axes.get_extent()
        self.axes.plot(conf_lons, conf_lats, color='black', transform=self.transform, gid='conf_ellipse')
        self.axes.set_extent(current_extent)

        self.fig.canvas.draw()
        self.repaint()

    def clear_plot(self, reset_zoom=True):
        for c in self.axes.get_children():
            c_gid = c.get_gid()
            if c_gid == 'detection_label':
                c.remove()
            elif c_gid == 'detection_marker':
                c.remove()
            elif c_gid == 'detection_line':
                c.remove()
            elif c_gid == 'bisl_result_marker':
                c.remove()
            elif c_gid == 'conf_ellipse':
                c.remove()
        if reset_zoom:
            self.axes.set_global()

        self.fig.canvas.draw()  # update matlabplot
        self.repaint()          # update widget

    @pyqtSlot(float)
    def update_range_max(self, new_range):
        # we want to preserve the line color when we replot, so lets try to figure out what that is
        # first
        line_color = 'white'
        for c in self.axes.get_children():
            if c.get_gid() == 'detection_line':
                line_color = c.get_color()
                break

        self.update_detections(self.detections, linecolor=line_color)

    

    def autoscale_plot(self, source_location=None):
        # make an attempt to scale the plot so all relavent info is shown

        if len(self.detections) < 1:
            # nothing to scale to, so set to global extent and exit
            self.axes.set_global()
            return

        lons = []
        lats = []

        for detection in self.detections:
            lons.append(detection.longitude)
            lats.append(detection.latitude)

        if source_location is not None:
            lons.append(source_location[0])
            lats.append(source_location[1])

        if self.parent.showgroundtruth.showGT_cb.isChecked():
            lons.append(self.parent.parent.eventWidget.event_lon_edit.value())
            lats.append(self.parent.parent.eventWidget.event_lat_edit.value())

        maxLat = max(lats + self.end_lats)
        minLat = min(lats + self.end_lats)

        maxLon = max(lons + self.end_lons)
        minLon = min(lons + self.end_lons)

        if maxLon != minLon:
            width = abs(maxLon - minLon)
            if width < 250:
                width = 250
        else:
            width = 250

        if maxLat != minLat:
            height = abs(maxLat - minLat)
            if height < 250:
                height = 250
        else:
            height = 250

        if maxLat == minLat and maxLon == minLon:
            # there is only one point, so behave accordingly
            self.axes.set_extent([minLon - width / 10.,
                                 maxLon + width / 10.,
                                 minLat - height / 10.,
                                 maxLat + height / 10.],
                                 crs=self.transform)
        else:
            # the normal case
            self.axes.set_extent([minLon - width / 10.,
                                 maxLon + width / 10.,
                                 minLat - height / 10.,
                                 maxLat + height / 10.],
                                 crs=self.transform)

    def motion_notify_callback(self, event):
        if event.xdata is None or event.inaxes != self.axes:
            self.axes.set_title('')
            self.fig.canvas.draw()
            return

        elif event.button == 1:     # make sure the left button is clicked for a drag
            extent = self.axes.get_extent()

            dx = event.xdata - self.start_mouse_loc[0]
            dy = event.ydata - self.start_mouse_loc[1]

            lo1 = extent[0] - dx
            lo2 = extent[1] - dx

            la1 = extent[2] - dy
            la2 = extent[3] - dy

            if lo1 < -180 or lo2 > 180:
                return

            if la1 < -90 or la2 > 90:
                return

            self.axes.set_extent([lo1, lo2, la1, la2], crs=self.transform)
            self.fig.canvas.draw_idle()
        else:
            self.mouse_moved = True
            # if event.button is None:
            self.axes.set_title('Lon = {:+f}, Lat = {:+f}'.format(event.xdata, event.ydata), loc='center', pad=20, fontsize=10)
            self.fig.canvas.draw()

    # matplotlib events are not to be confused with (py)Qt events
    def button_press_callback(self, event):
        # This is to handle the button click from within matplotlib...doesnt really do anything yet
        if event.button != 1:
            return
        else:
            self.start_mouse_loc = [event.xdata, event.ydata]
            # print("start = {}".format(self.start_mouse_loc))

    def button_release_callback(self, event):
        # This is to handle the button release from within matplotlib...undoes whatever the button press did
        if event.button != 1:
            return
        else:
            self.end_mouse_loc = [event.xdata, event.ydata]
            # print("end = {}".format(self.end_mouse_loc))
'''
    # Matplotlib callbacks go here_____________________

    def scroll_event_callback(self, event):
        # ZOOOOOOOOM
        extent = self.axes.get_extent()

        if event.button == 'down':
            zoom = 0.25
        elif event.button == 'up':
            zoom = -0.25

        # if the mouse position is unmoving, continue to zoom in on the initial position,
        # if the mouse moved, update the center of the zoom to the new position
        if self.mouse_moved:
            mousex = event.xdata
            mousey = event.ydata
            self.startx = mousex
            self.starty = mousey
        else:
            mousex = self.startx
            mousey = self.starty

        self.mouse_moved = False    # set it false, if the mouse is moved this will flip to True

        new_width = abs(extent[1] - extent[0]) * (1 + zoom)
        new_height = abs(extent[3] - extent[2]) * (1 + zoom)

        lo1 = mousex - new_width / 2.
        lo1 = lo1 if lo1 >= -180 else -180
        lo2 = mousex + new_width / 2.
        lo2 = lo2 if lo2 <= 180 else 180

        la1 = mousey - new_height / 2.
        la1 = la1 if la1 >= -90 else -90
        la2 = mousey + new_height / 2.
        la2 = la2 if la2 <= 90 else 90

        extent = [lo1, lo2, la1, la2]

        self.axes.set_extent(extent, crs=self.transform)

        self.fig.canvas.draw()

    

    def area_select_callback(self, eclick, erelease):
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
'''
'''
    @pyqtSlot(int)
    def set_central_longitude(self, cl):

        # if the function is called directly, not via a signal, then make sure the spinbox is correct
        if self.map_settings_widget.central_lon_cb.value() != cl:
            self.map_settings_widget.central_lon_cb.setValue(cl)

        # in order to preserve the zoom, lets first grab the current extent with the goal of resetting it after the new projection
        current_extent = self.axes.get_extent()
        width = current_extent[1] - current_extent[0]
        height = current_extent[3] - current_extent[2]

        self.projection = ccrs.PlateCarree(central_longitude=cl)
        self.transform = ccrs.PlateCarree()

        self.axes.remove()
        self.axes = self.fig.add_subplot(1, 1, 1, projection=self.projection)
        self.draw_map()

        self.update_detections(autoscale=False)
'''

    


class IPMapSettingsWidget(QWidget):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.buildUI()

    def buildUI(self):

        self.borders_checkbox = QCheckBox('Countries  ')
        self.states_checkbox = QCheckBox('States and Provinces  ')
        self.lakes_checkbox = QCheckBox('Lakes  ')
        self.rivers_checkbox = QCheckBox('Rivers  ')
        self.coast_checkbox = QCheckBox('Coastline  ')

        label_resolution = QLabel(self.tr('Resolution'))
        self.resolution_cb = QComboBox()
        self.resolution_cb.addItem('50m')
        self.resolution_cb.addItem('110m')
        self.resolution_cb.setCurrentIndex(1)

        hboxlayout = QHBoxLayout()
        hboxlayout.addWidget(self.borders_checkbox)
        hboxlayout.addWidget(self.states_checkbox)
        hboxlayout.addWidget(self.lakes_checkbox)
        hboxlayout.addWidget(self.rivers_checkbox)
        # hboxlayout.addWidget(self.coast_checkbox)
        hboxlayout.addStretch()
        hboxlayout.addWidget(label_resolution)
        hboxlayout.addWidget(self.resolution_cb)

        hboxlayout.setContentsMargins(0, 0, 0, 2)

        self.setLayout(hboxlayout)
