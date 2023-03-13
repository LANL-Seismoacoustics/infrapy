from PyQt5 import QtWidgets
from PyQt5.QtWidgets import (QWidget, QColorDialog, QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, 
                             QLineEdit, QToolBar, QToolButton, QVBoxLayout, QCheckBox, QComboBox, QLabel, QPushButton, QDoubleSpinBox)
from PyQt5.QtCore import QRect, QSize, Qt, pyqtSlot, pyqtSignal, QSettings

from PyQt5.QtGui import QPainter, QPaintEvent, QColor, QPalette

import matplotlib
from matplotlib import cm

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import urllib
import time
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from pyproj import Geod, transform

import numpy as np

# Make sure that we are using QT5
matplotlib.use('Qt5Agg')

from InfraView.widgets import IPUtils


class IPMapWidget(QWidget):

    fig = None
    axes = None
    transform = None
    projection = None
    detections = []
    resolution = ''
    extent = None

    current_linecolor='gray'

    gt_marker = None

    toolbar = None

    sta_lats = []
    sta_lons = []
    evt_lat = None
    evt_lon = None

    bisl_rslt = (None,None)  #(lat, lon)
    conf_ellipse = (None, None) #(dx, dy)

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

        self.map_settings_dialog = IPMapSettingsDialog()
        self.map_export_dialog = IPMapExportDialog(self, self.fig)
        self.missing_maps_dialog = IPMissingMapsDialog(self)

        self.extentWidget = IPExtentSettingsWidget(self)
        self.extentWidget.setVisible(False)

        self.toolbar = QToolBar()
        self.toolbar.setStyleSheet("QToolBar{background-color: silver; }")

        #self.toolbar.setStyleSheet("QToolBar { border-bottom: 1px solid; } ")
        self.tool_settings_button = QToolButton()
        self.tool_settings_button.setText("Settings...")

        self.tool_export_button = QToolButton()
        self.tool_export_button.setText("Export...")

        self.tool_extent_button = QToolButton()
        self.tool_extent_button.setText("Extent...")

        self.toolbar.addWidget(self.tool_settings_button)
        self.toolbar.addWidget(self.tool_extent_button)
        self.toolbar.addWidget(self.tool_export_button)

        main_layout = QVBoxLayout()
        main_layout.addWidget(self.toolbar)
        main_layout.addWidget(self.extentWidget)
        main_layout.addWidget(self.mapCanvas)

        self.setLayout(main_layout)

        self.compute_figure()
        self.draw_map()

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):

        self.map_settings_dialog.borders_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_dialog.states_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_dialog.lakes_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_dialog.rivers_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_dialog.coast_checkbox.stateChanged.connect(self.update_feature_visibilities)

        self.map_settings_dialog.signal_offline_directory_changed.connect(self.draw_map)

        self.map_settings_dialog.resolution_cb.currentTextChanged.connect(self.update_map)
        self.map_settings_dialog.signal_colors_changed.connect(self.update_map)
        self.map_settings_dialog.signal_background_changed.connect(self.update_map)
        self.map_settings_dialog.signal_map_settings_changed.connect(self.update_map)

        self.tool_settings_button.clicked.connect(self.map_settings_dialog.exec_)
        self.tool_export_button.clicked.connect(self.map_export_dialog.exec_)
        self.tool_extent_button.clicked.connect(self.showhide_extent_widget)

        self.extentWidget.hide_button.clicked.connect(self.hide_extent_widget)
        self.extentWidget.sig_extent_changed.connect(self.set_map_extent)
        self.extentWidget.sig_set_to_global.connect(self.set_map_extent_to_global)
        self.extentWidget.sig_autoscale.connect(self.autoscale_plot)

        # these technically aren't qt signals and slots, these are matplotlib callback connections
        #self.fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
        #self.fig.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.fig.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        # self.fig.canvas.mpl_connect('scroll_event', self.scroll_event_callback)

    def compute_figure(self):

        self.sph_proj = Geod(ellps='WGS84')
        self.projection = ccrs.PlateCarree()
        self.transform = ccrs.PlateCarree()

        self.axes = self.fig.add_subplot(1, 1, 1, projection=self.projection)

    @pyqtSlot()
    def draw_map(self, preserve_extent=False):
        if preserve_extent:
            current_extent = self.axes.get_extent()

        self.axes.clear()

        if self.map_settings_dialog.offline_checkbox.isChecked():
            # use the offline maps...
            cartopy.config['pre_existing_data_dir'] = self.map_settings_dialog.offline_directory_label.text()
        else:
            cartopy.config['pre_existing_data_dir'] = ""

        resolution = self.map_settings_dialog.resolution_cb.currentText()

        if self.map_settings_dialog.backgroud_image_checkbox.isChecked():
            land_facecolor = 'none'
            ocean_facecolor = 'none'
            self.axes.stock_img()
        else:
            land_facecolor = (self.map_settings_dialog.land_color_button.color().redF(), 
                            self.map_settings_dialog.land_color_button.color().greenF(),
                            self.map_settings_dialog.land_color_button.color().blueF())

            ocean_facecolor = (self.map_settings_dialog.ocean_color_button.color().redF(), 
                           self.map_settings_dialog.ocean_color_button.color().greenF(),
                           self.map_settings_dialog.ocean_color_button.color().blueF())

        land = cfeature.NaturalEarthFeature('physical',
                                            'land',
                                            scale=self.map_settings_dialog.resolution_cb.currentText(),
                                            edgecolor='face',
                                            facecolor=land_facecolor,
                                            linewidth=0.5)
        
        states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                        name='admin_1_states_provinces_lines',
                                                        scale=self.map_settings_dialog.resolution_cb.currentText(),
                                                        facecolor='none')

        self.land = self.axes.add_feature(land)

        self.oceans = self.axes.add_feature(cfeature.OCEAN.with_scale(resolution), facecolor=ocean_facecolor)

        self.states = self.axes.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
        self.lakes = self.axes.add_feature(cfeature.LAKES.with_scale(resolution))
        self.rivers = self.axes.add_feature(cfeature.RIVERS.with_scale(resolution))
        self.borders = self.axes.add_feature(cfeature.BORDERS.with_scale(resolution), linewidth=0.5)
        self.coast = self.axes.add_feature(cfeature.COASTLINE.with_scale(resolution), linewidth=0.5)

        #try:
        self.update_feature_visibilities()
        #except:
        #    IPUtils.errorPopup("There seems to be an issue with the map downloads. If you don't have access to the internet you can download the maps seperately, and use the offline maps setting in the Locations tab to point to the directory where they are downloaded to.")
        #    return

        if preserve_extent:
            self.set_map_extent(current_extent)
            self.extentWidget.set_extent_spin_values(current_extent)

    def showhide_extent_widget(self):
        if self.extentWidget.isVisible():
            self.extentWidget.setVisible(False)
        else:
            self.extentWidget.setVisible(True)

    def hide_extent_widget(self):
        self.extentWidget.setVisible(False)

    @pyqtSlot(list)
    def set_map_extent(self, extent):
        self.axes.set_extent(extent)
        self.fig.canvas.draw()  # update matlabplot

    @pyqtSlot()
    def set_map_extent_to_global(self):
        self.axes.set_global()
        global_extent = [-179.99, 180, -90, 90]
        self.extentWidget.set_extent_spin_values(global_extent)
        self.update_map()

    def update_feature_visibilities(self):
        try:
            # This shows/hides the various features shown on the map
            self.states.set_visible(self.map_settings_dialog.states_checkbox.isChecked())
            self.lakes.set_visible(self.map_settings_dialog.lakes_checkbox.isChecked())
            self.rivers.set_visible(self.map_settings_dialog.rivers_checkbox.isChecked())
            self.borders.set_visible(self.map_settings_dialog.borders_checkbox.isChecked())
            self.coast.set_visible(self.map_settings_dialog.coast_checkbox.isChecked())
            self.fig.canvas.draw()  # update matlabplot

        except urllib.error.URLError:

            return

    @pyqtSlot()
    def update_map(self, replot_bisl=True):
        self.draw_map(preserve_extent=True)
        self.update_detections(preserve_colors=True, autoscale=False)
        self.plot_ground_truth()
        if self.parent.dm_view.is_group_selected():
            self.plot_bisl_result(replot=replot_bisl)
            self.plot_conf_ellipse(replot=replot_bisl)
        self.draw_gridlines()
        self.fig.canvas.draw()  # update matlabplot

    def draw_gridlines(self, preserve_extent=True):
        
        self.gl = None
        if self.map_settings_dialog.show_grid_checkbox.isChecked():
            self.gl = self.axes.gridlines(draw_labels=True)
            

    def update_detections(self, line_color='gray', autoscale=True, preserve_colors=False):
        
        # trimmed_detections will hold either the entire set if not trimmed, or just the detections
        # chosen in the distance matrix.
        ip_detections = self.parent.get_trimmed_detections()
        
        if ip_detections is None:
            return

        if not autoscale:
            current_extent = self.axes.get_extent() # save this in case we are preserving the current extent
        
        if preserve_colors:
            linecolor = self.current_linecolor
        else:
            linecolor = line_color
            self.current_linecolor = linecolor

        self.clear_plot(reset_zoom=False)

        rng_max = self.parent.bislSettings.rng_max_edit.value() * 1000

        lons = []
        lats = []

        # for scaling purposes, lets keep a copy of the lons and lats in a seperate array
        for detection in ip_detections:
            lons.append(detection.longitude)
            lats.append(detection.latitude)

        # for scaling, lets keep track of the backaz line end points
        self.end_lats = []
        self.end_lons = []

        # this for loop draws the back azimuth lines. They will be length d (in degrees)
        for idx, detection in enumerate(ip_detections):

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

        for detection in ip_detections:
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

        if autoscale:
            self.autoscale_plot()
        else:
            # if we don't autoscale, then we want to return the plot to what it was when 
            # we entered this function
            self.set_map_extent(current_extent)
            self.extentWidget.set_extent_spin_values(current_extent)     # update extentWidget

        # draw it
        try:
            self.fig.canvas.draw()
        except:
            return

    @pyqtSlot()
    def clear_detections(self):
        #do i still need this?
        self.clear_plot()

    def plot_ground_truth(self):
        if self.gt_marker is not None: # clear old one, and make new one
            self.gt_marker.remove()

        lat = self.parent.showgroundtruth.event_widget.getLat()
        lon = self.parent.showgroundtruth.event_widget.getLon()
        
        current_extent = self.axes.get_extent() # plotting the event should not change the extent
        self.gt_marker, = self.axes.plot(lon, lat, 'X', color='red', transform=self.transform, markersize=16, gid='ground_truth_marker')

        self.set_map_extent(current_extent)
        self.extentWidget.set_extent_spin_values(current_extent)     # update extentWidget

        self.show_hide_ground_truth(self.parent.showgroundtruth.event_widget.showGT_cb.checkState())

    @pyqtSlot(int)
    def show_hide_ground_truth(self, show):

        if self.gt_marker is not None:
            if show == Qt.Checked:
                self.gt_marker.set_visible(True)
            else:
                self.gt_marker.set_visible(False)

        self.fig.canvas.draw()
        self.repaint()

    def plot_bisl_result(self, result_lon=None, result_lat=None, replot=False):
        # note that if replot is False, result_lat and result_lon are required

        #clear out previous marker
        self.remove_bisl_result()

        if replot: 
            # we just need to replot existing data
            result_lat = self.bisl_rslt[0]
            result_lon = self.bisl_rslt[1]
            if result_lat == None or result_lon == None:
                # nothing to plot
                return
        else:
            # we have a new result to plot
            self.bisl_rslt = (result_lat, result_lon)

        current_extent = self.axes.get_extent()
        self.axes.plot(result_lon, result_lat, 'o', markersize=7, color='blue', transform=self.transform, gid='bisl_result_marker')
        self.set_map_extent(current_extent)

        self.extentWidget.set_extent_spin_values(current_extent)     # update extentWidget

    def remove_bisl_result(self):
        for c in self.axes.get_children():
            if c.get_gid() == 'bisl_result_marker':
                c.remove()

    def plot_conf_ellipse(self, result_lons=None, result_lats=None, conf_dx=None, conf_dy=None, replot=False):
        # clear out previous ellipse
        self.remove_conf_ellipse()

        if replot:
            # we just need to replot existing data
            result_lats = self.bisl_rslt[0]
            result_lons = self.bisl_rslt[1]
            conf_dx = self.conf_ellipse[0]
            conf_dy = self.conf_ellipse[1]
            if result_lats == None or result_lons == None:
                # nothng to do
                return
        else:
            self.bisl_reslt = (result_lats, result_lons)
            self.conf_ellipse = (conf_dx, conf_dy)

        conf_lons, conf_lats = self.sph_proj.fwd(np.array([result_lons] * len(conf_dx)),
                                                 np.array([result_lats] * len(conf_dy)),
                                                 np.degrees(np.arctan2(conf_dx, conf_dy)),
                                                 np.sqrt(conf_dx**2 + conf_dy**2) * 1e3)[:2]

        current_extent = self.axes.get_extent()
        self.axes.plot(conf_lons, conf_lats, color='black', transform=self.transform, gid='conf_ellipse')
        self.set_map_extent(current_extent)
        self.extentWidget.set_extent_spin_values(current_extent)     # update extentWidget

        self.fig.canvas.draw()
        self.repaint()

    def remove_conf_ellipse(self):
        for c in self.axes.get_children():
            if c.get_gid() == 'conf_ellipse':
                c.remove()

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
            self.extentWidget.set_extent_spin_values([-179.99,180,-90,90])

        self.fig.canvas.draw()  # update matlabplot
        self.repaint()          # update widget

    #@pyqtSlot(float)
    def update_range_max(self, new_range):
        self.update_detections(preserve_colors=True)

    def autoscale_plot(self, source_location=None):
        # make an attempt to scale the plot so all relavent info is shown

        detections = self.parent.get_trimmed_detections()

        if len(detections) < 1:
            # nothing to scale to, so set to global extent and exit
            self.axes.set_global()
            self.extentWidget.set_extent_spin_values([-180,180,-90,90])
            return

        lons = []
        lats = []

        for detection in detections:
            lons.append(detection.longitude)
            lats.append(detection.latitude)

        """ if source_location is not None:
            lons.append(source_location[0])
            lats.append(source_location[1]) """

        if self.parent.showgroundtruth.event_widget.showGT_cb.isChecked():
            lons.append(self.parent.showgroundtruth.event_widget.event_lon_edit.value())
            lats.append(self.parent.showgroundtruth.event_widget.event_lat_edit.value())

        maxLat = max(lats + self.end_lats)
        minLat = min(lats + self.end_lats)

        maxLon = max(lons + self.end_lons)
        minLon = min(lons + self.end_lons)

        if maxLon != minLon:
            width = abs(maxLon - minLon)
        else:
            width = 20

        if maxLat != minLat:
            height = abs(maxLat - minLat)
        else:
            height = 20

        width_adj = width * 0.10
        height_adj = height * 0.10

        if maxLat == minLat and maxLon == minLon:
            # there is only one point, so behave accordingly
            new_extent = [minLon - width_adj, maxLon + width_adj, minLat - height_adj, maxLat + height_adj]
            self.set_map_extent(new_extent)
        else:
            # the normal case
            new_extent = [minLon - width_adj, maxLon + width_adj, minLat - height_adj, maxLat + height_adj]
            self.set_map_extent(new_extent)
        
        self.extentWidget.set_extent_spin_values(new_extent)     # update extentWidget

        # now redraw the gridlines since the extent has changed
        #self.update_map()

    def motion_notify_callback(self, event):
        if event.xdata is None or event.inaxes != self.axes:
            self.axes.set_title('')
            self.fig.canvas.draw()
            return

        elif event.button == 1:     # make sure the left button is clicked for a drag
            pass
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
            print("start = {}".format(self.start_mouse_loc))

    def button_release_callback(self, event):
        # This is to handle the button release from within matplotlib...undoes whatever the button press did
        if event.button != 1:
            return
        else:
            self.end_mouse_loc = [event.xdata, event.ydata]
            print("end = {}".format(self.end_mouse_loc))
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

        self.set_map_extent(extent)
        self.extentWidget.set_extent_spin_values(extent)

        self.fig.canvas.draw()

    

    def area_select_callback(self, eclick, erelease):
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
'''
'''
    @pyqtSlot(int)
    def set_central_longitude(self, cl):

        # if the function is called directly, not via a signal, then make sure the spinbox is correct
        if self.map_settings_dialog.central_lon_cb.value() != cl:
            self.map_settings_dialog.central_lon_cb.setValue(cl)

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


class IPMissingMapsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.buildUI()

    def buildUI(self):
        self.setWindowTitle("Infraview: Missing map files...")
        
        missing_maps_str = """When initially run, infrapy will attempt to download maps from 
        the internet. If you don't have an internet connection, it is most likely a proxy issue.
        In the rare case where you can't connect to the internet, then you can download the maps
        seperately and point infraview to that directory (see below)."""

        missing_maps_label = QLabel(missing_maps_str)
        map_dir_label = QLabel("Pre-existing maps directory: ")
        self.map_location_lineedit = QLineEdit()
        self.map_location_button = QPushButton("Browse...")
        self.map_location_button.clicked.connect(self.select_offline_maps_directory)

        select_map_dir_layout = QHBoxLayout()
        select_map_dir_layout.addWidget(map_dir_label)
        select_map_dir_layout.addWidget(self.map_location_lineedit)
        select_map_dir_layout.addWidget(self.map_location_button)

        ###   dialog buttons   ###
        buttons = QDialogButtonBox(QDialogButtonBox.Cancel,
                                   Qt.Horizontal,
                                   self)
        buttons.rejected.connect(self.reject)

        main_layout = QVBoxLayout()
        main_layout.addWidget(missing_maps_label)
        main_layout.addLayout(select_map_dir_layout)
        main_layout.addWidget(buttons)

        self.setLayout(main_layout)

    def select_offline_maps_directory(self):
        new_dir = QFileDialog.getExistingDirectory()
        
        self.map_location_lineedit.setText(new_dir) 
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('LocationWidget')
        settings.setValue('offline_maps_dir', new_dir)            


class IPMapExportDialog(QDialog):
    
    def __init__(self, parent=None, figure=None):
        super().__init__(parent)
        self.fig = figure
        self.buildUI()

    def buildUI(self):
        self.setWindowTitle("Infraview: Map Export")

        # export pdf...
        pdf_group_box = QGroupBox("Export to PDF")
        self.pdf_file_label = QLabel("")
        self.pdf_file_label.setMinimumWidth(200)
        self.pdf_button = QPushButton("Choose file..")
        self.pdf_export_button = QPushButton("Export")
        pdf_layout = QHBoxLayout()
        pdf_layout.addWidget(self.pdf_button)
        pdf_layout.addWidget(self.pdf_file_label)
        pdf_layout.addWidget(self.pdf_export_button)
        pdf_group_box.setLayout(pdf_layout)

        #export img...
        img_group_box = QGroupBox("Export to image file")
        self.img_file_label = QLabel("")
        self.img_file_label.setMinimumWidth(200)
        self.img_button = QPushButton("Choose file...")
        self.img_export_button = QPushButton("Export")
        img_layout = QHBoxLayout()
        img_layout.addWidget(self.img_button)
        img_layout.addWidget(self.img_file_label)
        img_layout.addWidget(self.img_export_button)
        img_group_box.setLayout(img_layout)

        ###   dialog buttons   ###
        buttons = QDialogButtonBox(QDialogButtonBox.Cancel,
                                   Qt.Horizontal,
                                   self)
        buttons.rejected.connect(self.reject)

        main_layout = QVBoxLayout()
        main_layout.addWidget(pdf_group_box)
        main_layout.addWidget(img_group_box)
        main_layout.addStretch()
        main_layout.addWidget(buttons)

        self.setLayout(main_layout)

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.pdf_button.clicked.connect(self.save_pdf)
        self.img_button.clicked.connect(self.save_img)
        self.img_export_button.clicked.connect(self.export_img)
        self.pdf_export_button.clicked.connect(self.export_pdf)


    def save_pdf(self):
        filename = QFileDialog.getSaveFileName(parent=self, caption="Save PDF", filter="PDF files (*.pdf)" )
        if filename[0] == '':
            # dialog was cancelled, just leave
            return

        if filename[0].endswith('.pdf'):
            new_filename = filename[0]
        else:
            if filename[0] != "":
                new_filename = filename[0] + '.pdf'

        self.pdf_file_label.setText(new_filename)

    def save_img(self):
        filename = QFileDialog.getSaveFileName(parent=self, caption="Save Image", filter="Images (*.png *.xpm *.jpg)" )
        if filename[0] == '':
            # dialog was cancelled, just leave
            return

        self.img_file_label.setText(filename[0])

    def export_img(self):
        if self.img_file_label.text() == "":
            IPUtils.errorPopup("Can't export image.  No image file selected.")
            return
        self.fig.savefig(self.img_file_label.text())
        time.sleep(0.5)
        self.close()

    def export_pdf(self):
        if self.pdf_file_label.text() == "":
            IPUtils.errorPopup("Can't export to pdf.  No pdf file selected.")
            return 
        self.fig.savefig(self.pdf_file_label.text())
        time.sleep(1.2)
        self.close()


class IPMapSettingsDialog(QDialog):

    signal_colors_changed = pyqtSignal()
    signal_background_changed = pyqtSignal()
    signal_offline_directory_changed = pyqtSignal()
    signal_map_settings_changed = pyqtSignal()
    
    ocean_color = QColor(0, 107, 166)
    land_color = QColor(222, 222, 222)
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.setMinimumWidth(400)
        self.setWindowTitle("InfraView: Map Settings")
        self.buildUI()

    def buildUI(self):

        ###   feature settings   ###
        features_gb = QGroupBox("Features")
        self.borders_checkbox = QCheckBox('Countries  ')
        self.states_checkbox = QCheckBox('States and Provinces  ')
        self.lakes_checkbox = QCheckBox('Lakes  ')
        self.rivers_checkbox = QCheckBox('Rivers  ')
        self.coast_checkbox = QCheckBox('Coastline  ')

        features_layout = QVBoxLayout()
        features_layout.addWidget(self.borders_checkbox)
        features_layout.addWidget(self.states_checkbox)
        features_layout.addWidget(self.lakes_checkbox)
        features_layout.addWidget(self.rivers_checkbox)
        features_layout.addWidget(self.coast_checkbox)
        
        features_gb.setLayout(features_layout)

        ###   color settings   ###
        colors_gb = QGroupBox("Colors")
        ocean_color_label = QLabel("Oceans: ")
        self.ocean_color_button = IPColorButton(self.ocean_color)
        land_color_label = QLabel("Land: ")
        self.land_color_button = IPColorButton(self.land_color)

        colors_layout = QFormLayout()
        colors_layout.addRow(ocean_color_label, self.ocean_color_button)
        colors_layout.addRow(land_color_label, self.land_color_button)

        colors_gb.setLayout(colors_layout)

        ###   grid settings
        # grid_gb = QGroupBox("Grid Lines")
        self.show_grid_checkbox = QCheckBox("Show Grid Lines ")

        ###   resolution settings   ###
        label_resolution = QLabel(self.tr('Resolution'))
        self.resolution_cb = QComboBox()
        self.resolution_cb.addItem('50m')
        self.resolution_cb.addItem('110m')
        self.resolution_cb.setCurrentIndex(1)

        resolution_layout = QHBoxLayout()
        resolution_layout.addWidget(label_resolution)
        resolution_layout.addWidget(self.resolution_cb)

        ### background image ##
        self.backgroud_image_checkbox = QCheckBox('Use background image  ')

        ###   offline maps settings   ###
        self.offline_checkbox = QCheckBox('Use offline maps  ')
        self.offline_directory_label = QLabel("Use offline maps")
        # read in the offline_director from settings if there is one
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('LocationWidget')
        odd = settings.value('offline_maps_dir', '')
        odd_isChecked_str = settings.value('use_offline_cb', 'False')
        if type(odd_isChecked_str) is str:
            odd_isChecked = odd_isChecked_str.lower() == 'true'
        else:
            odd_isChecked = odd_isChecked_str
        settings.endGroup()

        self.offline_directory_label.setText(odd)
        # for now, if there is a directory in the offline_directory_label, assume they want to use that, and activate checkbox
        self.offline_checkbox.setChecked(odd_isChecked)
        self.offline_directory_label.setEnabled(odd_isChecked)

        self.offline_directory_select_button = QPushButton("Select Folder...")
        self.offline_directory_select_button.setEnabled(odd_isChecked)

        self.offline_file_dialog = QFileDialog()
        self.offline_file_dialog.setFileMode(QFileDialog.Directory)

        offline_layout = QHBoxLayout()
        offline_layout.addWidget(self.offline_checkbox)
        offline_layout.addWidget(self.offline_directory_label)
        offline_layout.addWidget(self.offline_directory_select_button)

        ###   dialog buttons   ###
        buttons = QDialogButtonBox(QDialogButtonBox.Ok,
                                   Qt.Horizontal,
                                   self)
        buttons.accepted.connect(self.accept)

        ###   layouts   ###
        boxes_layout = QHBoxLayout()
        boxes_layout.addWidget(features_gb)
        boxes_layout.addWidget(colors_gb)

        main_layout = QVBoxLayout()
        main_layout.addLayout(boxes_layout)
        main_layout.addWidget(self.backgroud_image_checkbox)
        #main_layout.addWidget(self.show_grid_checkbox)
        main_layout.addLayout(resolution_layout)
        main_layout.addLayout(offline_layout)
        main_layout.addStretch()
        main_layout.addWidget(buttons)

        self.setLayout(main_layout)

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.offline_checkbox.clicked.connect(self.offline_directory_select_button.setEnabled)
        self.offline_checkbox.clicked.connect(self.offline_directory_label.setEnabled)
        self.offline_checkbox.clicked.connect(self.update_settings)

        self.ocean_color_button.clicked.connect(self.update_ocean_color)
        self.land_color_button.clicked.connect(self.update_land_color)

        self.backgroud_image_checkbox.clicked.connect(self.toggle_background_image)
        self.show_grid_checkbox.clicked.connect(self.update_grid_lines)
        self.offline_directory_select_button.clicked.connect(self.select_offline_maps_directory)

    def toggle_background_image(self):
        self.land_color_button.setDisabled(self.backgroud_image_checkbox.isChecked())
        self.ocean_color_button.setDisabled(self.backgroud_image_checkbox.isChecked())
        self.signal_background_changed.emit()

    def update_grid_lines(self):
        self.signal_map_settings_changed.emit()

    def update_ocean_color(self):
        new_color = QColorDialog.getColor(self.ocean_color_button.color())
        if new_color.isValid():
            self.ocean_color_button.set_color(new_color)
            self.signal_colors_changed.emit()

    def update_land_color(self):
        new_color = QColorDialog.getColor(self.land_color_button.color())
        if new_color.isValid(): 
            self.land_color_button.set_color(new_color)
            self.signal_colors_changed.emit()

    def select_offline_maps_directory(self):
        curr_dir = self.offline_directory_label.text()
        
        new_dir = QFileDialog.getExistingDirectory()
        
        self.offline_directory_label.setText(new_dir) 
        self.signal_offline_directory_changed.emit()
        
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('LocationWidget')
        settings.setValue('offline_maps_dir', new_dir)
        settings.endGroup()

    def update_settings(self):
        settings = QSettings('LANL', 'InfraView')
        settings.beginGroup('LocationWidget')
        settings.setValue('use_offline_cb', self.offline_checkbox.isChecked())
        settings.endGroup()


class IPColorButton(QPushButton):
    current_color = QColor(255, 0, 0)

    def __init__(self, color):
        super().__init__()
        self.current_color = color
        size = QSize(self.height(), self.height())
        self.setFixedSize(QSize(26,26))

    def paintEvent(self, a0: QPaintEvent) -> None:
        super().paintEvent(a0)
        r = QRect(0, 0, self.width() * 0.75, self.height() * 0.75)
        r.moveTo(self.rect().center() - r.center())
        painter = QPainter(self)
        painter.setBrush(self.current_color)
        painter.drawRect(r)

    def set_color(self, new_color):
        # new color should be a QColor type
        self.current_color = new_color

    def color(self):
        return QColor(self.current_color)


class IPExtentSettingsWidget(QWidget):

    sig_extent_changed = pyqtSignal(list)
    sig_set_to_global = pyqtSignal()
    sig_autoscale = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.buildUI()

    def buildUI(self):

        # pal = QPalette()
        # pal.setColor(QPalette.Window, Qt.lightGray)
        # self.setAutoFillBackground(True)
        # self.setPalette(pal)

        ll_groupbox = QGroupBox("Lower left coordinates")
        ur_groupbox = QGroupBox("Upper right coordinates")

        self.ll_lat_spin = QDoubleSpinBox()
        self.ll_lat_spin.setMaximumWidth(100)
        self.ll_lat_spin.setRange(-90.0, 90.0)
        self.ll_lat_spin.setValue(-90.0)
        self.ll_lat_spin.setPrefix("Lat: ")
        self.ll_lat_spin.valueChanged.connect(self.activate_update_button)

        self.ll_lon_spin = QDoubleSpinBox()
        self.ll_lon_spin.setMaximumWidth(100)
        self.ll_lon_spin.setRange(-179.99, 180.0)
        self.ll_lon_spin.setValue(-179.99)
        self.ll_lon_spin.setPrefix("Lon: ")
        self.ll_lon_spin.valueChanged.connect(self.activate_update_button)

        ll_layout = QHBoxLayout()
        ll_layout.addWidget(self.ll_lon_spin)
        ll_layout.addWidget(self.ll_lat_spin)
        ll_groupbox.setLayout(ll_layout)

        self.ur_lat_spin = QDoubleSpinBox()
        self.ur_lat_spin.setMaximumWidth(100)
        self.ur_lat_spin.setRange(-90., 90.0)
        self.ur_lat_spin.setValue(90.0)
        self.ur_lat_spin.setPrefix("Lat: ")
        self.ur_lat_spin.valueChanged.connect(self.activate_update_button)

        self.ur_lon_spin = QDoubleSpinBox()
        self.ur_lon_spin.setMaximumWidth(100)
        self.ur_lon_spin.setRange(-179.99, 180.0)
        self.ur_lon_spin.setValue(180.0)
        self.ur_lon_spin.setPrefix("Lon: ")
        self.ur_lon_spin.valueChanged.connect(self.activate_update_button)

        ur_layout = QHBoxLayout()
        ur_layout.addWidget(self.ur_lon_spin)
        ur_layout.addWidget(self.ur_lat_spin)
        ur_groupbox.setLayout(ur_layout)

        self.update_plot_button = QPushButton("Update")
        self.update_plot_button.setMaximumWidth(100)
        self.update_plot_button.setEnabled(False)
        self.update_plot_button.clicked.connect(self.deactivate_update_button)
        self.update_plot_button.clicked.connect(self.update_map_extent)

        self.set_to_global_button = QPushButton("Global")
        self.set_to_global_button.setMaximumWidth(100)
        self.set_to_global_button.clicked.connect(self.set_to_global)

        self.autoscale_button = QPushButton("Autoscale")
        self.autoscale_button.setMaximumWidth(100)
        self.autoscale_button.clicked.connect(self.autoscale_map)

        self.hide_button = QPushButton("Hide")
        self.hide_button.setMaximumWidth(60)

        h_layout = QHBoxLayout()
        h_layout.addWidget(ll_groupbox)
        h_layout.addWidget(ur_groupbox)
        h_layout.addWidget(self.update_plot_button)
        h_layout.addWidget(self.set_to_global_button)
        h_layout.addWidget(self.autoscale_button)
        h_layout.addStretch()
        h_layout.addWidget(self.hide_button)
        h_layout.setContentsMargins(0,0,0,0)
        self.setLayout(h_layout) 

    def set_extent_spin_values(self, extent):
        # ll_lon: lower left longitude
        # ur_lat: upper right latitude
        # etc

        self.ll_lon_spin.setValue(extent[0])
        self.ll_lat_spin.setValue(extent[2])
        self.ur_lon_spin.setValue(extent[1])
        self.ur_lat_spin.setValue(extent[3])

    def set_to_global(self):
        self.sig_set_to_global.emit()

    def autoscale_map(self):
        self.sig_autoscale.emit()

    def activate_update_button(self):
        self.update_plot_button.setEnabled(True)

    def deactivate_update_button(self):
        self.update_plot_button.setEnabled(False)

    def update_map_extent(self):
        extent = [self.ll_lon_spin.value(), self.ur_lon_spin.value(), self.ll_lat_spin.value(), self.ur_lat_spin.value()]
        self.sig_extent_changed.emit(extent)
