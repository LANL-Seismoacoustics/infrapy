from PyQt5 import QtWidgets
from PyQt5.QtWidgets import (QWidget, QColorDialog, QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, 
                             QLineEdit, QVBoxLayout, QCheckBox, QComboBox, QLabel, QPushButton, QDoubleSpinBox)
from PyQt5.QtCore import QRect, QSize, Qt, pyqtSlot, pyqtSignal, QSettings
from PyQt5.QtGui import QPainter, QPaintEvent, QColor, QPalette


from InfraView.widgets import IPBaseWidgets
from InfraView.widgets import IPUtils

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import urllib
import time
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from pyproj import Geod

import numpy as np

# Make sure that we are using QT5
matplotlib.use('Qt5Agg')

from InfraView.widgets import IPUtils
from InfraView.widgets import IPBaseWidgets


class IPMapWidget(QWidget):
    def __init__(self, parent):
        super().__init__()

        self.fig = None
        self.axes = None
        self.transform = None
        self.projection = None
        self.detections = []
        self.resolution = ''
        self.extent = None
  
        self.current_linecolor='gray'

        self.gt_marker = None

        self.toolbar = None

        self.sta_lats = []
        self.sta_lons = []
        self.evt_lat = None
        self.evt_lon = None

        self.bisl_rslt = (None,None)  #(lat, lon)
        self.conf_ellipse = (None, None) #(dx, dy)

        
        self.parent = parent
        self.buildUI()

    def buildUI(self):

        self.fig = Figure()
        self.zoom = 1
        self.mapCanvas = FigureCanvas(self.fig)

        self.map_export_dialog = IPMapExportDialog(self, self.fig)
        self.missing_maps_dialog = IPMissingMapsDialog(self)

        main_layout = QVBoxLayout()
        main_layout.addWidget(self.mapCanvas)

        self.setLayout(main_layout)

        self.compute_figure()

    def update_theme(self, t):
        if t == 'light':
            self.fig.patch.set_facecolor('w')
        elif t == 'dark':
            self.fig.patch.set_facecolor(IPUtils.ip_dark_grey_hex)
        plt.rcParams.update({'text.color': (0.5,0.5,0.5),
                             'axes.labelcolor': (0.5,0.5,0.5)})
        self.fig.canvas.draw()
       

    def connect_signals_and_slots(self):

        self.map_settings_widget.borders_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_widget.states_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_widget.lakes_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_widget.rivers_checkbox.stateChanged.connect(self.update_feature_visibilities)
        self.map_settings_widget.coast_checkbox.stateChanged.connect(self.update_feature_visibilities)

        self.map_settings_widget.signal_offline_directory_changed.connect(self.draw_map)

        self.map_settings_widget.resolution_cb.currentTextChanged.connect(self.update_map)
        self.map_settings_widget.signal_colors_changed.connect(self.update_map)
        self.map_settings_widget.signal_background_changed.connect(self.update_map)
        self.map_settings_widget.signal_map_settings_changed.connect(self.update_map)

        self.extentWidget.sig_extent_changed.connect(self.set_map_extent)
        self.extentWidget.sig_set_to_global.connect(self.set_map_extent_to_global)
        self.extentWidget.sig_autoscale.connect(self.autoscale_plot)

        # these technically aren't qt signals and slots, these are matplotlib callback connections
        #self.fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
        #self.fig.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.fig.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        # self.fig.canvas.mpl_connect('scroll_event', self.scroll_ev_callback)

    def compute_figure(self):

        self.sph_proj = Geod(ellps='WGS84')
        self.projection = ccrs.PlateCarree()
        self.transform = ccrs.PlateCarree()

        self.axes = self.fig.add_subplot(1, 1, 1, projection=self.projection)

        c = '0.6'
        self.axes.tick_params(axis='both', labelsize=8, colors=c)
        self.axes.title.set_color(c)
        for spine in ['top', 'right', 'bottom', 'left']:
            self.axes.spines[spine].set_color(c)
        self.axes.set_xlabel('Detection Number')
        self.axes.set_ylabel('Distance')
        self.axes.xaxis.label.set_size(8)
        self.axes.xaxis.label.set_color(c)
        self.axes.xaxis.label.set_color(c)
        self.axes.yaxis.label.set_size(8)
        self.axes.yaxis.label.set_color(c)

    @pyqtSlot()
    def draw_map(self, preserve_extent=False):
        if preserve_extent:
            current_extent = self.axes.get_extent()

        self.axes.clear()

        if self.map_settings_widget.offline_checkbox.isChecked():
            # use the offline maps...
            cartopy.config['pre_existing_data_dir'] = self.map_settings_widget.offline_directory_label.text()
        else:
            cartopy.config['pre_existing_data_dir'] = ""

        resolution = self.map_settings_widget.resolution_cb.currentText()

        if self.map_settings_widget.backgroud_image_checkbox.isChecked():
            land_facecolor = 'none'
            ocean_facecolor = 'none'
            self.axes.stock_img()
        else:
            land_facecolor = (self.map_settings_widget.land_color_button.color().redF(), 
                            self.map_settings_widget.land_color_button.color().greenF(),
                            self.map_settings_widget.land_color_button.color().blueF())

            ocean_facecolor = (self.map_settings_widget.ocean_color_button.color().redF(), 
                           self.map_settings_widget.ocean_color_button.color().greenF(),
                           self.map_settings_widget.ocean_color_button.color().blueF())

        land = cfeature.NaturalEarthFeature('physical',
                                            'land',
                                            scale=self.map_settings_widget.resolution_cb.currentText(),
                                            edgecolor='face',
                                            facecolor=land_facecolor,
                                            linewidth=0.5)
        
        states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                        name='admin_1_states_provinces_lines',
                                                        scale=self.map_settings_widget.resolution_cb.currentText(),
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
        # toggle the visibility
        self.extentWidget.setVisible(self.extentWidget.isHidden())
        # set button color and hide other widgets
        if self.extentWidget.isVisible():
            style = "color: rgba(0,0,180,255);"
            self.hide_map_settings_widget()
        else:
            style = "color: rgba(20,20,20,255);"
            
    def showhide_map_settings_widget(self):
        # toggle the visibility 
        self.map_settings_widget.setVisible(self.map_settings_widget.isHidden())
        # set button color and hide other widgets
        if self.map_settings_widget.isVisible():
            style = "color: rgba(0,0,180,255);"
            self.hide_extent_widget()
        else:
            style = "color: rgba(20,20,20,255);"
        self.parent.toolButton_settings.setStyleSheet(style)

    def hide_map_settings_widget(self):
        # hide widget and reset the button color
        self.map_settings_widget.setVisible(False)

    def hide_extent_widget(self):
        # hide widget and reset the button color
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
            self.states.set_visible(self.map_settings_widget.states_checkbox.isChecked())
            self.lakes.set_visible(self.map_settings_widget.lakes_checkbox.isChecked())
            self.rivers.set_visible(self.map_settings_widget.rivers_checkbox.isChecked())
            self.borders.set_visible(self.map_settings_widget.borders_checkbox.isChecked())
            self.coast.set_visible(self.map_settings_widget.coast_checkbox.isChecked())
            self.fig.canvas.draw()  # update matlabplot

        except urllib.error.URLError:
            return
        
        except Exception as e:
            print(e)

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
        if self.map_settings_widget.show_grid_checkbox.isChecked():
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

        lat = self.parent.showgroundtruth.ev_widget.getLat()
        lon = self.parent.showgroundtruth.ev_widget.getLon()
        
        current_extent = self.axes.get_extent() # plotting the event should not change the extent
        self.gt_marker, = self.axes.plot(lon, lat, 'X', color='red', transform=self.transform, markersize=16, gid='ground_truth_marker')

        self.set_map_extent(current_extent)
        self.extentWidget.set_extent_spin_values(current_extent)     # update extentWidget

        self.show_hide_ground_truth(self.parent.showgroundtruth.ev_widget.showGT_cb.checkState())

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


    def update_range_max(self):
        self.update_detections(preserve_colors=True)

    def autoscale_plot(self):
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


        if self.parent.showgroundtruth.ev_widget.showGT_cb.isChecked():
            lons.append(self.parent.showgroundtruth.ev_widget.ev_lon_edit.value())
            lats.append(self.parent.showgroundtruth.ev_widget.ev_lat_edit.value())

        maxLat = max(lats + self.end_lats)
        minLat = min(lats + self.end_lats)
        center_lat = minLat + (maxLat - minLat)/2.

        maxLon = max(lons + self.end_lons)
        minLon = min(lons + self.end_lons)
        center_lon = minLon + (maxLon - minLon)/2.

        if maxLon != minLon:
            width = abs(maxLon - minLon)
        else:
            width = 50

        if maxLat != minLat:
            height = abs(maxLat - minLat)
        else:
            height = 50

        width, height = self.fix_aspect(width, height)

        width_adj = width * 0.10
        height_adj = height * 0.10

        minLat = center_lat - height/2. - height_adj
        maxLat = center_lat + height/2. + height_adj
        minLon = center_lon - width/2. - width_adj
        maxLon = center_lon + width/2. + width_adj

        # if there is only one detection, or they are all in line, the map could come up very thin, so 
        # we want to jump through some hoops to make sure the map is at least a little normal
        self.set_map_extent([minLon, maxLon, minLat, maxLat])

        self.extentWidget.set_extent_spin_values([minLon, maxLon, minLat, maxLat])     # update extentWidget


    def fix_aspect(self, w: float, h: float) -> tuple[float, float]:
        golden_ratio = 1.618
        #try to make the extent aspect ration to be something normal
        if abs(w/h) >= golden_ratio:
            new_h = w/golden_ratio
            new_w = w
        elif abs(w/h) <= 1./golden_ratio:
            new_w = h/golden_ratio
            new_h = h
        else:
            return w, h
        return new_w, new_h

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
            self.axes.set_title('Lon = {:+f}, Lat = {:+f}'.format(event.xdata, event.ydata), loc='center', pad=20, fontsize=10, color='0.6')
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

    def scroll_ev_callback(self, event):
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

