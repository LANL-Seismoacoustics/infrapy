import numpy as np
import math

import pyqtgraph as pg

from PyQt5.QtWidgets import QComboBox, QLabel, QPushButton, QHBoxLayout, QFormLayout, QSpinBox
from PyQt5.QtCore import Qt, pyqtSignal

from InfraView.widgets import IPBaseWidgets

class IPSlownessSettingsWidget(IPBaseWidgets.IPSettingsWidget):

    def __init__(self, parent):
        super().__init__(parent)

        self.beamformingWidget = parent
        
        self.buildUI()

    def buildUI(self):
        self.update_button = QPushButton("Update")
        self.update_button.setMaximumWidth(100)
        self.update_button.setEnabled(False)
        self.update_button.clicked.connect(self.deactivate_update_button)
        self.update_button.clicked.connect(self.beamformingWidget.slownessPlot.update)

        colormap_label = QLabel("Color Map: ")
        self.colormap_cb = QComboBox()

        available_maps = pg.colormap.listMaps(source='matplotlib')
        self.colormap_cb.addItems(available_maps)
        self.colormap_cb.setCurrentText('jet')
        self.colormap_cb.currentTextChanged.connect(self.activate_update_button)

        resolution_label = QLabel("Resolution:")
        self.resolution_spin = QSpinBox()
        self.resolution_spin.setRange(10,1000)
        self.resolution_spin.setMaximumWidth(70)
        self.resolution_spin.setValue(300)
        self.resolution_spin.setToolTip("Number of points (horizontal and vertical) that make up the slowness image.\nIf you want to 'smooth' the plot, reduce the size of the trace velocity step \nsize and the azimuth step size in the beamformer settings.")
        self.resolution_spin.valueChanged.connect(self.activate_update_button)

        form1_layout = QFormLayout()
        form1_layout.addRow(colormap_label, self.colormap_cb)

        form2_layout = QFormLayout()
        form2_layout.addRow(resolution_label, self.resolution_spin)

        main_layout = QHBoxLayout()
        main_layout.addLayout(form1_layout)
        main_layout.addLayout(form2_layout)

        main_layout.addWidget(self.update_button)
        self.setLayout(main_layout)

    def settings(self):
        '''returns the current settings'''
        settings = {'cmap': self.colormap_cb.currentText()}

        return settings

    def activate_update_button(self):
        self.update_button.setEnabled(True)

    def deactivate_update_button(self):
        self.update_button.setEnabled(False)

class IPSlownessImageItem(pg.ImageItem):
    sig_info_changed = pyqtSignal(str)

    def __init__(self, parent):
        super().__init__()

        self.parent = parent
        self.resolution = 0
        self.hr = 1.
        self.pps = 1.
        self.tracev_range = ()

    def set_params(self, resolution, tracev_range, pps):
        self.resolution = resolution
        self.hr = resolution/2.
        self.pps = pps
        self.tracev_range = tracev_range

    def hoverEvent(self, event):
        if not event.isExit():
            pos = event.pos()

            x = (pos.x() - self.hr)/self.pps
            y = (pos.y() - self.hr)/self.pps

            slow = np.sqrt(x**2 + y**2)
            vel = 1./slow

            az = np.degrees(np.arctan2(x, y))

            if vel > self.tracev_range[0] and vel < self.tracev_range[1]:
                info_str = 'Velocity: {:3.2f} m/s \n Azimuth: {:3.2f} deg.'.format(vel, az)
                self.sig_info_changed.emit(info_str)

class IPSlownessPlot(pg.PlotItem):

    def __init__(self, parent):
        super().__init__()

        self.parent = parent

        self.resolution = 0
        self.hr = 1.
        self.image_item = IPSlownessImageItem(self)
        self.image_item.sig_info_changed.connect(self.update_info_label)
        self.tracev_range = ()

        self.addItem(self.image_item)

        self.showAxis('top')
        self.showAxis('right')

        self.vb.setAspectLocked(lock=True, ratio=1)

        # initialize circles
        self.i_circle = pg.QtWidgets.QGraphicsEllipseItem(0 , 0, 0, 0)
        self.i_circle.setPen(pg.mkPen(width=3, color='k'))
        self.addItem(self.i_circle)

        self.o_circle = pg.QtWidgets.QGraphicsEllipseItem(0 , 0, 0, 0)
        self.o_circle.setPen(pg.mkPen(width=3, color='k'))
        self.addItem(self.o_circle)

        self.info_label = pg.TextItem(color=(0, 0, 0), html=None, anchor=(1, 0))
        self.addItem(self.info_label, ignoreBounds=True)

        self.radial_list = []

        ax = self.getAxis('bottom')
        ax.setTicks([])
        ax = self.getAxis('top')
        ax.setTicks([])
        ax = self.getAxis('right')
        ax.setTicks([])
        ax = self.getAxis('left')
        ax.setTicks([])

        cmap = pg.colormap.get('jet', source='matplotlib')
        self.image_item.setColorMap(cmap)
    
    def set_image(self, image, resolution, tracev_range):

        self.resolution = resolution
        self.hr  = self.resolution/2.
        self.tracev_range = tracev_range

        self.pps = self.resolution / (2./self.tracev_range[0]) # points per 1/vel unit (points per slowness)

        self.setXRange(0, resolution)
        self.setYRange(0, resolution)

        self.image_item.setImage(image)
        self.image_item.set_params(self.resolution, self.tracev_range, self.pps)

        self.draw_radials()
        self.draw_circles()

        self.setAutoVisible(y=True, x=True)
        self.getViewBox().autoRange()

        myRange = self.viewRange()
        self.info_label.setPos(myRange[0][1], myRange[1][1])

    def update_info_label(self, info_str):
        self.info_label.setText(info_str)

    def draw_radials(self):
        count = 8

        angles = np.arange(0, 360, 360./count)

        # clear old radials
        for rline in self.radial_list:
            self.removeItem(rline)

        for angle in angles:
            r_line = pg.InfiniteLine(pos=(self.hr, self.hr), angle=angle, pen=pg.mkPen((100,100,100), width=1, style=Qt.DotLine))
            self.radial_list.append(r_line)
            self.addItem(r_line)

    def draw_circles(self):

        # outer circle
        self.o_circle.setRect(0, 0, self.resolution, self.resolution)

        # inner circle
        lx =  (1/self.tracev_range[0] - 1/self.tracev_range[1]) * self.pps   # lower x coord
        w = 2. * self.pps/self.tracev_range[1]                               # circle width
        self.i_circle.setRect(lx,lx,w,w)

        # circle_count = 5
        # for r in range(1, circle_count+1):
        #     circle = pg.QtWidgets.QGraphicsEllipseItem(-r / (min_trace_vel*circle_count), -r / (min_trace_vel*circle_count), r * 2 / (min_trace_vel*circle_count), r * 2 / (min_trace_vel*circle_count))
        #     circle.setPen(pg.mkPen(0.75, width=2))
        #     circle.setZValue(10)
        #     self.addItem(circle)

    def update(self):
        settings = self.parent.slownessSettings.settings()
        cmap = pg.colormap.get(settings['cmap'], source='matplotlib')
        self.image_item.setColorMap(cmap)

    def clear_slowness(self):
        self.image_item.clear()
        self.resolution = 0
        self.hr = 1.
        self.tracev_range = ()