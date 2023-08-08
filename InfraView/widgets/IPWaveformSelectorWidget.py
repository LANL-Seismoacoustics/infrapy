from PyQt5.QtWidgets import (QWidget, QAbstractButton, QButtonGroup, QPushButton, 
                             QCheckBox, QGridLayout, QLabel, QVBoxLayout)

from PyQt5.QtCore import pyqtSlot, pyqtSignal, Qt

from obspy.core.stream import Stream


class IPWaveformSelectorWidget(QWidget):

    sig_checkbox_clicked = pyqtSignal(list)
    sig_remove_trace_by_id = pyqtSignal(str)
    sig_remove_station_by_name = pyqtSignal(str)

    checkbox_list = []
    del_button_list = []
    name_list = []
    value_list = []

    # button ids are used to keep track of buttons and corresponding checkboxes
    # the pushbutton will have an id that is the negative of the checkbox.  Since
    # this id can't be -1, we start counting at 2.  Maybe a better way to do this?
    button_id = 2 

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent

        self.buildUI()
    
    def buildUI(self):

        self.showhide_button_group = QButtonGroup()
        self.showhide_button_group.setExclusive(False)
        self.showhide_button_group.idClicked.connect(self.update_what_is_displayed)

        self.del_button_group = QButtonGroup()
        self.del_button_group.idClicked.connect(self.del_button_clicked)

        self.grid_layout = QGridLayout()

        self.tooltip_label = QLabel()
        tooltip_text = 'Select which waveforms to show/hide with the checkbox on the left. Remove the waveforms from the plot with the buttons on the right.'
        self.tooltip_label.setText(tooltip_text)
        self.tooltip_label.setWordWrap(True)
        self.tooltip_label.setVisible(False)

        self.main_layout = QVBoxLayout()
        self.main_layout.addLayout(self.grid_layout)
        self.main_layout.addStretch()
        self.main_layout.addWidget(self.tooltip_label)
        
        self.setLayout(self.main_layout)

    def update_selections(self, new_stream):
        # This is called when waveforms are added or removed
        # It updates the items in the checkbox list
        
        # before we clear everything, lets first record the current elements in the name dictionary
        # so we can keep track of what is currently visible/not visible.  That way we can preserve 
        # their settings

        if new_stream is None:
            return  # nothing to do

        previous_name_list = self.name_list.copy()
        previous_value_list = self.value_list.copy()
 
        # now clear everything
        self.name_list.clear()
        self.value_list.clear()
        self.clear_layout()

        self.grid_layout.addWidget(QLabel("Show/Hide"), 0, 0)
        self.grid_layout.addWidget(QLabel("Remove"), 0, 1)

        for i, trace in enumerate(new_stream):
            # All new traces are automatically set to display
            val = True # default value for a checkbox is True
            if trace.id in previous_name_list:
                # if entry was previously loaded, keep that setting
                idx = previous_name_list.index(trace.id)
                val = previous_value_list[idx]

            self.value_list.append(val)
            # make a new checkbox and add it to the list
            showhide_checkbox = QCheckBox(trace.id)
            showhide_checkbox.setChecked(val)

            self.checkbox_list.append(showhide_checkbox)
            self.showhide_button_group.addButton(showhide_checkbox)
            self.showhide_button_group.setId(showhide_checkbox, self.button_id)

            del_button = QPushButton('x')
            del_button.setToolTip('Remove {} from plot'.format(trace.id))
            del_button.setMaximumWidth(50) # makes it square?
            self.del_button_group.addButton(del_button)
            self.del_button_group.setId(del_button, -self.button_id)
            self.del_button_list.append(del_button)

            # increment the button id for the next one.
            self.button_id += 1

            self.grid_layout.addWidget(showhide_checkbox, i+1 , 0)   # the plus one is because the title is in row 0
            self.grid_layout.addWidget(del_button, i+1, 1)

            self.tooltip_label.setVisible(True)

    def clear_layout(self):
        self.checkbox_list.clear()
        self.tooltip_label.setVisible(False)

        for idx in reversed(range(self.grid_layout.count())):
            item = self.grid_layout.takeAt(idx)
            item.widget().deleteLater()

    @pyqtSlot(int)
    def del_button_clicked(self, id):
        clicked_checkbox = self.showhide_button_group.button(-id)
        self.sig_remove_trace_by_id.emit(clicked_checkbox.text())
        print("removing station {}".format(clicked_checkbox.text()))
        self.sig_remove_station_by_name.emit(clicked_checkbox.text())

    @pyqtSlot(int)
    def update_what_is_displayed(self, id):
        # the value list is used by pl_widget to keep track of what is clicked
        checked_button = self.showhide_button_group.button(id)
        idx = self.checkbox_list.index(checked_button)
        self.value_list[idx] = checked_button.isChecked()
        
        # This tells the plotviewer to update which plots to show
        self.parent.pl_widget.draw_plots()

    def get_checkbox_list(self):
        return self.checkbox_list
    
    def get_value_list(self):
        return self.value_list
