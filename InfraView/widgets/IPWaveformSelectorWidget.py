from PyQt5.QtWidgets import (QWidget, QAbstractButton, QButtonGroup, 
                             QCheckBox, QFormLayout, QLabel, QVBoxLayout)

from PyQt5.QtCore import pyqtSlot, pyqtSignal

from obspy.core.stream import Stream


class IPWaveformSelectorWidget(QWidget):

    sig_checkbox_clicked = pyqtSignal(list)

    checkbox_list = []
    name_list = []
    value_list = []

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent

        self.buildUI()
    
    def buildUI(self):

        self.button_group = QButtonGroup()
        self.button_group.setExclusive(False)
        self.button_group.buttonClicked.connect(self.update_what_is_displayed)
        
        self.title_label = QLabel("Show/Hide")

        self.form_layout = QFormLayout()

        self.tooltip_label = QLabel()
        tooltip_text = 'Note: This is for displaying/hiding plots only.  When you run the beamformer, all traces loaded will be used in the calculation. If you need to delete a trace, do it using the tabs in the trace viewer.'
        self.tooltip_label.setText(tooltip_text)
        self.tooltip_label.setWordWrap(True)

        self.main_layout = QVBoxLayout()
        self.main_layout.addWidget(self.title_label)
        self.main_layout.addLayout(self.form_layout)
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
        self.clear_form()

        for trace in new_stream:
            # All new traces are automatically set to display
            val = True # default value for a checkbox is True
            if trace.id in previous_name_list:
                # if entry was previously loaded, keep that setting
                idx = previous_name_list.index(trace.id)
                val = previous_value_list[idx]

            # append new info to the name and value lists
            self.name_list.append(trace.id)
            self.value_list.append(val)

            # make a new checkbox and add it to the list
            new_checkbox = QCheckBox()
            self.button_group.addButton(new_checkbox)
            new_checkbox.setChecked(val)
            self.checkbox_list.append(new_checkbox)
            self.form_layout.addRow(trace.id, new_checkbox)

    def clear_form(self):
        self.checkbox_list.clear()
        for row in reversed(range(self.form_layout.count())):
            item = self.form_layout.takeAt(row)
            item.widget().deleteLater()

    @pyqtSlot(QAbstractButton)
    def update_what_is_displayed(self, button):
        # the value list is used by pl_widget to keep track of what is clicked
        for idx, checkbox in enumerate(self.checkbox_list):
            if checkbox is button:
                self.value_list[idx] = checkbox.isChecked()
        
        # This tells the plotviewer to update which plots to show
        self.parent.pl_widget.draw_plots()

    def get_checkbox_list(self):
        return self.checkbox_list
    
    def get_value_list(self):
        return self.value_list
