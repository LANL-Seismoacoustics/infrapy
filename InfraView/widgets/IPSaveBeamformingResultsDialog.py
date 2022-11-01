
from PyQt5.QtWidgets import (QCheckBox, QDialog, QFileDialog, QLabel, QFormLayout, QHBoxLayout, 
                             QVBoxLayout, QPushButton, QDialogButtonBox, QLineEdit, QMessageBox,
                             QWidget, QSizePolicy)
from PyQt5.QtCore import Qt

from pathlib import Path

class IPSaveBeamformingResultsDialog(QDialog):

    path = None

    def __init__(self, parent, directory=None):
        super(IPSaveBeamformingResultsDialog, self).__init__(parent)
        
        self.buildUI()

    def buildUI(self):
        self.setWindowTitle(self.tr('InfraView: Save Results'))

        main_label = QLabel("Export beamforming results to a csv file")

        self.filepath_line = IPFileBrowseLine(self, self.path)
        form_layout = QFormLayout()
        form_layout.setLabelAlignment(Qt.AlignCenter)
        
        file_label = QLabel("Filename: ")
        file_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        form_layout.addRow(file_label, self.filepath_line)

        self.export_waveform_checkbox = QCheckBox()
        form_layout.addRow(QLabel('Export waveform csv: '), self.export_waveform_checkbox )
        self.waveform_line = IPFileBrowseLine(self, self.path)
        self.waveform_line.setEnabled(False)
        wave_label = QLabel("Waveform filename:")
        wave_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        form_layout.addRow(wave_label, self.waveform_line)
        

        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal, self)
        buttons.button(QDialogButtonBox.Ok).setText("Export");
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(main_label)
        layout.addLayout(form_layout)
        layout.addWidget(buttons)

        self.setLayout(layout)

        self.export_waveform_checkbox.stateChanged.connect(self.waveform_line.setEnabled)

        self.resize(600, self.height())

    def exec_(self, directory=None):
        if directory is None:
            self.path = Path.home()     
        else:
            self.path = directory

        if not self.path.is_dir():
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Question)
            msg.setText("The current save directory for beamforming results doesn't seem to exist. \nThis could happen if you are in an older Project.\n  This folder will be created now.\n (" + str(self.path) + ")")
            msg.exec_()
            self.path.mkdir(parents=True, exist_ok=True)

        path_text = str(self.path) + "/results.csv"
        wave_path_text = str(self.path) + "/waveform.csv"

        self.filepath_line.setText(path_text)
        self.filepath_line.setDirectory(self.path)

        self.waveform_line.setText(wave_path_text)
        self.waveform_line.setDirectory(self.path)

        return super().exec_()

    def wavefileIsChecked(self):
        return self.export_waveform_checkbox.isChecked()

    def getWaveFilename(self):
        return self.waveform_line.getText()

    def getFilename(self):
        return self.filepath_line.getText()

class IPFileBrowseLine(QWidget):
    # This is a combination of a line edit and a button to launch a file browser
    def __init__(self, parent, directory):
        super().__init__(parent)

        self.directory = directory

        self.text_edit = QLineEdit()
        self.text_edit.setMaximumWidth(600)
        self.text_edit.setFixedWidth(300)
        browse_button = QPushButton("Browse...")
        browse_button.clicked.connect(self.run_filedialog)

        layout = QHBoxLayout()
        layout.addWidget(self.text_edit)
        layout.addWidget(browse_button)
        self.setLayout(layout)

        self.my_file_dialog = QFileDialog()

    def setText(self, text):
        self.text_edit.setText(text)

    def getText(self):
        return self.text_edit.text()

    def setDirectory(self, directory):
        self.directory = directory

    def run_filedialog(self):
        fname = self.my_file_dialog.getSaveFileName(self, str(self.directory), filter="CSV (*.csv)")
        self.text_edit.setText(fname[0])

