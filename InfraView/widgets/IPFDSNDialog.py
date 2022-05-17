from PyQt5.QtWidgets import QDialog, QVBoxLayout, QDialogButtonBox
from PyQt5.QtCore import Qt

from InfraView.widgets import IPFDSNWidget


class IPFDSNDialog(QDialog):

    fdsnWidget = None

    def __init__(self, parent):
        super(IPFDSNDialog, self).__init__(parent)
        self.buildUI()

    def buildUI(self):
        self.setWindowTitle(self.tr('FDSN Import'))

        self.fdsnWidget = IPFDSNWidget.IPFDSNWidget()

        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(self.fdsnWidget)
        layout.addWidget(buttons)

        self.setLayout(layout)

        self.resize(400, self.height())

    def getStreams(self):
        # pass through to get the stream info
        return self.fdsnWidget.getStreams()

    def getInventory(self):
        # pass through to get the inventory info out
        return self.fdsnWidget.getInventory()
