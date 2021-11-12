from PyQt5.QtWidgets import QWidget, QVBoxLayout
from InfraView.widgets import IPDatabaseConnectWidget
from InfraView.widgets import IPDatabaseQueryWidget

class IPDatabaseWidget(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent  # reference to the IPBeamformingWidget to which this belongs
        self.ipdatabase_connect_widget = None
        self.ipdatabase_query_widget = None
        self.buildUI()

    def buildUI(self):
        self.ipdatabase_connect_widget = IPDatabaseConnectWidget.IPDatabaseConnectWidget(self)
        self.ipdatabase_query_widget = IPDatabaseQueryWidget.IPDatabaseQueryWidget(self)
        layout = QVBoxLayout()
        layout.addWidget(self.ipdatabase_connect_widget)
        layout.addWidget(self.ipdatabase_query_widget)
        layout.addStretch()

        self.setLayout(layout)