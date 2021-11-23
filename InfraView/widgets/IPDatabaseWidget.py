from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout
from InfraView.widgets import IPDatabaseConnectWidget
from InfraView.widgets import IPDatabaseQueryWidget
from InfraView.widgets import IPDatabaseQueryResultsTable

class IPDatabaseWidget(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent  # reference to the IPBeamformingWidget to which this belongs
        self.ipdatabase_connect_widget = None
        self.ipdatabase_query_widget = None
        self.ipdatabase_query_results_table = None

        self.buildUI()

    def buildUI(self):
        self.ipdatabase_connect_widget = IPDatabaseConnectWidget.IPDatabaseConnectWidget(self)
        self.ipdatabase_query_widget = IPDatabaseQueryWidget.IPDatabaseQueryWidget(self)
        self.ipdatabase_query_results_table = IPDatabaseQueryResultsTable.IPDatabaseQueryResultsTable(self)
        
        hlayout = QHBoxLayout()
        hlayout.addWidget(self.ipdatabase_connect_widget)
        hlayout.addWidget(self.ipdatabase_query_widget)
        hlayout.addStretch()

        vlayout = QVBoxLayout()
        vlayout.addLayout(hlayout)
        vlayout.addWidget(self.ipdatabase_query_results_table)
        vlayout.addStretch()

        self.setLayout(vlayout)