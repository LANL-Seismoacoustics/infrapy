from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QSplitter
from PyQt5.QtCore import Qt
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
        self.ipevent_query_widget = None
        self.ipevent_query_results_table = None

        self.buildUI()

    def buildUI(self):

        self.ipdatabase_connect_widget = IPDatabaseConnectWidget.IPDatabaseConnectWidget(self)
        self.ipdatabase_connect_widget.setMinimumWidth(400)
        self.ipdatabase_query_widget = IPDatabaseQueryWidget.IPDatabaseQueryWidget(self)
        self.ipevent_query_widget = IPDatabaseQueryWidget.IPEventQueryWidget(self)
        self.ipdatabase_query_results_table = IPDatabaseQueryResultsTable.IPDatabaseQueryResultsTable(self)
        
        self.ipevent_query_results_table = IPDatabaseQueryResultsTable.IPEventQueryResultsTable(self)
        hlayout = QHBoxLayout()
        hlayout.addWidget(self.ipdatabase_connect_widget)
        hlayout.addWidget(self.ipdatabase_query_widget)
        hlayout.addWidget(self.ipevent_query_widget)
        hlayout.addWidget(self.ipevent_query_results_table)
        
        #QSplitter only accepts widgets, so we need to put the hlayout into one
        top_widget = QWidget()
        top_widget.setLayout(hlayout)

        vertical_splitter = QSplitter(Qt.Vertical)
        vertical_splitter.addWidget(top_widget)
        vertical_splitter.addWidget(self.ipdatabase_query_results_table)

        #vlayout = QVBoxLayout()
        #vlayout.addLayout(hlayout)
        #vlayout.addWidget(self.ipdatabase_query_results_table)
        
        main_layout = QVBoxLayout()
        main_layout.addWidget(vertical_splitter)
        self.setLayout(main_layout)