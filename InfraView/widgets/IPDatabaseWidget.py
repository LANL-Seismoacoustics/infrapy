from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout
from PyQt5.QtCore import Qt

from InfraView.widgets import IPBaseWidgets
from InfraView.widgets import IPDatabaseQueryWidget
from InfraView.widgets import IPDatabaseQueryResultsTable
from InfraView.widgets import IPUtils

class IPDatabaseWidget(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent  
        self.ipdatabase_settings_widget = None
        self.ipdatabase_query_widget = None
        self.ipdatabase_query_results_table = None
        self.ipev_query_widget = None
        self.ipev_query_results_table = None

        self.buildUI()

    def buildUI(self):
        self.ipdatabase_query_widget = IPDatabaseQueryWidget.IPDatabaseQueryWidget(self)
        self.ipev_query_widget = IPDatabaseQueryWidget.IPEventQueryWidget(self)

        self.ipdatabase_query_results_table = IPDatabaseQueryResultsTable.IPDatabaseQueryResultsTable(self)
        self.ipev_query_results_table = IPDatabaseQueryResultsTable.IPEventQueryResultsTable(self)
        
        hlayout = QHBoxLayout()
        hlayout.addWidget(self.ipev_query_widget)
        hlayout.addWidget(self.ipev_query_results_table)
        
        #IPSplitter only accepts widgets, so we need to put the hlayout into one
        top_widget = QWidget()
        top_widget.setLayout(hlayout)

        wave_widget = QWidget()
        wave_layout = QHBoxLayout()
        wave_layout.addWidget(self.ipdatabase_query_widget)
        wave_layout.addWidget(self.ipdatabase_query_results_table)
        wave_widget.setLayout(wave_layout)

        vertical_splitter = IPBaseWidgets.IPSplitter(Qt.Vertical)
        vertical_splitter.addWidget(top_widget)
        vertical_splitter.addWidget(wave_widget)   
        
        main_layout = QVBoxLayout()
        main_layout.addWidget(vertical_splitter)
        self.setLayout(main_layout)

    def set_controlling_widget(self, widget):
        self.ipdatabase_settings_widget = widget
        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.ipdatabase_settings_widget.connect_widget.sig_session_created.connect(self.ipdatabase_query_widget.set_session)
        self.ipdatabase_settings_widget.connect_widget.sig_session_created.connect(self.ipev_query_widget.set_session)

        self.ipev_query_results_table.sig_origin_changed.connect(self.ipdatabase_query_widget.update_time)                                                       