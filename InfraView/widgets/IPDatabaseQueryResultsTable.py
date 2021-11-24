from PyQt5.QtWidgets import QHBoxLayout, QPushButton, QTableView, QVBoxLayout, QWidget
from PyQt5.QtCore import Qt, QAbstractTableModel

import pandas as pd


class IPPandasModel(QAbstractTableModel):
    """
    class to populate a tableview with a pandas dataframe
    """
    def __init__(self, dataframe, parent=None):
        super().__init__()
        self.dataframe = dataframe

    def rowCount(self, parent=None):
        return self.dataframe.shape[0]

    def columnCount(self, parent=None):
        return self.dataframe.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self.dataframe.iloc[index.row()][index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.dataframe.columns[col]


class IPDatabaseQueryResultsTable(QWidget):
    """
    Table widget to view results from an sql query
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.buildUI()

    def buildUI(self):
        self.tableView = QTableView(self)

        self.clear_button = QPushButton("Clear Table")

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.clear_button)
        button_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addLayout(button_layout)
        main_layout.addWidget(self.tableView)
        self.setLayout(main_layout)

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.clear_button.clicked.connect(self.clearTable)

    def setData(self, data):
        '''This takes a pandas dataframe, and converts it for display in our tableView'''
        self.model = IPPandasModel(data)
        self.tableView.setModel(self.model)
        self.tableView.reset()

    def clearTable(self):
        # maybe do some additional clean-up here?
        self.model.deleteLater()


    
