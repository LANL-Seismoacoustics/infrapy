from PyQt5.QtWidgets import QHBoxLayout, QPushButton, QTableView, QVBoxLayout, QWidget, QAbstractItemView
from PyQt5.QtCore import Qt, QAbstractTableModel, QVariant, QItemSelectionModel

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

    def data(self, index, role):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self.dataframe.iloc[index.row()][index.column()])
            #elif role == Qt.CheckStateRole:
            #    if (index.row() == 1 and index.column() == 2):
            #        return Qt.Checked
            elif role == Qt.EditRole:
                return str(self.dataframe.iloc[index.row()][index.column()])
        return None

    def setData(self, index, value, role):
        if index.isValid():
            if role == Qt.EditRole:
                self.dataframe.iat[index.row(), index.column()] = value
                self.editCompleted.emit(value)
                return True
            else:
                return False
        else:
            return False

    def flags(self, index):
        # return Qt.ItemIsEditable | QAbstractTableModel.flags(index)
        #return Qt.ItemIsEditable | Qt.ItemIsSelectable | Qt.ItemIsEnabled
        return Qt.ItemIsSelectable | Qt.ItemIsEnabled
        
    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.dataframe.columns[section]
            elif orientation == Qt.Vertical:
                return self.dataframe.index[section]
        return QVariant()


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
        self.tableView.setSelectionBehavior(QAbstractItemView.SelectRows)

        self.clear_button = QPushButton("Clear Table")
        self.select_all_button = QPushButton("Select All")
        self.select_none_button = QPushButton("Select None")

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.clear_button)
        button_layout.addWidget(self.select_all_button)
        button_layout.addWidget(self.select_none_button)
        button_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addLayout(button_layout)
        main_layout.addWidget(self.tableView)
        self.setLayout(main_layout)

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.clear_button.clicked.connect(self.clearTable)
        self.select_none_button.clicked.connect(self.selectNone)
        self.select_all_button.clicked.connect(self.selectNone)

    def setData(self, data):
        '''This takes a pandas dataframe, and converts it for display in our tableView'''
        self.model = IPPandasModel(data)
        self.tableView.setModel(self.model)
        self.tableView.reset()

    def clearTable(self):
        # maybe do some additional clean-up here?
        self.model.deleteLater()

    def selectAll(self):
        self.tableView.selectAll()

    def selectNone(self):
        self.tableView.clearSelection()


    
