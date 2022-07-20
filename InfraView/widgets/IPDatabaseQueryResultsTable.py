from PyQt5.QtWidgets import QHBoxLayout, QPushButton, QTableView, QVBoxLayout, QWidget, QAbstractItemView, QFrame, QLabel, QSizePolicy
from PyQt5.QtCore import Qt, QAbstractTableModel, QVariant

from obspy.core import read as obsRead


class IPWfdiscModel(QAbstractTableModel):
    """
    class to populate a tableview with rows of Wfdisc results
    """

    def __init__(self, wfs, parent=None):
        super().__init__()

        self.wfs = wfs
        self.col_headers = [c.name for c in self.wfs[0].__table__.columns]

    def rowCount(self, parent=None):
        return len(self.wfs)

    def columnCount(self, parent=None):
        return len(self.wfs[0])

    def data(self, index, role):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self.wfs[index.row()][index.column()])
            elif role == Qt.EditRole:
                return str(self.wfs[index.row()][index.column()])

        return None

    def setData(self, index, value, role):
        if index.isValid():
            if role == Qt.EditRole:
                self.wfs[index.row()][index.column()] = value
                self.editCompleted.emit(value)
                return True
            return False
        return False

    def flags(self, index):
        return Qt.ItemIsSelectable | Qt.ItemIsEnabled

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.col_headers[section]
            elif orientation == Qt.Vertical:
                return section
        return QVariant()


class IPPandasModel(QAbstractTableModel):
    """
    NOT CURRENTLY USED!
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
            # elif role == Qt.CheckStateRole:
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
            return False
        return False

    def flags(self, index):
        # return Qt.ItemIsEditable | QAbstractTableModel.flags(index)
        # return Qt.ItemIsEditable | Qt.ItemIsSelectable | Qt.ItemIsEnabled
        return Qt.ItemIsSelectable | Qt.ItemIsEnabled

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.dataframe.columns[section]
            elif orientation == Qt.Vertical:
                return self.dataframe.index[section]
        return QVariant()


class IPDatabaseQueryResultsTable(QFrame):
    """
    Table widget to view results from an sql query
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        self.data = None
        self.model = None
        self.parent = parent
        self.buildUI()

    def buildUI(self):
        
        title_label = QLabel("\tDatabase Query Results")
        title_label.setStyleSheet("QLabel {font-weight:bold; color: white; background-color: black}")

        self.tableView = QTableView(self)
        self.tableView.setSelectionBehavior(QAbstractItemView.SelectRows)

        self.clear_button = QPushButton("Clear Table")
        self.select_all_button = QPushButton("Select All")
        self.select_none_button = QPushButton("Select None")
        self.get_selected_rows_button = QPushButton("Get Selected")

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.clear_button)
        button_layout.addWidget(self.select_all_button)
        button_layout.addWidget(self.select_none_button)
        button_layout.addWidget(self.get_selected_rows_button)
        button_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addWidget(title_label)
        main_layout.addLayout(button_layout)
        main_layout.addWidget(self.tableView)
        self.setLayout(main_layout)

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        """I find it useful to group as many connections in one place as possible"""

        self.clear_button.clicked.connect(self.clearTable)
        self.select_none_button.clicked.connect(self.selectNone)
        self.select_all_button.clicked.connect(self.selectAll)
        self.get_selected_rows_button.clicked.connect(self.getSelected)

    def setData(self, data):
        '''This takes Wfdisc rows, and converts it for display in our tableView'''

        self.data = data
        self.model = IPWfdiscModel(data)
        self.tableView.setModel(self.model)
        self.tableView.reset()

    def clearTable(self):
        # maybe do some additional clean-up here?
        if self.model:
            self.model.deleteLater()

    def selectAll(self):
        self.tableView.selectAll()

    def selectNone(self):
        self.tableView.clearSelection()

    def getSelected(self):
        # if nothing is selected, then selectionModel() will return None
        if self.tableView.selectionModel():
            rows = self.tableView.selectionModel().selectedRows()
        else:
            return None

        selected_rows = []
        for row in rows:
            selected_rows.append(row.row())

        selected_wds = []
        for idx, wd in enumerate(self.data):
            if idx in selected_rows:
                selected_wds.append(wd)

        print("selected wfdisc length = {}".format(len(selected_wds)))



class IPEventsModel(QAbstractTableModel):
    """
    class to populate a tableview with rows of event results
    """

    def __init__(self, evs, parent=None):
        super().__init__()

        self.evs = evs
        self.col_headers = [c.name for c in self.evs[0].__table__.columns]

    def rowCount(self, parent=None):
        return len(self.evs)

    def columnCount(self, parent=None):
        return len(self.evs[0])

    def data(self, index, role):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self.evs[index.row()][index.column()])
            elif role == Qt.EditRole:
                return str(self.evs[index.row()][index.column()])
        return None

    def setData(self, index, value, role):
        if index.isValid():
            if role == Qt.EditRole:
                self.evs[index.row()][index.column()] = value
                self.editCompleted.emit(value)
                return True
            return False
        return False

    def flags(self, index):
        return Qt.ItemIsSelectable | Qt.ItemIsEnabled

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.col_headers[section]
            elif orientation == Qt.Vertical:
                return section
        return QVariant()


class IPEventQueryResultsTable(QFrame):
    """
    Table widget to view results from an sql event query (evid or name)
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        self.size_policy = self.sizePolicy()
        self.size_policy.setHorizontalPolicy(QSizePolicy.Expanding)
        self.setSizePolicy(self.size_policy)
        self.data = None
        self.model = None
        self.parent = parent
        self.buildUI()

    def buildUI(self):
        
        title_label = QLabel("\tEvent Query Results")
        title_label.setStyleSheet("QLabel {font-weight:bold; color: white; background-color: black}")

        self.clear_button = QPushButton("Clear Table")
        self.clear_button.clicked.connect(self.clearTable)

        self.tableView = QTableView(self)
        self.tableView.setSelectionBehavior(QAbstractItemView.SelectRows)

        horiz_layout_0 = QHBoxLayout()
        horiz_layout_0.addWidget(self.clear_button)
        horiz_layout_0.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addWidget(title_label)
        main_layout.addLayout(horiz_layout_0)
        main_layout.addWidget(self.tableView)
        self.setLayout(main_layout)

    def setData(self, data):
        '''This takes Wfdisc rows, and converts it for display in our tableView'''

        self.data = data
        self.model = IPEventsModel(data)
        self.tableView.setModel(self.model)
        self.tableView.reset()

    def clearTable(self):
        # maybe do some additional clean-up here?
        if self.model:
            self.model.deleteLater()