from PyQt5.QtWidgets import QTableView
from PyQt5.QtCore import Qt, QAbstractTableModel


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


class IPDatabaseQueryResultsTable(QTableView):
    """
    Table to view (hopefully) generic results from a sql query
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent

    def setData(self, data):
        model = IPPandasModel(data)
        self.setModel(model)

    
