import pandas as pd

from PyQt5.QtWidgets import QMenu, QTableView
from PyQt5 import QtCore
from PyQt5.QtCore import QAbstractTableModel, Qt, pyqtSignal, pyqtSlot


class IPDetectionTableView(QTableView):

    signal_delete_detections = pyqtSignal(list)

    columns = ['Name',
               'Time (UTC)',
               'F Stat.',
               'Trace Vel. (m/s)',
               'Back Azimuth',
               'Latitude',
               'Longitude',
               'Elevation (m)',
               'Start',
               'End',
               'Freq Range',
               'Array Dim.',
               'Method',
               'Event',
               'Note']

    def __init__(self, parent=None):
        super().__init__(parent)

        self.parent = parent

        # initialize with empty dataframe
        self.pandas_table_model = IPPandasTableModel(pd.DataFrame(columns=self.columns))
        self.setModel(self.pandas_table_model)

        self.setContextMenuPolicy(Qt.CustomContextMenu)

        # hide the last column since it has reference to the line object
        # self.hideColumn(self.pandas_table_model.columnCount() - 1)

        self.horizontalHeader().setStretchLastSection(True)

        # we want to have a custom menu pop up when someone clicks on the row headers
        self._vertHeader = self.verticalHeader()
        self._vertHeader.setContextMenuPolicy(Qt.CustomContextMenu)
        self._vertHeader.setToolTip('You can click, Ctrl+click, or Shift+click the row\n'
                                    'number to select row(s). Then right click the row \n'
                                    'header to delete the detections')

        # set up a selection model for the table
        self.selectionModel = self.selectionModel()

        # connect signals and slots
        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):

        # decide if you want to be able to sort columns...
        # self.horizontalHeader().sectionClicked.connect(self.sort)

        self.signal_delete_detections.connect(self.pandas_table_model.delete_rows)
        self.customContextMenuRequested.connect(self.showContextMenu)
        self._vertHeader.customContextMenuRequested.connect(self.showHeaderContextMenu)

    def showContextMenu(self, position):
        return

    def showHeaderContextMenu(self, position):
        row = self._vertHeader.logicalIndexAt(position)

        menu = QMenu()

        deleteRows = menu.addAction("Delete Selected")
        ret = menu.exec_(self.mapToGlobal(position))
        if ret == deleteRows:
            # get the indexes of the selected rows
            rows_to_delete = [row]
            if self.selectionModel.hasSelection():
                selection = self.selectionModel.selectedRows()
                for item in selection:
                    if item.row() not in rows_to_delete:
                        rows_to_delete.append(item.row())
            # tell the pandas model to drop those rows, and update
            self.signal_delete_detections.emit(rows_to_delete)
            self.clearSelection()

    @pyqtSlot(int)
    def sort(self, idx):
        self.pandas_table_model.sort(idx, Qt.DescendingOrder)

    def get_model(self):
        return self.pandas_table_model

    def get_dataframe(self):
        return self.pandas_table_model._df

    def set_data(self, new_data):
        self.pandas_table_model.set_data(new_data)


class IPPandasTableModel(QAbstractTableModel):

    _current_sort_column = 1
    _current_sort_order = Qt.AscendingOrder
    _column_order = None

    _df = None

    def __init__(self, df=pd.DataFrame, parent=None):
        QAbstractTableModel.__init__(self, parent=parent)
        self._df = df
        self._column_order = list(self._df.columns)     # this sets the column order to the default order defined in the detectionview

    def headerData(self, section, orientation=QtCore.Qt.Horizontal, role=QtCore.Qt.DisplayRole):
        if role != QtCore.Qt.DisplayRole:
            return None

        if orientation == QtCore.Qt.Horizontal:
            try:
                return self._df.columns.tolist()[section]
            except(IndexError, ):
                return None

        elif orientation == QtCore.Qt.Vertical:
            try:
                return self._df.index.tolist()[section]
            except(IndexError):
                return None

        if not index.isValid():
            return None

        if self._df is None:
            return None

        return str(self._df.iloc[index.row(), index.column()])

    def data(self, index, role=QtCore.Qt.DisplayRole):

        if role != QtCore.Qt.DisplayRole:
            return None

        if not index.isValid():
            return None

        if self._df is None:
            return None

        return str(self._df.iloc[index.row(), index.column()])

    @pyqtSlot(pd.DataFrame)
    def set_data(self, new_df):
        self.layoutAboutToBeChanged.emit()

        self._df = new_df
        # self.sort(self._current_sort_column, self._current_sort_order)
        # self.setColumnOrder()

        self.layoutChanged.emit()

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self._df.index)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self._df.columns)

    def sort(self, column, order):

        colnames = self._df.columns.tolist()[int(column)]

        # flip the sort order if a column is reclicked
        if self._current_sort_column == column:
            if self._current_sort_order == Qt.AscendingOrder:
                self._current_sort_order = Qt.DescendingOrder
            elif self._current_sort_order == Qt.DescendingOrder:
                self._current_sort_order = Qt.AscendingOrder
        else:
            self._current_sort_column = column
            self._current_sort_order = order

        self.layoutAboutToBeChanged.emit()

        self._df.sort_values(colnames, ascending=self._current_sort_order, inplace=True)
        self._df.reset_index(inplace=True, drop=True)

        self.layoutChanged.emit()

    def append(self, arrival):
        # build a temp dataframe from the arrival
        df_new = pd.DataFrame([arrival], columns=arrival.keys())

        self.layoutAboutToBeChanged.emit()

        self._df = pd.concat([self._df, df_new], axis=0).reset_index()

        # re-sort so new value will be in correct place
        # self.sort(self._current_sort_column, self._current_sort_order)
        # self.setColumnOrder()

        self.layoutChanged.emit()

    @pyqtSlot(list)
    def delete_rows(self, rows):

        self.layoutAboutToBeChanged.emit()

        self._df = self._df.drop(rows)
        self._df.reset_index(inplace=True, drop=True)

        self.layoutChanged.emit()

    def setColumnOrder(self):
        self._df = self._df[self._column_order]
