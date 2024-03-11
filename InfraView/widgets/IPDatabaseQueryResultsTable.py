from PyQt5.QtWidgets import QHBoxLayout, QPushButton, QTableView, QVBoxLayout, QAbstractItemView, QFrame, QLabel, QSizePolicy
from PyQt5.QtCore import Qt, QAbstractTableModel, QModelIndex, QVariant, pyqtSignal, pyqtSlot

from obspy.core.stream import Stream
from obspy.core.event.origin import Origin

from obspy.core import UTCDateTime

from infrapy.utils import database
from InfraView.widgets import IPUtils


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

    def column_headers(self):
        return self.col_headers

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

    signal_new_stream_from_db = pyqtSignal(Stream, bool)      

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
        button_font = self.clear_button.font()
        button_font.setPointSize(10)
        self.clear_button.setFont(button_font)

        self.select_all_button = QPushButton("Select All")
        self.select_all_button.setFont(button_font)

        self.select_none_button = QPushButton("Select None")
        self.select_none_button.setFont(button_font)

        self.get_selected_rows_button = QPushButton("Get Selected")
        self.get_selected_rows_button.setFont(button_font)

        self.append_selected_button = QPushButton("Append Selected")
        self.append_selected_button.setFont(button_font)

        self.replace_with_selected_button = QPushButton("Replace with Selected")
        self.replace_with_selected_button.setFont(button_font)

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.clear_button)
        button_layout.addWidget(self.select_all_button)
        button_layout.addWidget(self.select_none_button)
        # button_layout.addWidget(self.get_selected_rows_button)
        button_layout.addWidget(self.append_selected_button)
        button_layout.addWidget(self.replace_with_selected_button)
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
        self.append_selected_button.clicked.connect(self.getSelected_append)
        self.replace_with_selected_button.clicked.connect(self.getSelected_replace)


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
    
    @pyqtSlot()
    def getSelected_append(self):
        st = self.getSelected()
        # this signal will connect to a slot in ApplicationWindow to assemble the streams and inventories and put them on the waveform widget.
        if st is not None:
            self.signal_new_stream_from_db.emit(st, True)

    def getSelected_replace(self):
        st = self.getSelected()
        # this signal will connect to a slot in ApplicationWindow to assemble the streams and inventories and put them on the waveform widget.
        if st is not None:
            self.signal_new_stream_from_db.emit(st, False)
    
    def getSelected(self):
        # if nothing is selected, then selectionModel() will return None
        # maybe this should be in the EventQueryWidget?

        if self.tableView.selectionModel():
            rows = self.tableView.selectionModel().selectedRows()
        else:
            return None

        if len(rows) == 0:
            IPUtils.errorPopup("No stations selected")
            return

        selected_rows = []
        for row in rows:
            selected_rows.append(row.row())

        selected_wds = []
        for idx, wd in enumerate(self.data):
            if idx in selected_rows:
                selected_wds.append(wd)

        # now assemble the output stream
        st = Stream(traces=None)
        starttime, stoptime = self.parent.ipdatabase_query_widget.get_startstop_times()
        tables= self.get_tables()

        db_tables = database.make_tables_from_dict(tables=tables, schema=self.get_schema())

        # we need the indexes of the stations and channel columns in the wfdisc table
        headers = self.model.column_headers()
        sta_idx = headers.index('sta')
        chn_idx = headers.index('chan')

        for wd in selected_wds:
            new_stream, _ = database.wvfrms_from_db(self.get_session(),
                                                    db_tables=db_tables,
                                                    stations=wd[sta_idx],
                                                    channel=wd[chn_idx],
                                                    starttime=starttime,
                                                    endtime=stoptime)

            st += new_stream

        return st
        

    def get_session(self):
        return self.parent.ipdatabase_connect_widget.session
    
    def get_schema(self):
        return self.parent.ipdatabase_connect_widget.schema_type_combo.currentText()

    def get_tables(self):
        table_dictionary = self.parent.ipdatabase_connect_widget.table_dialog.get_tables_from_text()
        return table_dictionary
        

class IPEventsModel(QAbstractTableModel):

    def __init__(self, origins, prefor, parent=None):
        super().__init__()

        self.origins = origins
        self.prefor = prefor
        self.col_headers = [c.name for c in self.origins[0].__table__.columns]

        print(self.prefor)

    def rowCount(self, parent=None):
        return len(self.origins)

    def columnCount(self, parent=None):
        return len(self.origins[0])

    def data(self, index, role):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self.origins[index.row()][index.column()])
            elif role == Qt.EditRole:
                return str(self.origins[index.row()][index.column()])

        return None

    def setData(self, index, value, role):
        if index.isValid():
            if role == Qt.EditRole:
                self.origins[index.row()][index.column()] = value
                self.editCompleted.emit(value)
                return True
            return False
        return False

    def get_data(self):
        return self.origins

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

    sig_origin_changed = pyqtSignal(dict)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        self.size_policy = self.sizePolicy()
        self.size_policy.setHorizontalPolicy(QSizePolicy.Expanding)
        self.setSizePolicy(self.size_policy)
        self.model = None
        self.parent = parent
        self.buildUI()

    def buildUI(self):
        
        title_label = QLabel("\tEvent Query Results")
        title_label.setStyleSheet("QLabel {font-weight:bold; color: white; background-color: black}")

        self.use_selected_button = QPushButton("Get Selected")
        self.clear_button = QPushButton("Clear Table")
        button_font = self.clear_button.font()
        button_font.setPointSize(10)
        self.use_selected_button.setFont(button_font)
        self.clear_button.setFont(button_font)

        self.tableView = QTableView(self)
        self.tableView.setSelectionMode(QAbstractItemView.SingleSelection)
        self.tableView.setSelectionBehavior(QAbstractItemView.SelectRows)

        horiz_layout_0 = QHBoxLayout()
        horiz_layout_0.addWidget(self.clear_button)
        horiz_layout_0.addWidget(self.use_selected_button)
        horiz_layout_0.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addWidget(title_label)
        main_layout.addLayout(horiz_layout_0)
        main_layout.addWidget(self.tableView)
        self.setLayout(main_layout)
        
        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.clear_button.clicked.connect(self.clearTable)
        self.use_selected_button.clicked.connect(self.useSelected)

    def setData(self, data, prefor):
        '''This takes Wfdisc rows, and converts it for display in our tableView'''
        self.model = IPEventsModel(data, prefor)
        self.tableView.setModel(self.model)
        self.tableView.reset()

        print(prefor[0].orid)
        start_idx = self.model.index(0,0)
        # print(self.model.data(start_idx, Qt.DisplayRole))
        #mdl_idxs = self.model.match(QModelIndex(), Qt.DisplayRole, prefor[0].orid)
        #print(mdl_idxs)
        for idx, origin in enumerate(data):
            if origin.orid == prefor[0].orid:
                self.tableView.selectRow(idx)
                break
        

    def clearTable(self):
        # maybe do some additional clean-up here?
        if self.model:
            self.model.deleteLater()

    def useSelected(self):
        if self.tableView.selectionModel():
            model_idx = self.tableView.selectionModel().selectedRows()
        else:
            return None
        
        for idx in model_idx:
            selected_origin = self.model.origins[idx.row()]
            print(type(selected_origin))
            new_origin ={}
            new_origin['Latitude'] = selected_origin.lat
            new_origin['Longitude'] = selected_origin.lon
            new_origin['UTC Date'] = UTCDateTime(selected_origin.time).date
            new_origin['UTC Time'] = UTCDateTime(selected_origin.time).time
            new_origin['Evid'] = selected_origin.evid
            new_origin['Orid'] = selected_origin.orid

            self.sig_origin_changed.emit(new_origin)



        # selected = self.model.get_data()[row.row()]
        # lat = selected[0]
        # lon = selected[1]
        # datetime = UTCDateTime(selected[3])
        # date = datetime.date
        # time = datetime.time
        # evid = selected[5]

        # event = {'Name': evid, 'UTC Date': date, 'UTC Time': time, 'Latitude': lat, 'Longitude':lon}
        # self.parent.parent.eventWidget.setEvent(event)

        # self.parent.parent.mainTabs.setCurrentIndex(4)




