
import configparser, yaml, traceback

from PyQt5.QtWidgets import (QComboBox, QFileDialog, QFrame, QHBoxLayout, QTextEdit,
                             QLabel, QLineEdit, QPushButton, QVBoxLayout, QDialog, QDialogButtonBox,)

from PyQt5.QtCore import pyqtSlot, pyqtSignal, QTimer, Qt

from sqlalchemy.orm import Session
from pathlib import Path

from InfraView.widgets import IPUtils
from infrapy.utils import database

from InfraView.widgets import IPBaseWidgets

class IPDatabaseSettingsWidget(IPBaseWidgets.IPSettingsWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self.connect_widget = IPDatabaseConnectWidget(parent=self, title='Connection')

        layout = QHBoxLayout()
        layout.addWidget(self.connect_widget)
        layout.addStretch()

        self.setLayout(layout)

    def to_dict(self):
        s_dict = {}
        s_dict['connect'] = self.connect_widget.to_dict()
        return s_dict

    def from_dict(self, _):
        # TODO: add this to manually add settings
        pass

class IPDatabaseConnectWidget(IPBaseWidgets.IPSettingsGroupBox):

    sig_session_created = pyqtSignal(Session)
        
    def __init__(self, parent=None, title=""):
        super().__init__(parent=parent, title=title)

        self.session = None
        self.cnfgname = ""

        self.buildUI()
        
    def buildUI(self):

        self.load_config_button = QPushButton("Load Config...")

        self.save_current_button = QPushButton("Save Config...")
        self.save_current_button.setEnabled(False)

        self.save_default_button = QPushButton("Save Defaults...")
        self.save_default_button.setToolTip("Save current db settings to the $HOME/.lanl_network_config.yml file")

        self.table_dialog = IPTableDialog(self)
        self.show_tables_button = QPushButton("Tables...")

        self.env_vars_dialog = IPEnvVarDialog(self)
        self.show_env_vars_button = QPushButton("Env Vars...")

        self.schema_type_combo = QComboBox()
        self.schema_type_combo.addItem("KBCore")
        self.schema_type_combo.addItem("CSS3")

        self.url_edit = QLineEdit()
        self.url_edit.setPlaceholderText("URL")
        self.url_edit.setMinimumWidth(500)
        self.url_edit.setMaximumWidth(500)

        self.create_session_button = QPushButton("Create Session")

        self.close_session_button = QPushButton("Clear Session")

        self.test_connection_button = QPushButton("Test Connection")

        row1_layout = QHBoxLayout()
        row1_layout.addWidget(self.load_config_button)
        row1_layout.addWidget(self.save_current_button)
        row1_layout.addWidget(self.save_default_button)
        row1_layout.addStretch()

        row2_layout = QHBoxLayout()
        row2_layout.addWidget(self.show_tables_button)
        row2_layout.addWidget(self.show_env_vars_button)
        row2_layout.addWidget(self.schema_type_combo)
        row2_layout.addStretch()

        row3_layout = QHBoxLayout()
        row3_layout.addWidget(self.url_edit)
        row3_layout.addWidget(self.create_session_button)
        row3_layout.addWidget(self.test_connection_button)
        row3_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addLayout(row1_layout)
        main_layout.addLayout(row2_layout)
        main_layout.addLayout(row3_layout)
        main_layout.addStretch()

        self.setLayout(main_layout)

        # create dialogs here
        self.cnfg_dialog = QFileDialog()
        self.cnfg_dialog.setFileMode(QFileDialog.ExistingFile)
        self.cnfg_dialog.setNameFilter("(*.ini)")

        self.table_dialog = IPTableDialog(self)

        self.env_vars_dialog = IPEnvVarDialog(self)

        # connect signals and slots
        self.connect_signals_and_slots()

        # We provide the option of having a default net config file for database and fdsn settings in the
        # users home directory.  If it exists, load the information in to save the user from having to do it
        # each time.
        try:
            self.load_default_net_config()
        except FileNotFoundError:
            # dont load defaults, just go on with life i guess
            pass

    def to_dict(self):
        s_dict = {}
        s_dict['schema'] = self.schema_type_combo.currentText() 
        s_dict['url'] = self.url_edit.text()
        return s_dict

    def from_dict(self, s_dict):
        self.schema_type_combo.setCurrentText(s_dict['schema'])
        self.url_edit.setText(s_dict['url'])

    def connect_signals_and_slots(self):
        self.load_config_button.clicked.connect(self.load_cnfg)
        # self.save_current_button.clicked.connect(self.save_current_config)
        self.save_default_button.clicked.connect(self.save_default_net_config)
        self.show_tables_button.clicked.connect(self.show_tables_dialog)
        self.show_env_vars_button.clicked.connect(self.show_env_vars_dialog)
        self.create_session_button.clicked.connect(self.create_session)
        self.test_connection_button.pressed.connect(self.check_connection)

    def load_default_net_config(self):
        '''
        some applications can store database and fdsn information in a file $HOME/.lanl_network_config.yml.  This function will look to see if that file 
        exists, and if it does will open it and read in the information into a dictionary.
        '''
        file_path = Path.home() / ".lanl_network_config.yml"
        
        if file_path.is_file():

            with  open(str(file_path)) as ifile:
                try:
                    net_dict =  yaml.safe_load(ifile)
                except yaml.YAMLError as e:
                    traceback.print_exc()
                    IPUtils.errorPopup('Error reading {}.\nNetwork configuration skipped.'.format(str(file_path)))
                    return
        else:
            raise FileNotFoundError
        
        # the default network config file is a yml file and is different from the normal ini config files, so
        # we can just handle it here for now
        try:
            self.url_edit.setText(net_dict['database']['url'])
            self.schema_type_combo.setCurrentText(net_dict['database']['schema'])
            self.table_dialog.set_text_from_table_dict(net_dict['database']['tables'])
            self.env_vars_dialog.set_text_from_vars_dict(net_dict['database']['environment_vars'])
            
        except KeyError as e:
            print(f"Key error: {e}")
            # IPUtils.errorPopup("Poorly formed config file.\nMissing some information", title="Key Error")

    def save_default_net_config(self):
        # Method to save the current db settings to the default network configuration file $HOME/.lanl_network_config.yml
        file_path = Path.home() / ".lanl_network_config.yml"

        # new dict will contain the settings to be written
        new_dict = dict()
        new_dict['url'] = self.url_edit.text()
        new_dict['schema'] = self.schema_type_combo.currentText()
        new_dict['tables'] = self.table_dialog.get_tables_from_text()
        new_dict['environmentvars'] = self.env_vars_dialog.get_vars_from_text()

        if file_path.is_file():
            # There is an existing config file. First read it in, we will only overwrite the database portion here.

            # first verify the user wants to do this
            check_dialog = IPBaseWidgets.IPContinueDialog(self, "This will overwrite the database section in the default configuration file. Are you sure?", "Overwrite Default Config")

            if check_dialog.exec() == QDialog.Accepted:
                # read in existing settings, we will only overwrite the db section
                with open(str(file_path)) as ifile:
                    try:
                        net_dict = yaml.safe_load(ifile)
                    except yaml.YAMLError as e:
                        IPUtils.errorPopup('Error reading {}. Bailing out'.format(str(file_path)))
                        return
                net_dict['database'] = new_dict

        else:
            check_dialog = IPBaseWidgets.IPContinueDialog(self, "Default Config file ({}) does not currently exist.\n Would you like to create it?".format(str(file_path)))
            if check_dialog.exec() == QDialog.Accepted:
                try:
                    def_config_path = Path(__file__).resolve().parent.parent.parent.parent / 'infrapy' / 'resources' / 'default_lanl_network_config.yml'
                    with open(str(def_config_path)) as ifile:
                        net_dict = yaml.safe_load(ifile)
                        net_dict['database'] = new_dict
                except yaml.YAMLError as e:
                    IPUtils.errorPopup('Error reading {}.\nNo default configuration file created.'.format(str(def_config_path)))
                    return
            else:
                return

        with open(str(file_path), 'w') as ofile:
            yaml.dump(net_dict, ofile, default_flow_style=False)


    @pyqtSlot()
    def load_cnfg(self):
        # This loads configuration settings from the standard infrapy ini files.
        if self.cnfg_dialog.exec_():
            self.cnfgname = self.cnfg_dialog.selectedFiles()[0]
            try:
                config = configparser.ConfigParser()
                config.read(self.cnfgname)
                self.schema_type_combo.setCurrentText(config['DATABASE']['schema'])
                self.url_edit.setText(config['DATABASE']['url'])

                self.table_dialog.set_text_from_table_dict(config['DBTABLES'])

                self.env_vars_dialog.set_text_from_vars_dict(config['DBENVIRONMENTVARS'])

                self.save_current_button.setEnabled(False)

            except Exception as e:
                IPUtils.errorPopup("Error reading config file \n{}".format(str(e)))

    @pyqtSlot()
    def show_tables_dialog(self):
        if self.table_dialog.exec_():
            self.save_current_button.setEnabled(True)

    @pyqtSlot()
    def show_env_vars_dialog(self):
        if self.env_vars_dialog.exec_():
            self.save_current_button.setEnabled(True)

    @pyqtSlot()
    def create_session(self):
        # first, if there is already an active session, close it...
        self.close_session()

        # make sure any environment variables are loaded.
        env_vars = self.env_vars_dialog.get_vars_from_text()
        if env_vars:
            database.set_db_env_variables(env_vars)

        url = self.url_edit.text()

        try:
            self.session = database.db_connect_url(url)
            self.sig_session_created.emit(self.session)
            self.url_edit.setStyleSheet("color: green")
        except ValueError as e:
            self.session = None
            self.url_edit.setStyleSheet("color: red")
            IPUtils.errorPopup("Error creating session.  \nMake sure you can reach the database and that the displayed url is correct.")
    
    def close_session(self):
        if self.session is not None:
            self.session.close()
            self.url_edit.setStyleSheet("color: black")
            self.session = None

    @pyqtSlot()
    def clear_form(self):
        if self.session is not None:
            self.close_session()

        self.schema_type_combo.setCurrentIndex(0)
        self.url_edit.setText("")

    @pyqtSlot()
    def reset_connection_colors(self):
        self.test_connection_button.setText('Test Connection')
        self.test_connection_button.setStyleSheet('QPushButton {color: black}')

    @pyqtSlot()
    def check_connection(self):
        if self.session is not None:
            if database.check_connection(self.session):
                self.test_connection_button.setText("Good Connection")
                self.test_connection_button.setStyleSheet('QPushButton {color: green}')
                QTimer.singleShot(3000, self.reset_connection_colors)
            else:
                self.test_connection_button.setText("Bad Connection")
                self.test_connection_button.setStyleSheet('QPushButton {color: red}')

                QTimer.singleShot(3000, self.reset_connection_colors)
        else:
            self.test_connection_button.setText("No active session")
            QTimer.singleShot(3000, self.reset_connection_colors)




class IPTableDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.buildUI()
    
    def buildUI(self):
        self.setWindowTitle("InfraView: Table Editor")
        self.tables_textEdit = QTextEdit()

        descriptor_label = QLabel('''The format for the tables should take the form of\nsite_descriptor: owner.tablename\n\n\tsite: global.site\n\twfdisc: global.wfdisc_raw\n\torigin: myorigin.origin\n\tevent: myevent.event\n''')
        
        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        vlayout = QVBoxLayout()
        vlayout.addWidget(self.tables_textEdit)
        vlayout.addWidget(descriptor_label)
        vlayout.addWidget(buttons)

        self.setLayout(vlayout)

    def exec_(self):
        self.initial_text = self.tables_textEdit.toPlainText().rstrip()       
        return super().exec_()

    def reset(self):
        self.tables_textEdit.setText(self.initial_text)
    
    def set_text_from_table_dict(self, tables):
        # set the text of the table editor from a dictionary of tables
        text = ""
        for key, value in tables.items():
            text += key + ':' + str(value) + '\n'
            
        self.tables_textEdit.setText(text)

    def get_tables_from_text(self):
        # Convert the QTextEdit contents into individual lines, then seperate into a dictionary
        lines = self.tables_textEdit.toPlainText().split("\n")

        table_dict = dict()

        for line in lines:
            if line.strip():    # if line is blank, this will skip it
                key_val = line.split(':')
                try:
                    table_dict[key_val[0]] = key_val[1]
                except IndexError:
                    table_dict[key_val[0]] = ""

        return table_dict 

    def reject(self):
        self.reset()
        super().reject()

class IPEnvVarDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.buildUI()
    
    def buildUI(self):
        self.setWindowTitle("InfraView: Database Env Variables")
        self.vars_textEdit = QTextEdit()

        descriptor_label = QLabel("Enter database specific environment variables here.\nThese can be entered here or in the appropriate shell rc file.\nIf entered here, they will last only for the duration of the session.")
        entry_label = QLabel("Entries should have the variable and value seperated by a colon:\n\n\tTNS_ADMIN: $ORACLE_HOME/network/admin\n\tORACLE_PORT: 1523\n")
        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        vlayout = QVBoxLayout()
        vlayout.addWidget(self.vars_textEdit)
        vlayout.addWidget(descriptor_label)
        vlayout.addWidget(entry_label)
        vlayout.addWidget(buttons)

        self.setLayout(vlayout)

    def exec_(self):
        self.initial_text = self.vars_textEdit.toPlainText().rstrip()       
        return super().exec_()

    def reset(self):
        self.vars_textEdit.setText(self.initial_text)
    
    def set_text_from_vars_dict(self, vars):
        if vars is None:
            return  # Nothing to do
        
        # set the text of the table editor from a dictionary of tables
        text = ""
        for key, value in vars.items():
            text += key + ':' + str(value) + '\n'

        self.vars_textEdit.setText(text)

    def get_vars_from_text(self):
        lines = self.vars_textEdit.toPlainText().split("\n")      # the rstrip removes trailing newlines etc
        
        vars_dict = dict()

        for line in lines:
            if line.strip():
                key_val = line.split(':')
                try:
                    vars_dict[key_val[0]] = key_val[1].strip()
                except IndexError as e:
                    vars_dict[key_val[0]] = ""
        return vars_dict

    def reject(self):
        self.reset()
        super().reject()