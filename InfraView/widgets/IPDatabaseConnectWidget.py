from ssl import OP_NO_RENEGOTIATION
from tkinter import N
import urllib.parse
import configparser

from PyQt5.QtWidgets import (QComboBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QFrame, QHBoxLayout, 
                             QLabel, QLineEdit, QMessageBox, QPushButton, QTextEdit, QVBoxLayout, QSizePolicy)
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtCore import QRegExp, pyqtSlot, QTimer, Qt
from infrapy.utils import database

class IPTableDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.buildUI()
    
    def buildUI(self):
        self.setWindowTitle("InfraView: Table Editor")
        self.tables_textEdit = QTextEdit()

        descriptor_label = QLabel('''The format for the tables should look something like...\nsite: global.site\nwfdisc: global.wfdisc_raw\norigin: myorigin.origin\nevent: myevent.event''')
        
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
            text += key + ':' + value + '\n'

        self.tables_textEdit.setText(text.rstrip())

    def get_tables_from_text(self):
        text = self.tables_textEdit.toPlainText().rstrip()      # the rstrip removes trailing newlines etc
        lines = text.split("\n")
        table_dict = {}
        for line in lines:
            key_val = line.split(':')
            table_dict[key_val[0]] = key_val[1]
        return table_dict

    def reject(self):
        self.reset()
        super().reject()

class IPDatabaseConnectWidget(QFrame):
    def __init__(self, parent):
        super().__init__()
        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        size_policy = self.sizePolicy()
        size_policy.setHorizontalPolicy(QSizePolicy.Fixed)
        self.setSizePolicy(size_policy)
        self.parent = parent  # reference to the IPBeamformingWidget to which this belongs
        self.session = None
        self.config_filename = ""
        self.buildUI()
        
    def buildUI(self):
        
        self.title_label = QLabel("\tDatabase Connection")
        self.title_label.setStyleSheet("QLabel {font-weight: bold; color: white; background-color: black}")

        self.load_config_button = QPushButton("Load Config File...")
        self.load_config_button.setMaximumWidth(200)

        self.save_current_button = QPushButton("Save Config File...")
        self.save_current_button.setMaximumWidth(200)
        self.save_current_button.setEnabled(False)

        self.table_dialog = IPTableDialog(self)
        self.show_tables_button = QPushButton("Tables...")

        self.schema_type_combo = QComboBox()
        self.schema_type_combo.addItem("KBCore")
        self.schema_type_combo.addItem("CSS3")

        self.hostname_edit = QLineEdit()
        self.hostname_edit.setMaximumWidth(400)

        self.username_edit = QLineEdit()
        self.username_edit.setMaximumWidth(300)

        self.password_edit = QLineEdit()
        self.password_edit.setMaximumWidth(300)
        self.password_edit.setEchoMode(QLineEdit.PasswordEchoOnEdit)

        self.dialect_combo = QComboBox()
        self.dialect_combo.setMaximumWidth(100)
        self.dialect_combo.addItems(database.DIALECT_LIST)

        self.driver_edit = QLineEdit()
        self.driver_edit.setMaximumWidth(100)

        self.database_name = QLineEdit()
        self.database_name.setMaximumWidth(300)

        # I want to send an empty string if the portnum is not used, so I'm going to use a lineedit instead of a spinbox
        self.portnum_edit = QLineEdit()
        rxv = QRegExpValidator(QRegExp("\\d*"))
        self.portnum_edit.setValidator(rxv)
        self.portnum_edit.setMaximumWidth(100)

        self.url_edit = QLineEdit()
        self.url_edit.setPlaceholderText("URL")
        self.url_edit.setToolTip("You may edit this by hand if the url generated from the pieces doesn't work")

        self.create_session_button = QPushButton("Create Session")
        self.create_session_button.setMaximumWidth(200)

        self.close_session_button = QPushButton("Clear Session")
        self.close_session_button.setMaximumWidth(200)

        self.test_connection_button = QPushButton("Test Connection")

        horiz_layout_0 = QHBoxLayout()
        horiz_layout_0.addWidget(self.load_config_button)
        horiz_layout_0.addWidget(self.save_current_button)
        horiz_layout_0.addWidget(self.show_tables_button)
        horiz_layout_0.addStretch()
        
        form_layout = QFormLayout()
        form_layout.addRow("Schema: ", self.schema_type_combo)
        form_layout.addRow("Hostname: ", self.hostname_edit)
        form_layout.addRow("Username: ", self.username_edit)
        form_layout.addRow("Password: ", self.password_edit)
        form_layout.addRow("Dialect: ", self.dialect_combo)
        form_layout.addRow("Driver: ", self.driver_edit)
        form_layout.addRow("Database Name: ", self.database_name)
        form_layout.addRow("Port Number: ", self.portnum_edit)

        horiz_layout_1 = QHBoxLayout()
        horiz_layout_1.addStretch()
        horiz_layout_1.addWidget(self.create_session_button)
        horiz_layout_1.addWidget(self.close_session_button)
        horiz_layout_1.addWidget(self.test_connection_button)
        horiz_layout_1.addStretch()

        horiz_layout_2 = QHBoxLayout()
        horiz_layout_2.addLayout(form_layout)
        horiz_layout_2.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addWidget(self.title_label)
        main_layout.addLayout(horiz_layout_0)
        main_layout.addLayout(horiz_layout_2)
        main_layout.addWidget(self.url_edit)
        main_layout.addStretch()
        main_layout.addLayout(horiz_layout_1)
        self.setLayout(main_layout)

        self.config_file_dialog = QFileDialog()
        self.config_file_dialog.setFileMode(QFileDialog.ExistingFile)
        self.config_file_dialog.setNameFilter("(*.ini)")

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.hostname_edit.textEdited.connect(self.update_url)
        self.username_edit.textEdited.connect(self.update_url)
        self.password_edit.textEdited.connect(self.update_url)
        self.dialect_combo.currentIndexChanged.connect(self.update_url)
        self.schema_type_combo.currentIndexChanged.connect(self.update_url)

        self.url_edit.textEdited.connect(self.url_manually_edited)

        self.driver_edit.textEdited.connect(self.update_url)
        self.database_name.textEdited.connect(self.update_url)
        self.portnum_edit.textEdited.connect(self.update_url)

        self.load_config_button.pressed.connect(self.load_config_file)
        self.save_current_button.pressed.connect(self.save_current_config)
        self.show_tables_button.clicked.connect(self.show_tables_dialog)
        self.create_session_button.pressed.connect(self.create_session)
        self.close_session_button.pressed.connect(self.close_session)

        self.test_connection_button.pressed.connect(self.check_connection)

    def show_tables_dialog(self):
        if self.table_dialog.exec_():
            self.save_current_button.setEnabled(True)

    def update_url(self):
        
        dialect = self.dialect_combo.currentText()
        driver = self.driver_edit.text()
        if driver:
            driver = '+' + driver
        username = self.username_edit.text()
        # special characters in the password need to be correctly parsed into url strings
        password = self.password_edit.text()
        hostname = self.hostname_edit.text()
        port = self.portnum_edit.text()
        db_name = self.database_name.text()

        self.url_edit.setText(dialect + driver + "://" + username + ":" + password + "@" + hostname + ":" + port + "/" + db_name)

        # if the url changes... so has some of the information.  This might need to be saved, so activate the save_current button
        self.save_current_button.setEnabled(True)

    def url_manually_edited(self):
        self.save_current_button.setEnabled(True)

    def create_session(self):
        # first, if there is already an active session, close it...
        self.close_session()

        # now proceed as usual
        dialect = self.dialect_combo.currentText()
        driver = self.driver_edit.text()

        username = self.username_edit.text()
        # special characters in the password need to be correctly parsed into url strings
        password = urllib.parse.quote_plus(self.password_edit.text())
        hostname = self.hostname_edit.text()
        port = self.portnum_edit.text()
        db_name = self.database_name.text()

        url = self.url_edit.text()

        try:
            #self.session = database.db_connect_url(url)
            self.session = database.db_connect(dialect=dialect, 
                                               hostname=hostname, 
                                               db_name=db_name, 
                                               port=port, 
                                               username=username, 
                                               password=password, 
                                               driver=driver)
            self.url_edit.setStyleSheet("color: green")
        except Exception as e:
            self.url_edit.setStyleSheet("color: red")
            self.errorPopup("Error connecting to url:\n{}\n{}".format(url, e))
    
    def close_session(self):
        if self.session is not None:
            self.session.close()
            self.url_edit.setStyleSheet("color: black")
            self.session = None

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

    @pyqtSlot()
    def reset_connection_colors(self):
        self.test_connection_button.setText('Test Connection')
        self.test_connection_button.setStyleSheet('QPushButton {color: black}')


    def load_config_file(self):
        if self.config_file_dialog.exec_():
            self.config_filename = self.config_file_dialog.selectedFiles()[0]
            try:
                config = configparser.ConfigParser()
                config.read(self.config_filename)
                self.schema_type_combo.setCurrentText(config['DATABASE']['schema'])
                self.hostname_edit.setText(config['DATABASE']['hostname'])
                self.username_edit.setText(config['DATABASE']['username'])
                self.password_edit.setText("")
                self.database_name.setText(config['DATABASE']['database_name'])
                self.portnum_edit.setText(config['DATABASE']['port'])
                self.driver_edit.setText(config['DATABASE']['driver'])
                self.dialect_combo.setCurrentText(config['DATABASE']['dialect'])
                if config.has_option('DATABASE', 'url'):
                    self.url_edit.setText(config['DATABASE']['url'])
                else:
                    my_url = database.assemble_db_url(self.dialect_combo.currentText(), 
                                                      self.hostname_edit.text(), 
                                                      self.database_name.text(), 
                                                      self.portnum_edit.text(), 
                                                      self.username_edit.text(),
                                                      self.password_edit.text(),
                                                      self.driver_edit.text())
                    self.url_edit.setText(my_url)

                self.table_dialog.set_text_from_table_dict(config['DBTABLES'])

                self.save_current_button.setEnabled(False)

            except Exception as e:
                self.errorPopup("Error reading config file \n{}".format(str(e)))

    def save_current_config(self):
        save_filename = QFileDialog.getSaveFileName(self, "Save db configuration", "/home/jwebster/IPProjects", "(*.ini)")[0]
        if save_filename:
            try:
                config = configparser.ConfigParser()
                config['DATABASE'] = {}
                config['DATABASE']['schema'] = self.schema_type_combo.currentText()
                config['DATABASE']['hostname'] = self.hostname_edit.text()
                config['DATABASE']['username'] = self.username_edit.text()
                config['DATABASE']['database_name'] = self.database_name.text()
                config['DATABASE']['port'] = self.portnum_edit.text()
                config['DATABASE']['driver'] = self.driver_edit.text()
                config['DATABASE']['dialect'] = self.dialect_combo.currentText()
                config['DATABASE']['url'] = self.url_edit.text()
                
                table_dict = self.table_dialog.get_tables_from_text()
                config['DBTABLES'] = table_dict

                with open(save_filename, 'w') as config_file:
                    config.write(config_file)

                self.save_current_button.setEnabled(False)
            except Exception as e:
                self.errorPopup("Error saving config file. \n{}".format(e))

    @pyqtSlot(str, str)
    def errorPopup(self, message, title="Oops..."):
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setText(message)
        msg_box.setWindowTitle(title)
        msg_box.exec_()