import urllib.parse
import configparser

from PyQt5.QtWidgets import (QWidget, QComboBox, QFileDialog, QFormLayout, QFrame, QHBoxLayout, QLabel, QLineEdit, QMessageBox, 
                             QPushButton, QVBoxLayout)
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtCore import QRegExp, pyqtSlot, QTimer
from infrapy.utils import database


class IPDatabaseConnectWidget(QFrame):
    def __init__(self, parent):
        super().__init__()
        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        self.parent = parent  # reference to the IPBeamformingWidget to which this belongs
        self.session = None
        self.buildUI()

    def buildUI(self):
        self.form_layout = QFormLayout()

        self.title_label = QLabel("\tDatabase Connection")
        self.title_label.setStyleSheet("QLabel {font-weight: bold; color: white; background-color: black}")

        self.load_config_button = QPushButton("Load Config File...")
        self.load_config_button.setMaximumWidth(200)

        self.save_current_button = QPushButton("Save Config File...")
        self.save_current_button.setMaximumWidth(200)
        self.save_current_button.setEnabled(False)

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

        self.url_label = QLabel("")

        self.create_session_button = QPushButton("Create Session")
        self.create_session_button.setMaximumWidth(200)

        self.close_session_button = QPushButton("Close Session")
        self.close_session_button.setMaximumWidth(200)

        self.test_connection_button = QPushButton("Test Connection")
        self.test_connection_label = QLabel("")

        self.form_layout.addWidget(self.load_config_button)
        self.form_layout.addWidget(self.save_current_button)
        self.form_layout.addRow("Hostname: ", self.hostname_edit)
        self.form_layout.addRow("Username: ", self.username_edit)
        self.form_layout.addRow("Password: ", self.password_edit)
        self.form_layout.addRow("Dialect: ", self.dialect_combo)
        self.form_layout.addRow("Driver: ", self.driver_edit)
        self.form_layout.addRow("Database Name: ", self.database_name)
        self.form_layout.addRow("Port Number: ", self.portnum_edit)
        self.form_layout.addRow("URL: ", self.url_label)

        horiz_layout_1 = QHBoxLayout()
        horiz_layout_1.addStretch()
        horiz_layout_1.addWidget(self.create_session_button)
        horiz_layout_1.addWidget(self.close_session_button)
        horiz_layout_1.addWidget(self.test_connection_button)
        horiz_layout_1.addWidget(self.test_connection_label)
        horiz_layout_1.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addWidget(self.title_label)
        main_layout.addLayout(self.form_layout)
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

        self.driver_edit.textEdited.connect(self.update_url)
        self.database_name.textEdited.connect(self.update_url)
        self.portnum_edit.textEdited.connect(self.update_url)

        self.load_config_button.pressed.connect(self.load_config_file)
        self.save_current_button.pressed.connect(self.save_current_config)
        self.create_session_button.pressed.connect(self.create_session)
        self.close_session_button.pressed.connect(self.close_session)

        self.test_connection_button.pressed.connect(self.check_connection)

    def update_url(self):
        self.url_label.setStyleSheet("QLabel {color:black}")
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

        self.url_label.setText(dialect + driver + "://" + username + ":" + password + "@" + hostname + ":" + port + "/" + db_name)

        # if the url changes... so has some of the information.  This might need to be saved, so activate the save_current button
        self.save_current_button.setEnabled(True)

    def create_session(self):
        # first, if there is already an active session, close it...
        self.close_session()

        # now proceed as usual
        dialect = self.dialect_combo.currentText()
        driver = self.driver_edit.text()
        if driver:
            driver = '+' + driver
        username = self.username_edit.text()
        # special characters in the password need to be correctly parsed into url strings
        password = urllib.parse.quote_plus(self.password_edit.text())
        hostname = self.hostname_edit.text()
        port = self.portnum_edit.text()
        db_name = self.database_name.text()

        url = dialect + driver + "://" + username + ":" + password + "@" + hostname + ":" + port + "/" + db_name

        try:
            self.session = database.db_connect_url(url)
            self.url_label.setStyleSheet("QLabel {color: green}")
        except Exception as e:
            self.url_label.setStyleSheet("QLabel {color: red}")
            self.errorPopup("Error connecting to url:\n{}\n{}".format(url, e))
    
    def close_session(self):
        if self.session is not None:
            self.session.close()
            self.url_label.setStyleSheet("QLabel {color: black}")
            self.session = None

    def check_connection(self):
        if self.session is not None:
            if database.check_connection(self.session):
                self.test_connection_label.setStyleSheet("QLabel {color: green}")
                self.test_connection_label.setText("Connection is good")
                QTimer.singleShot(3000, self.reset_connection_colors)
            else:
                self.test_connection_label.setStyleSheet("QLabel {color: red}")
                self.test_connection_label.setText("Bad Connection")
                QTimer.singleShot(3000, self.reset_connection_colors)
        else:
            self.test_connection_label.setText("No active session")
            QTimer.singleShot(3000, self.reset_connection_colors)

    @pyqtSlot()
    def reset_connection_colors(self):
        self.test_connection_label.setStyleSheet("QLabel {color: black}")
        self.test_connection_label.setText("")

    def load_config_file(self):
        if self.config_file_dialog.exec_():
            config_filename = self.config_file_dialog.selectedFiles()[0]
            try:
                config = configparser.ConfigParser()
                config.read(config_filename)
                self.hostname_edit.setText(config['DATABASE']['hostname'])
                self.username_edit.setText(config['DATABASE']['username'])
                self.database_name.setText(config['DATABASE']['database_name'])
                self.portnum_edit.setText(config['DATABASE']['port'])
                self.driver_edit.setText(config['DATABASE']['driver'])
                self.dialect_combo.setCurrentText(config['DATABASE']['dialect'])
                self.update_url()
                self.save_current_button.setEnabled(False)
            except Exception as e:
                self.errorPopup("Error reading config file \n{}".format(e))

    def save_current_config(self):
        save_filename = QFileDialog.getSaveFileName(self, "Save db configuration", "/home/jwebster/IPProjects", "(*.ini)")[0]
        if save_filename:
            try:
                config = configparser.ConfigParser()
                config['DATABASE'] = {}
                config['DATABASE']['hostname'] = self.hostname_edit.text()
                config['DATABASE']['username'] = self.username_edit.text()
                config['DATABASE']['database_name'] = self.database_name.text()
                config['DATABASE']['port'] = self.portnum_edit.text()
                config['DATABASE']['driver'] = self.driver_edit.text()
                config['DATABASE']['dialect'] = self.dialect_combo.currentText()
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