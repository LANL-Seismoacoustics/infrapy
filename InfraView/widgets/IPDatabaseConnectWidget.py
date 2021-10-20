from PyQt5.QtWidgets import (QWidget, QComboBox, QFormLayout, QLabel, QLineEdit, QMessageBox, 
                             QPushButton, QSpinBox, QVBoxLayout)
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtCore import QRegExp, pyqtSignal, pyqtSlot
from infrapy.utils import database

import urllib.parse


class IPDatabaseConnectWidget(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent  # reference to the IPBeamformingWidget to which this belongs
        self.buildUI()

    def buildUI(self):
        self.form_layout = QFormLayout()

        self.load_config_button = QPushButton("Load Configuration File")
        self.load_config_button.setMaximumWidth(200)

        self.save_current_button = QPushButton("Save Current")
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
        self.dialect_combo.addItems(database.dialect_list)

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

        self.connect_button = QPushButton("Connect")
        self.connect_button.setMaximumWidth(200)

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
        self.form_layout.addWidget(self.connect_button)

        main_layout = QVBoxLayout()
        main_layout.addLayout(self.form_layout)
        self.setLayout(main_layout)

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.hostname_edit.textEdited.connect(self.wake_up_save_button)
        self.username_edit.textEdited.connect(self.wake_up_save_button)
        self.password_edit.textEdited.connect(self.wake_up_save_button)
        self.dialect_combo.currentIndexChanged.connect(self.wake_up_save_button)
        self.driver_edit.textEdited.connect(self.wake_up_save_button)
        self.database_name.textEdited.connect(self.wake_up_save_button)
        self.portnum_edit.textEdited.connect(self.wake_up_save_button)

        self.hostname_edit.textEdited.connect(self.update_url)
        self.username_edit.textEdited.connect(self.update_url)
        self.password_edit.textEdited.connect(self.update_url)
        self.dialect_combo.currentIndexChanged.connect(self.update_url)
        self.driver_edit.textEdited.connect(self.update_url)
        self.database_name.textEdited.connect(self.update_url)
        self.portnum_edit.textEdited.connect(self.update_url)

        self.connect_button.pressed.connect(self.connect_to_database)

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

        self.url_label.setText(dialect + driver + "://" + username + ":" + password + "@" + hostname + ":" + port + "/" + db_name)

    def wake_up_save_button(self):
        self.save_current_button.setEnabled(True)

    def connect_to_database(self):
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
        except ModuleNotFoundError as err:
            self.errorPopup("Missing module... Infrapy doesn't automatically install database modules and drivers into the infrapy_env environment, so that will to be done manually\n\n {}".format(err))

        if self.test_connection():
            print("good connection")
        else:
            print("bad connection")

    def test_connection(self):
        try:
            my_engine = self.session.get_bind()
            conn = my_engine.connect()
            conn.close()
            return True
        except:
            return False

    @pyqtSlot(str, str)
    def errorPopup(self, message, title="Oops..."):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setText(message)
        msgBox.setWindowTitle(title)
        msgBox.exec_()