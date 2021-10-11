from PyQt5.QtWidgets import (QWidget, QFormLayout, QLineEdit, QPushButton, QVBoxLayout)


class IPDatabaseConnectWidget(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent  # reference to the IPBeamformingWidget to which this belongs
        self.buildUI()

    def buildUI(self):
        form_layout = QFormLayout()
        url_edit = QLineEdit()
        
        username_edit = QLineEdit()
        password_edit = QLineEdit()
        password_edit.setEchoMode(QLineEdit.Password)

        form_layout.addRow("Database URL: ", url_edit)
        form_layout.addRow("Username: ", username_edit)
        form_layout.addRow("Password: ", password_edit)

        connect_button = QPushButton("Connect")
        
        main_layout = QVBoxLayout()
        main_layout.addLayout(form_layout)
        main_layout.addWidget(connect_button)
        self.setLayout(main_layout)

    def connect_to_database(self):
        pass