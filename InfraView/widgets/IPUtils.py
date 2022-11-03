from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtGui import QValidator
from PyQt5.QtCore import pyqtSlot    
    
    
@pyqtSlot(str, str)
def errorPopup(message, title="Oops..."):
    title = "InfraView: " + title 
    msgBox = QMessageBox()
    msgBox.setIcon(QMessageBox.Warning)
    msgBox.setText(message)
    msgBox.setWindowTitle(title)
    msgBox.exec_()

class CapsValidator(QValidator):
    # since most of the fields will require capitalized values only, here is a validator for the
    # lineEdits
    def validate(self, string, pos):
        return QValidator.Acceptable, string.upper(), pos