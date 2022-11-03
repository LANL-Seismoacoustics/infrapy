from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtGui import QValidator
from PyQt5.QtCore import pyqtSlot    
    
    
@pyqtSlot(str, str)
def errorPopup(message, title="Oops..."):
    # generic error popup dialog for generic errors
    title = "InfraView: " + title 
    msgBox = QMessageBox()
    msgBox.setIcon(QMessageBox.Warning)
    msgBox.setText(message)
    msgBox.setWindowTitle(title)
    msgBox.exec_()


class CapsValidator(QValidator):
    # since many text fields require capitalized values only, here is a validator for the lineEdits etc
    # general usage is something like...
    #     my_validator = IPUtils.CapsValidator()
    #     my_lineedit.setValidator(my_validator)

    def validate(self, string, pos):
        return QValidator.Acceptable, string.upper(), pos
