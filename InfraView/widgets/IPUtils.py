from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtGui import QValidator
from PyQt5.QtCore import pyqtSlot    

import pyqtgraph as pg
    

reb_blue = pg.mkColor(80, 159, 250)
reb_red = pg.mkColor(255, 71, 71)   

# arrays of QColors
# blue to red... one white color removed
blue_to_red = [pg.mkColor("#1984c5"), pg.mkColor("#c23728"), pg.mkColor("#22a7f0"), pg.mkColor("#63bff0"), pg.mkColor("#a7d5ed"), pg.mkColor("#e1a692"), pg.mkColor("#de6e56"), pg.mkColor("#e14b31")]

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
