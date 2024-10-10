from PyQt5.QtWidgets import QMessageBox, QSplitter
from PyQt5.QtGui import QValidator
from PyQt5.QtCore import Qt, pyqtSlot

import pyqtgraph as pg


    

# Define some useful colors here
reb_blue = pg.mkColor(80, 159, 250)
reb_red = pg.mkColor(255, 71, 71)   

lanl_primary = pg.mkColor("#000F7E")
lanl_gradient = pg.mkColor("#090238")
lanl_secondary = pg.mkColor("#0070C1")
lanl_secondary_tint = pg.mkColor("#3296DC")
lanl_screen_text_black = pg.mkColor("#0C0D17")
lanl_dark_grey = pg.mkColor("#555962")
lanl_light_grey = pg.mkColor("#CDD1E2")
lanl_background_accent = pg.mkColor("#F1EFF7")

lanl_blue = pg.mkColor("#0070C1")
lanl_blue_tint = pg.mkColor("#3296DC")
lanl_red = pg.mkColor("#EB0F1E")
lanl_red_tint = pg.mkColor("#FF474D")
lanl_orange = pg.mkColor("#E17800")
lanl_orange_tint = pg.mkColor("#FF9129")
lanl_green = pg.mkColor("#00AA64")
lanl_green_tint = pg.mkColor("#2CC486")

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

class IPSplitter(QSplitter):

    def __init__(self, orientation, parent=None):
        super().__init__(orientation, parent)

        self.setStyleSheet("QSplitter::handle{ background-color: #DDD}")
        self.setHandleWidth(20)


class CapsValidator(QValidator):
    ''' 
    since many text fields require capitalized values only, here is a validator for the lineEdits etc
    general usage is something like...
    my_validator = IPUtils.CapsValidator()
    my_lineedit.setValidator(my_validator)
    '''

    def validate(self, string, pos):
        return QValidator.Acceptable, string.upper(), pos
