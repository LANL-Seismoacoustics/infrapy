from PyQt5.QtCore import QLine


class IPLine(QLine):

    my_zValue = 20

    def __init__(self):
        super().__init__()
        self.setZValue(self.my_zValue)
