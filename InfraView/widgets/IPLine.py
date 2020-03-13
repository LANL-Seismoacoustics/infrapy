from PyQt5.QtCore import QLine


class IPLine(QLine):

    my_zValue = 5

    def __init__(self):
        super().__init__()

    def zValue(self):
        return self.my_zValue

    def setZValue(self, z):
        self.zValue = z
