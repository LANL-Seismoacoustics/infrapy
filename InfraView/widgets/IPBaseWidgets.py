from PyQt5.QtWidgets import QWidget

class IPSettingsWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        font = self.font()
        font.setPointSize(11)
        self.setFont(font)

