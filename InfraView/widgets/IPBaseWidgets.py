from PyQt5.QtWidgets import QWidget, QSizePolicy

class IPSettingsWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setSizePolicy(QSizePolicy.Maximum,
                           QSizePolicy.Maximum)
        font = self.font()
        fontsize = font.pointSize()
        if fontsize > 10:
            font.setPointSize(fontsize-2)
        # font.setFamily("monospace")
        self.setFont(font)

