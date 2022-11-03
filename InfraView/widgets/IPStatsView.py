from PyQt5.QtWidgets import QTextEdit, QTabWidget
from PyQt5.QtCore import pyqtSignal, pyqtSlot


class IPStatsView(QTabWidget):

    removeTrace = pyqtSignal(str)

    # class to ease the process of displaying the obspy trace stats
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setTabsClosable(True)
        self.tabCloseRequested.connect(self.removeTraceByName)
        self.show()

    def setStats(self, streamList):

        self.clear()

        for idx, trace in enumerate(streamList):
            newTab = QTextEdit()
            self.addTab(newTab, trace.id)

            # build the string to show...
            statString = ''
            for key, value in trace.stats.items():
                if str(key) == 'mseed':
                    statString = statString + '<b>MSEED info:</b>'
                    for mkey, mvalue in value.items():
                        statString = statString + '<li>' + str(mkey) + ': ' + str(mvalue) + '</li>'
                    statString = statString + '</ul>'
                elif str(key) == 'sac':
                    statString = statString + '<b>SAC info:</b>'
                    statString = statString + '<ul>'
                    for mkey, mvalue in value.items():
                        statString = statString + '<li>' + str(mkey) + ': ' + str(mvalue) + '</li>'
                    statString = statString + '</ul>'
                else:
                    statString = statString + '<b>' + str(key) + '</b>: ' + str(value) + '<br>'

                # add the newly built string to the view
                # First get the QTextEdit in the current tab
                textEdit = self.widget(idx)
                textEdit.setHtml(statString)

        return

    @pyqtSlot(int)
    def removeTraceByName(self, idx):
        title = self.tabText(idx)
        self.removeTrace.emit(title)

