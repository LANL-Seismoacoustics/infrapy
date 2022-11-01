#!/usr/bin/env python

from PyQt5 import QtWidgets
from PyQt5.QtGui import QIcon

import os
import sys
import platform

from InfraView.widgets.IPApplicationWindow import *


def main():
    progname = "InfraView"
    progversion = "0.4.0"

    here = os.path.abspath(__file__)
    here = os.path.dirname(here)

    qApp = QtWidgets.QApplication(sys.argv)

    #icon_file = here + '/graphics/infraViewIcon.svg'
    icon_file = here + '/graphics/icons/start_64.png'
    qApp.setWindowIcon(QIcon(icon_file))

    my_system = platform.system()
    my_release = platform.release()

    aw = IPApplicationWindow(qApp, progname, progversion)

    aw.setWindowTitle("%s" % progname)
    aw.show()

    sys.exit(qApp.exec())


if __name__ == '__main__':
    main()
