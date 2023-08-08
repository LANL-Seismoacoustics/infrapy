import os
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl
from PyQt5.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QPushButton, QLabel, QToolBar, QLineEdit, QAction



class IPDocViewer(QWidget):
    def __init__(self, parent):
        super().__init__(parent)

        self.webView = QWebEngineView()
        self.webView.urlChanged.connect(self.on_url_changed)

        toolbar = self.create_toolbar()

        app_dir = parent.get_app_directory()
        base_dir = os.path.abspath(os.path.join(app_dir, os.pardir))
        self.doc_dir = base_dir + '/docs/'
        self.html_home_dir = self.doc_dir + 'build/html/'
        self.index_file = self.html_home_dir + 'index.html'
        os.makedirs(self.html_home_dir, exist_ok=True)

        self.compileLabel = QLabel("It appears you haven't compiled the documentation yet. Would you like to do it now?")
        self.compileButton = QPushButton("Compile Documentation")
        self.compileButton.setFlat(True)
        self.compileButton.clicked.connect(self.compile_documentation)

        if self.check_for_index_file():
            self.compileLabel.setVisible(False)
            self.compileButton.setVisible(False)
            self.load_html_file(self.index_file)
        else:
            self.compileLabel.setVisible(True)
            self.compileButton.setVisible(True)

        hLayout = QHBoxLayout()
        hLayout.addWidget(self.compileLabel)
        hLayout.addWidget(self.compileButton)
        hLayout.addStretch()

        vLayout = QVBoxLayout()
        vLayout.addWidget(toolbar)
        vLayout.addLayout(hLayout)
        vLayout.addWidget(self.webView)

        self.setLayout(vLayout)

    def create_toolbar(self):
        toolbar = QToolBar()

        back_action = QAction("Back", self)
        back_action.triggered.connect(self.webView.back)
        toolbar.addAction(back_action)

        forward_action = QAction("Forward", self)
        forward_action.triggered.connect(self.webView.forward)
        toolbar.addAction(forward_action)

        reload_action = QAction()
        reload_action.triggered.connect(self.webView.reload)
        toolbar.addAction(reload_action)

        self.search_bar = QLineEdit()
        self.search_bar.returnPressed.connect(self.load_url)
        toolbar.addWidget(self.search_bar)

        return toolbar

    def check_for_index_file(self):
        return os.path.isfile(self.html_home_dir + 'index.html')

    def on_url_changed(self, url):
        print("new url: {}".format(url))

    def compile_documentation(self):
        # run the make file and hope for the best
        command = 'make -C ' + self.doc_dir + ' html'
        print(command)
        os.system(command)
        self.load_html_file(self.index_file)

        if self.check_for_index_file():
            self.compileLabel.setVisible(False)
            self.compileButton.setVisible(False)
            self.load_html_file(self.index_file)
        else:
            self.compileLabel.setVisible(True)
            self.compileButton.setVisible(True)

    def load_url(self):
        url = self.search_bar.text()
        if not url.startswith("http"):
            url = "https://" + url
        self.webView.load(QUrl(url))

    def load_html_file(self, file):
        with open(file, 'r') as f:
            html = f.read()
            self.webView.setHtml(html)

