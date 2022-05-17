
from PyQt5.QtWidgets import QGridLayout, QWidget

from pyqtgraph.parametertree import Parameter, ParameterTree


class IPDisplaySettingsWidget(QWidget):

    params = [
        {'name': 'Axes Settings', 'type': 'group', 'children': [
            {'name': 'X-Axis Labels: ', 'type': 'list', 'values': {"Show All": 1,
                                                                   "Show only on bottom plot": 2, "Hide All": 3}, 'value': 2},
            {'name': 'X-Axis Values: ', 'type': 'list', 'values': {"Show All": 1,
                                                                   "Show only on bottom plot": 2, "Hide All": 3}, 'value': 2},
            {'name': 'Y-Axis Labels: ', 'type': 'list', 'values': {"Show All": 1, "Hide All": 2}, 'value': 2},
            {'name': 'Y-Axis Values: ', 'type': 'list', 'values': {"Show All": 1, "Hide All": 2}, 'value': 2}
        ], 'expanded': False},

        {'name': 'Beamforming Tab', 'type': 'group', 'children': [
            {'name': 'Show Grid on Results Plots: ', 'type': 'bool', 'value': True}]}

        # {'name': 'Pick Lines', 'type': 'group', 'children':[
        #     {'name': 'Default Color:', 'type': 'color', 'value': "F00", 'tip': "This is the default color for pick lines"},
        # ], 'expanded': False},

        # {'name': 'Crosshairs', 'type': 'group', 'children':[
        #     {'name': 'Mirror on all plots: ', 'type': 'bool', 'value': True},
        #     {'name': 'Show Horizontal Marker: ', 'type': 'bool', 'value': True},
        #     {'name': 'Show Vertical Marker: ', 'type': 'bool', 'value': True}

        # ], 'expanded': False},

        # {'name': 'Blah Balh', 'type': 'group', 'children':[
        #     {'name': 'String', 'type': 'str', 'value': "hi"},
        #     {'name': 'List', 'type': 'list', 'values': [1,2,3], 'value': 2},
        #     {'name': 'Named List', 'type': 'list', 'values': {"one": 1, "two": "twosies", "three": [3,3,3]}, 'value': 2},
        #     {'name': 'Boolean', 'type': 'bool', 'value': True, 'tip': "This is a checkbox"},
        #     {'name': 'Default Color', 'type': 'color', 'value': "FF0", 'tip': "This is a color button"},
        #     {'name': 'Gradient', 'type': 'colormap'},
        #     {'name': 'Subgroup', 'type': 'group', 'children': [
        #         {'name': 'Sub-param 1', 'type': 'int', 'value': 10},
        #         {'name': 'Sub-param 2', 'type': 'float', 'value': 1.2e6},
        #     ]},
        #     {'name': 'Text Parameter', 'type': 'text', 'value': 'Some text...'},
        #     {'name': 'Action Parameter', 'type': 'action'},
        # ], 'expanded': False},

        # {'name': 'Numerical Parameter Options', 'type': 'group', 'children': [
        #     {'name': 'Units + SI prefix', 'type': 'float', 'value': 1.2e-6, 'step': 1e-6, 'siPrefix': True, 'suffix': 'V'},
        #     {'name': 'Limits (min=7;max=15)', 'type': 'int', 'value': 11, 'limits': (7, 15), 'default': -6},
        #     {'name': 'DEC stepping', 'type': 'float', 'value': 1.2e6, 'dec': True, 'step': 1, 'siPrefix': True, 'suffix': 'Hz'},

        # ], 'expanded': False},

        # {'name': 'Save/Restore functionality', 'type': 'group', 'children': [
        #     {'name': 'Save State', 'type': 'action'},
        #     {'name': 'Restore State', 'type': 'action', 'children': [
        #         {'name': 'Add missing items', 'type': 'bool', 'value': True},
        #         {'name': 'Remove extra items', 'type': 'bool', 'value': True},
        #     ]},
        # ], 'expanded': False},

        # {'name': 'Extra Parameter Options', 'type': 'group', 'children': [
        #     {'name': 'Read-only', 'type': 'float', 'value': 1.2e6, 'siPrefix': True, 'suffix': 'Hz', 'readonly': True},
        #     {'name': 'Renamable', 'type': 'float', 'value': 1.2e6, 'siPrefix': True, 'suffix': 'Hz', 'renamable': True},
        #     {'name': 'Removable', 'type': 'float', 'value': 1.2e6, 'siPrefix': True, 'suffix': 'Hz', 'removable': True},
        # ], 'expanded': False},
    ]

    def __init__(self):
        super().__init__()
        self.__buildUI__()

    def __buildUI__(self):

        mainLayout = QGridLayout()
        self.settingsTree1 = ParameterTree(showHeader=True)
        self.p = Parameter.create(name='params', type='group', children=self.params)
        self.settingsTree1.setParameters(self.p, showTop=False)
        mainLayout.addWidget(self.settingsTree1)
        self.setLayout(mainLayout)

        # for child in self.p.children():
        #    child.sigValueChanging.connect(self.valueChanging)
        #    for ch2 in child.children():
        #        ch2.sigValueChanging.connect(self.valueChanging)

        beamformingSettings = self.p.child('Beamforming Tab')
        # for child in beamformingSettings.children():
        #    print(child)

    def printParams(self):
        print(self.p.child('Axes Settings'))

    def valueChanging(self, param, value):
        print("Value changing (not finalized): {} {}".format(param, value))

    def getP(self):
        return self.p
