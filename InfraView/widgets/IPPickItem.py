from infrapy.propagation.likelihoods import InfrasoundDetection


class IPPickItem(InfrasoundDetection):

    _associatedPickLine = None    # reference to the pick line associated with this pick

    def __init__(self, name=""):
        super().__init__()

        self.set_name(name)

    def getAssociatedPickLine(self):
        return self._associatedPickLine

    def setAssociatedPickLine(self, pickline):
        self._associatedPickLine = pickline

    def to_InfrasoundDetection(self):
        return InfrasoundDetection(lat_loc=self.latitude, lon_loc=self.longitude, time=self.peakF_UTCtime,
                                   azimuth=self.back_azimuth, f_stat=self.peakF_value, array_d=self.array_dim)
