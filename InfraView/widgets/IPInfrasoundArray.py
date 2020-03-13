class IPInfrasoundArray(object):

    def __init(self, inventory=None, a_name=''):
        self._array_name = a_name
        self._inv = inventory

    def set_name(self, name):
        # set the arrays name
        self._array_name = name

    def name(self):
        # return the arrays name
        return self._array_name

    def set_inventory(self, inventory):
        self._inv = inventory

    def station_center(self):
        # This will return the lat,lon center of the stations in the array
        center_lat = None
        center_lon = None

        count = 0   # The number of stations in the array

        for network in self._inv:
            for station in network:
                if center_lat is None:
                    center_lat = station.latitude
                else:
                    center_lat += station.latitude

                if center_lon is None:
                    center_lon = station.longitude
                else:
                    center_lon += station.longitude

                count += 1

        center_lat = center_lat / count
        center_lon = center_lon / count

        return center_lat, center_lon

    def avg_elevation(self):
        # This will calculate and return the average elevation of the sensors in the array

        avg_elevation = None

        count = 0   # The number of stations in the array

        for network in self._inv:
            for station in network:
                if avg_elevation is None:
                    avg_elevation = station.elevation
                else:
                    avg_elevation += station.elevation

                count += 1

        if avg_elevation is not None:
            return avg_elevation / count
        else:
            return None
