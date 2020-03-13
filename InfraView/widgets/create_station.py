import obspy
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.clients.nrl import NRL


# We'll first create all the various objects. These strongly follow the
# hierarchy of StationXML files.
inv = Inventory(
    # We'll add networks later.
    networks=[],
    # The source should be the id whoever create the file.
    source="LANL_Webster")

net = Network(
    # This is the network code according to the SEED standard.
    code="LN",
    # A list of stations. We'll add one later.
    stations=[],
    description="Los Alamos.",
    # Start-and end dates are optional.
    start_date=obspy.UTCDateTime(2000, 1, 1))

sta = Station(
    # This is the station code according to the SEED standard.
    code="STGN",
    latitude=37.016061502,
    longitude=-113.616194995,
    elevation=818,
    creation_date=obspy.UTCDateTime(2018, 12, 12),
    site=Site(name="Saint George North"))

cha = Channel(
    # This is the channel code according to the SEED standard.
    code="BDF",
    # This is the location code according to the SEED standard.
    location_code="00",
    # Note that these coordinates can differ from the station coordinates.
    latitude=37.016061502,
    longitude=-113.616194995,
    elevation=818,
    depth=0.0,
    azimuth=0.0,
    dip=-90.0,
    sample_rate=200)

# By default this accesses the NRL online. Offline copies of the NRL can
# also be used instead
#nrl = NRL()
# The contents of the NRL can be explored interactively in a Python prompt,
# see API documentation of NRL submodule:
# http://docs.obspy.org/packages/obspy.clients.nrl.html
# Here we assume that the end point of data logger and sensor are already
# known:
# response = nrl.get_response( # doctest: +SKIP
#     sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
#     datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])


# # Now tie it all together.
# cha.response = response
sta.channels.append(cha)
net.stations.append(sta)


sta = Station(
    # This is the station code according to the SEED standard.
    code="STGS",
    latitude=37.015069043,
    longitude=-113.616940566,
    elevation=817,
    creation_date=obspy.UTCDateTime(2018, 12, 12),
    site=Site(name="Saint George South"))

cha = Channel(
    # This is the channel code according to the SEED standard.
    code="BDF",
    # This is the location code according to the SEED standard.
    location_code="00",
    # Note that these coordinates can differ from the station coordinates.
    latitude=37.015069043,
    longitude=-113.616940566,
    elevation=817,
    depth=0.0,
    azimuth=0.0,
    dip=-90.0,
    sample_rate=200)

sta.channels.append(cha)
net.stations.append(sta)

sta = Station(
    # This is the station code according to the SEED standard.
    code="STGE",
    latitude=37.015146,
    longitude=-113.616152,
    elevation=819,
    creation_date=obspy.UTCDateTime(2018, 12, 12),
    site=Site(name="Saint George East"))

cha = Channel(
    # This is the channel code according to the SEED standard.
    code="BDF",
    # This is the location code according to the SEED standard.
    location_code="00",
    # Note that these coordinates can differ from the station coordinates.
    latitude=37.015146,
    longitude=-113.616152,
    elevation=819,
    depth=0.0,
    azimuth=0.0,
    dip=-90.0,
    sample_rate=200)

sta.channels.append(cha)
net.stations.append(sta)

sta = Station(
    # This is the station code according to the SEED standard.
    code="STGW",
    latitude=37.015581303,
    longitude=-113.617117339,
    elevation=817,
    creation_date=obspy.UTCDateTime(2018, 12, 12),
    site=Site(name="Saint George West"))

cha = Channel(
    # This is the channel code according to the SEED standard.
    code="BDF",
    # This is the location code according to the SEED standard.
    location_code="00",
    # Note that these coordinates can differ from the station coordinates.
    latitude=37.015581303,
    longitude=-113.617117339,
    elevation=817,
    depth=0.0,
    azimuth=0.0,
    dip=-90.0,
    sample_rate=200)

sta.channels.append(cha)
net.stations.append(sta)


inv.networks.append(net)

# And finally write it to a StationXML file. We also force a validation against
# the StationXML schema to ensure it produces a valid StationXML file.
#
# Note that it is also possible to serialize to any of the other inventory
# output formats ObsPy supports.
inv.write("saintGeorgeStation.xml", format="stationxml", validate=True)