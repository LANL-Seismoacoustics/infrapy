# infrapy.propagation.likelihoods.py
#
# This file contains the likelihood probablility distribution
# functions (pdf's) for the source characteristics given a
# detection.  In general each detection type is expected to
# contain a function det.pdf(...) that returns the likelihood
# distribution.
#
# Author            Philip Blom (pblom@lanl.gov)

from datetime import datetime
import itertools

import json

import numpy as np

import obspy
from obspy import UTCDateTime

from scipy.integrate import simps
from scipy.interpolate import interp1d, interp2d
from scipy.special import i0

from pyproj import Geod

from . import infrasound
from . import seismic
from ..utils import prog_bar


# ################################ #
#    Set Integration Parameters    #
#    and Spherical Earth Model     #
# ################################ #
np.seterr(over='ignore')
int_opts = {'limit': 100, 'epsrel': 1.0e-3}
sph_proj = Geod(ellps='sphere')

# ################################# #
#   InfrasoundDetection class      #
#       Likelihood                  #
# ################################# #


class CustomError(Exception):
    pass


class BadLatitudeError(CustomError):
    """ Rasied when the latitude is out of range """
    pass


class BadLongitudeError(CustomError):
    """ Raised when the longitude is out of range """
    pass


class BadArrayDimensionError(CustomError):
    """ Raised when someone tries to set the array dimension to something that isnt an integer """


class InfrasoundDetection(object):

    __lat = None            # latitude of the center of the array (degrees range from -90 to 90)
    __lon = None            # longitude of the center of the array (degrees range from -180 to 180)
    __ele = None            # average elevation in meters of the array elements
    __peakF_UTCtime = None  # this is the UTC time of the peak F statistic
    __peakF_value = None    # value of the F statistics at the peakF_UTCtime
    __back_az = None        # back azimuth in degrees (range from -180 to 180)
    __start = None          # time in seconds before the peakF_UTCtime that the signal starts
    __end = None            # time in seconds after the peakF_UTCtime that the signal ends
    __freq_range = None     # tuple with the frequency range of the measurement
    __trace_vel = None      # trace velocity in m/s
    __array_dim = None      # integer number of elements in the detecting array
    __method = ''           # method used in the beamforming (eg bartlett)
    __name = ''             # string used to identify the detection (maybe only used in the GUI)
    __event_id = ''         # string containing the event id the detection is associated with
    __note = ''             # optional note that can be added for a reminder etc...
    __network = ''          # network code for the detecting array
    __station = ''          # station code for the detecting array

    # minimum and maximum measured beam widths
    __min_wdth = np.radians(4.0)
    __max_wdth = np.radians(12.0)

    # increase in beam width for propagation effects
    __prop_width = np.radians(4.0)  

    def __init__(self, lat_loc=None, lon_loc=None, time=None, azimuth=None, f_stat=None, array_d=None, f_range=None, start_end=None, note=None, traceV=None, network=None, station=None, method=None):

        self.set_lat(lat_loc)
        self.set_lon(lon_loc)
        self.set_peakF_UTCtime(time)
        self.set_array_dim(array_d)
        self.set_back_azimuth(azimuth)
        self.set_peakF_value(f_stat)
        self.set_freq_range(f_range)
        self.set_startend(start_end)
        self.set_note(note)
        self.set_method(method)
        self.set_trace_velocity(traceV)

        if self.__peakF_value is None or self.__array_dim is None:
            # This happens when creating a InfrasoundDetection object using the json_to_detection method, where an object is created,
            # and then populated later
            return
        else:
            self.calc_kappa_etc()

    def calc_kappa_etc(self):
        # this bit calculates several values used in the methods below. array_dim and peakF_value must have values when called.
        self.kappa0 = 2.0 * (self.__array_dim - 1) * self.__peakF_value
        self.kappa0 = min(self.kappa0, np.log(2.0) / (1.0 - np.cos(self.__min_wdth)))
        self.kappa0 = max(self.kappa0, np.log(2.0) / (1.0 - np.cos(self.__max_wdth)))

        self.cos_half = 1.0 - np.log(2.0) / self.kappa0
        self.sin_half = np.sqrt(1.0 - self.cos_half**2)

        self.kappa = np.log(2.0) / (1.0 - self.cos_half * np.cos(self.__prop_width) +
                                    self.sin_half * np.sin(self.__prop_width))

        self.vm_norm = 2.0 * np.pi * i0(self.kappa)

    def is_equal_to(self, other_detection):
        if self.__lat == other_detection.get_lat() and self.__lon == other_detection.get_lon() and self.__peakF_UTCtime == other_detection.get_peakF_UTCtime():
            # The two detections are pretty much the same, so return True
            return True
        else:
            return False

    def az_pdf(self, lat, lon, path_geo_model=None):
        if len(np.atleast_1d(lat)) == 1:
            temp = sph_proj.inv(self.__lon, self.__lat, lon, lat, radians=False)
            az_diff = np.radians(temp[0] - self.__back_az)

            if path_geo_model:
                rng = temp[2] / 1000.0
                az_diff -= np.radians(path_geo_model.eval_az_dev_mn(rng, self.__back_az - 180.0))
                width_diff = 2.0 * np.radians(path_geo_model.eval_az_dev_vr(rng, self.__back_az - 180.0))
                kappa = np.log(2.0) / (1.0 - self.cos_half * np.cos(width_diff) + self.sin_half * np.sin(width_diff))
                norm = 2.0 * np.pi * i0(kappa)
            else:
                kappa, norm = self.kappa, self.vm_norm

        else:
            temp = np.asarray(sph_proj.inv([self.__lon] * len(lon), [self.__lat] * len(lon), lon, lat, radians=False))
            az_diff = np.radians(temp[0] - [self.__back_az] * len(lon))

            if path_geo_model:
                rng = temp[2] / 1000.0
                az_diff -= np.radians(path_geo_model.eval_az_dev_mn(rng, np.array([self.__back_az] * len(lon)) - 180.0))
                width_diff = 2.0 * np.radians(path_geo_model.eval_az_dev_vr(rng, np.array([self.__back_az] * len(lon)) - 180.0))
                kappa = np.log(2.0) / (1.0 - self.cos_half * np.cos(width_diff) + self.sin_half * np.sin(width_diff))
                norm = 2.0 * np.pi * i0(kappa)
            else:
                kappa, norm = self.kappa, self.vm_norm

        return np.exp(kappa * np.cos(az_diff)) / norm

    def rng_pdf(self, lat, lon, t, path_geo_model=None):
        if len(np.atleast_1d(lat)) == 1:
            rng = sph_proj.inv(self.__lon, self.__lat, lon, lat, radians=False)[2] / 1000.0
            if path_geo_model:
                result = path_geo_model.eval_rcel_gmm(rng, (self.__peakF_UTCtime - t).astype('m8[s]').astype(float) / rng, self.__back_az - 180.0)
            else:
                result = infrasound.canonical_rcel((self.__peakF_UTCtime - t).astype('m8[s]').astype(float) / rng)
        else:
            rng = np.asarray(sph_proj.inv([self.__lon] * len(lon), [self.__lat] * len(lon), lon, lat, radians=False))[2] / 1000.0
            if path_geo_model:
                result = path_geo_model.eval_rcel_gmm(rng, (self.__peakF_UTCtime - t).astype('m8[s]').astype(float) / rng, np.asarray([self.__back_az - 180.0] * len(lon)))
            else:
                result = infrasound.canonical_rcel((self.__peakF_UTCtime - t).astype('m8[s]').astype(float) / rng)
        return result / rng

    def pdf(self, lat, lon, t, path_geo_model=None):
        if self.back_azimuth is None:
            return self.rng_pdf(lat, lon, t, path_geo_model)
        elif self.__peakF_UTCtime > UTCDateTime("9999-01-01T00:00:00"):
            return self.az_pdf(lat, lon, path_geo_model)
        else:
            if len(np.atleast_1d(lat)) == 1:
                temp = sph_proj.inv(self.__lon, self.__lat, lon, lat, radians=False)
                az_diff = np.radians(temp[0] - self.__back_az)
                rng = temp[2] / 1000.0

                if path_geo_model:
                    az_diff -= np.radians(path_geo_model.eval_az_dev_mn(rng, self.__back_az - 180.0))
                    width_diff = 2.0 * np.radians(path_geo_model.eval_az_dev_vr(rng, self.__back_az - 180.0))
                    kappa = np.log(2.0) / (1.0 - self.cos_half * np.cos(width_diff) + self.sin_half * np.sin(width_diff))
                    norm = 2.0 * np.pi * i0(kappa)

                    az_pdf = np.exp(kappa * np.cos(az_diff)) / norm
                    rng_pdf = path_geo_model.eval_rcel_gmm(rng, (self.__peakF_UTCtime - t).astype('m8[s]').astype(float) / rng, self.__back_az - 180.0)
                else:
                    az_pdf = np.exp(self.kappa * np.cos(az_diff)) / self.vm_norm
                    rng_pdf = infrasound.canonical_rcel((self.get_peakF_UTCtime() - t).astype('m8[s]').astype(float) / rng)
            else:
                temp = np.asarray(sph_proj.inv([self.__lon] * len(lon), [self.__lat] * len(lon), lon, lat, radians=False))
                az_diff = np.radians(temp[0] - [self.__back_az] * len(lon))
                rng = temp[2] / 1000.0

                if path_geo_model:
                    az_diff -= np.radians(path_geo_model.eval_az_dev_mn(rng, np.array([self.__back_az] * len(lon)) - 180.0))
                    width_diff = 2.0 * np.radians(path_geo_model.eval_az_dev_vr(rng, np.array([self.__back_az] * len(lon)) - 180.0))
                    kappa = np.log(2.0) / (1.0 - self.cos_half * np.cos(width_diff) + self.sin_half * np.sin(width_diff))
                    norm = 2.0 * np.pi * i0(kappa)

                    az_pdf = np.exp(kappa * np.cos(az_diff)) / norm
                    rng_pdf = path_geo_model.eval_rcel_gmm(rng, (self.__peakF_UTCtime - t).astype('m8[s]').astype(float) / rng, np.asarray([self.__back_az - 180.0] * len(lon)))
                else:
                    az_pdf = np.exp(self.kappa * np.cos(az_diff)) / self.vm_norm
                    rng_pdf = infrasound.canonical_rcel((self.__peakF_UTCtime - t).astype('m8[s]').astype(float) / rng)

            return az_pdf * rng_pdf / rng

    def src_spec_pdf(self, lat, lon, freqs, src_spec, smn_spec, tloss_models):
        """Defines the probability of detection being produced by a given source

            Computes the probability of a given source spectral amplitude, P(f),
            producing this detection given its signal-minus-noise spectrum and
            a transmission loss model

            Parameters
            ----------
            lat : float
                Latitude of the hypothetical source
            lon : float
                Longitue of the hypothetical source
            freqs : 1darray
                Frequencies at which to evaluate the source spectrum
            src_spec : 1darray
                Values at which to evaluate the source spectrum
            smn_spec : 2darray
                Numpy array containing frequencies and signal-minus-noise spectral amplitudes of the detection
            tloss_models : frequencies and TLossModel instances
                Transmission loss models to use and frequencies at which they are computed

            Returns:
            ----------
            freq_grid : 1darray
                Frequencies at which the pdf is evaluated
            src_spec_grid : 1darray
                Spectral amplitudes at which the pdf is evaluated
            pdf : 2darray
                Grid of probability for source spectrum at frequencies and source spectrum amplitudes specified by input
            """

        # define the source-receiver range in km and propagation azimuth from
        # detection back azimuths - 180.0 degrees
        rng = sph_proj.inv(self.longitude, self.latitude, lon, lat, radians=False)[2] / 1000.0

        prop_az = self.back_azimuth - 180.0
        if prop_az < -180.0:
            prop_az += 360.0
        elif prop_az > 180.0:
            prop_az -= 360.0

        # Interpolate the transmission loss at the source-receiver range and
        # azimuth to define tloss_interp(f, tloss)
        tloss_freqs = np.asarray(tloss_models[0])
        tloss_vals =  np.linspace(-100.0, 0.0, 100)
        tloss_grid = np.empty((len(tloss_freqs), len(tloss_vals)))

        for nf in range(len(tloss_freqs)):
            tloss_grid[nf] = tloss_models[1][nf].eval(np.array([rng] * len(tloss_vals)), tloss_vals, np.array([prop_az] * len(tloss_vals)))

        tloss_interp = interp2d(tloss_freqs, tloss_vals, tloss_grid.T, bounds_error=False, fill_value=0.0)

        smn_interp = interp1d(smn_spec[0], smn_spec[1])

        freq_grid, spec_grid = np.meshgrid(freqs, src_spec)
        F = freq_grid.flatten()
        SPEC = spec_grid.flatten()

        pdf = np.array([tloss_interp(F[n], smn_interp(F[n]) - SPEC[n]) for n in range(len(F))]).T[0]

        return F, SPEC, pdf

    # ################################
    # Accessor methods.  Useful to check the type and values before setting the local variables
    #
    # It's recommended using these instead of setting the values directly
    #
    # ################################

    def get_lat(self):
        # returns the latitude of the center of the detecting array
        return self.__lat

    def set_lat(self, lat):
        if lat is None:
            self.__lat = None
            return

        # lat is the latitude in degrees from -90 to 90
        if lat < -90 or lat > 90:
            raise BadLatitudeError
        else:
            self.__lat = lat

    latitude = property(get_lat, set_lat, doc="Latitude of the array")

    def get_lon(self):
        # returns the longitude of the center of the detecting array
        return self.__lon

    def set_lon(self, lon):
        if lon is None:
            self.__lon = None
            return

        # lon is the longitude in degrees from -180 to 180
        if lon < -180 or lon > 180:
            raise BadLongitudeError
        else:
            self.__lon = lon

    longitude = property(get_lon, set_lon, doc="Longitude of the detecting array")

    def get_name(self):
        return self.__name

    def set_name(self, n):
        self.__name = n

    name = property(get_name, set_name, doc="(optional) Name of the detection")

    def set_ele(self, ele):
        self.__ele = ele

    def get_ele(self):
        # returns the elevation
        return self.__ele

    elevation = property(get_ele, set_ele, doc="Elevation of the detecting array")

    def get_trace_velocity(self):
        return self.__trace_vel

    def set_trace_velocity(self, tv):
        self.__trace_vel = tv

    trace_velocity = property(get_trace_velocity, set_trace_velocity, doc="Trace velocity of the detection")

    def set_start(self, start):
        self.__start = start

    def get_start(self):
        # returns the start of the signal in seconds relative to the peak F utc time
        return self.__start

    start = property(get_start, set_start, doc="(float) Start time in seconds relative to the peak F utc time of the detection")

    def set_end(self, end):
        self.__end = end

    def get_end(self):
        # returns the end of the signal in seconds relative to the peak F utc tme
        return self.__end

    end = property(get_end, set_end, doc="(float) End time in seconds relative to the peak F utc time of the detection")

    def set_startend(self, start_end):
        #start_end is a tuple containing the start and end values of the detection
        if start_end is None:
            self.start = None
            self.end = None
        else:
            self.start = start_end[0]
            self.end = start_end[1]

    def get_startend(self):
        return (self.start, self.end)

    start_end = property(get_startend, set_startend, doc="(tuple) starting and ending values for the detection")

    def set_freq_range(self, frange):
        self.__freq_range = frange

    def get_freq_range(self):
        return self.__freq_range

    frequency_range = property(get_freq_range, set_freq_range, doc="(tuple) Frequency range of the measurement")

    def get_array_dim(self):
        return self.__array_dim

    def set_array_dim(self, ad):
        self.__array_dim = ad

    array_dim = property(get_array_dim, set_array_dim, doc="(int) Number of elements in the detecting array")

    def set_back_azimuth(self, ba):
        self.__back_az = ba

    def get_back_azimuth(self):
        return self.__back_az

    back_azimuth = property(get_back_azimuth, set_back_azimuth, doc="(float) Back azimuth of the detected signal")

    def set_note(self, note):
        self.__note = note

    def get_note(self):
        return self.__note

    note = property(get_note, set_note, doc="(optional) Note regarding the detection")

    def set_event_id(self, evt):
        self.__event_id = evt

    def get_event_id(self):
        return self.__event_id

    event_id = property(get_event_id, set_event_id, doc="(optional) ID of the event")

    def set_method(self, method):
        self.__method = method

    def get_method(self):
        return self.__method

    method = property(get_method, set_method, doc="(string) Beamforming method used to generate the detection")

    def get_peakF_UTCtime(self, type='numpy'):
        # this returns the peak F UTC time with the chosen type
        # types can be 'numpy', 'datetime', or 'obspy'
        if type == 'numpy':
            return self.__peakF_UTCtime
        elif type == 'datetime':
            return self.__peakF_UTCtime.astype(datetime)
        elif type == 'obspy':
            # first generate a datetime object
            dt = self.__peakF_UTCtime.astype(datetime)
            # then create a obspy UTCDateTime from that and return it
            return UTCDateTime(dt)

    def set_peakF_UTCtime(self, pf_UTCtime):
        if pf_UTCtime is None:
            self.__peakF_UTCtime = None
            return

        # pf_UTCtime is the peak F time in UTC.  This can be either a datetime.datetime object, a obspy.UTCDateTime object, or a numpy.datetime object.
        # This method will detect which type it is and convert it to a numpy.datetime object for internal use
        if type(pf_UTCtime) is np.datetime64:
            self.__peakF_UTCtime = pf_UTCtime
        elif type(pf_UTCtime) is datetime:
            self.__peakF_UTCtime = np.datetime64(pf_UTCtime)
        elif type(pf_UTCtime) is obspy.core.utcdatetime.UTCDateTime:
            # converting froma obspy datetime object to a numpy datetime object...
            # first convert from obspy to datetime...
            dt = pf_UTCtime._get_datetime()
            # then create a numpy datetime from that...
            self.__peakF_UTCtime = np.datetime64(dt)

    peakF_UTCtime = property(get_peakF_UTCtime, set_peakF_UTCtime, doc="time of the peak F value, in UTC")

    def set_peakF_value(self, pfv):
        self.__peakF_value = pfv

    def get_peakF_value(self):
        # returns the peak F value
        return self.__peakF_value

    peakF_value = property(get_peakF_value, set_peakF_value, doc="Peak F values")

    def set_network(self, network):
        self.__network = network

    def get_network(self):
        return self.__network

    network = property(get_network, set_network, doc="(optional) Network code")

    def set_station(self, station):
        self.__station = station

    def get_station(self):
        return self.__station

    station = property(get_station, set_station, doc="(optional) Station code")


    # #### End Accessor methods #################

    def generateDict(self):
        detectionDict = {'Name': self.__name,
                         'Time (UTC)': self.__peakF_UTCtime.astype(datetime).isoformat(),
                         'F Stat.': self.__peakF_value,
                         'Trace Vel. (m/s)': self.__trace_vel,
                         'Back Azimuth': self.__back_az,
                         'Latitude': self.__lat,
                         'Longitude': self.__lon,
                         'Elevation (m)': self.__ele,
                         'Start': self.__start,
                         'End': self.__end,
                         'Freq Range': self.__freq_range,
                         'Array Dim.': self.__array_dim,
                         'Method': self.__method,
                         'Event': self.__event_id,
                         'Note': self.__note,
                         'Network': self.__network,
                         'Station': self.__station}
        return detectionDict

    def fillFromDict(self, dict):
        if 'Latitude' in dict:
            self.set_lat(dict['Latitude'])
        if 'Longitude' in dict:
            self.set_lon(dict['Longitude'])
        if 'Elevation' in dict:  # This is here for legacy data
            self.set_ele(dict['Elevation'])
        if 'Elevation (m)' in dict:
            self.set_ele(dict['Elevation (m)'])
        if 'Name' in dict:
            self.set_name(dict['Name'])
        if 'Peak F-Stat Time (UTC)' in dict:    # This is here for legacy data
            if dict['Peak F-Stat Time (UTC)'] is not None:
                self.set_peakF_UTCtime(UTCDateTime(dict['Peak F-Stat Time (UTC)']))
            else:
                self.set_peakF_UTCtime(UTCDateTime("9999-01-02T00:00:00"))
        if 'Time (UTC)' in dict:
            if dict['Time (UTC)'] is not None:
                self.set_peakF_UTCtime(UTCDateTime(dict['Time (UTC)']))
            else:
                self.set_peakF_UTCtime(UTCDateTime("9999-01-02T00:00:00"))
        if 'F Statistic' in dict:   # This is here for legacy data
            self.set_peakF_value(dict['F Statistic'])
        if 'F Stat.' in dict:
            self.set_peakF_value(dict['F Stat.'])
        if 'Trace Velocity (m/s)' in dict:  # This is here for legacy data
            self.set_trace_velocity(dict['Trace Velocity (m/s)'])
        if 'Trace Vel. (m/s)' in dict:
            self.set_trace_velocity(dict['Trace Vel. (m/s)'])
        if 'Back Azimuth' in dict:
            self.set_back_azimuth(dict['Back Azimuth'])
        if 'Note' in dict:
            self.set_note(dict['Note'])
        if 'Event' in dict:
            self.set_event_id(dict['Event'])
        if 'Start' in dict:
            self.set_start(dict['Start'])
        if 'End' in dict:
            self.set_end(dict['End'])
        if 'Freq Range' in dict:
            self.set_freq_range(dict['Freq Range'])
        if 'Array Dimension' in dict:   # This is here for legacy data
            self.set_array_dim(dict['Array Dimension'])
        if 'Array Dim.' in dict:
            self.set_array_dim(dict['Array Dim.'])
        if 'Method' in dict:
            self.set_method(dict['Method'])
        if 'Network' in dict:
            self.set_network(dict['Network'])
        if 'Station' in dict:
            self.set_station(dict['Station'])

        # At this point, the array_dim, and peakF_value have probably changed, so we need to recalculate
        # some of the variables...
        self.calc_kappa_etc()

    def __str__(self):
        return '{}'.format(self.generateDict())


# ############################# #
#   SeismicDetection class      #
#       Likelihood              #
# ############################# #
class SeismicDetection(object):
    def __init__(self, lat, lon, time, phase="p", slowness=None):
        self.__lat = lat
        self.__lon = lon
        self.__peakF_UTCtime = time
        self.__phase = phase
        self.__slowness = slowness

    def trvl_tm(self, rng):
        if self.__phase == "p" or self.__phase == "P":
            return seismic.ak135_p_tr_time(rng)
        elif self.__phase == "s" or self.__phase == "S":
            return seismic.ak135_s_tr_time(rng)
        else:
            msg = "Unrecognized phase identifier in seismic detection. Current options are 'p' and 's'."
            raise ValueError(msg)

    def sigma(self):
        if self.__phase == "p" or self.__phase == "P":
            return seismic.ak135_p_sigma
        elif self.__phase == "s" or self.__phase == "S":
            return seismic.ak135_s_sigma
        else:
            msg = "Unrecognized phase identifier in seismic detection. Current options are 'p' and 's'."
            raise ValueError(msg)

    def pdf(self, lat, lon, t):
        if len(np.atleast_1d(lat)) == 1:
            rng = sph_proj.inv(self.__lon, self.__lat, lon, lat)[2] * 1.0e-3
        else:
            rng = sph_proj.inv([self.__lon] * len(lon), [self.__lat] * len(lat), lon, lat)[2] * 1.0e-3

        dt = (self.__peakF_UTCtime - t).astype('m8[ms]').astype(float) * 1.0e-3
        return np.exp(-1.0 / 2.0 * ((dt - self.trvl_tm(rng)) / self.sigma())**2) / np.sqrt(2.0 * np.pi * self.sigma()**2)


# ################################ #
#        Methods to Combine        #
#      Infrasound Likelihoods      #
# ################################ #
def joint_pdf(lat, lon, t, det_list, path_geo_model=None, prog_step=0):
    result = np.array([det.pdf(lat, lon, t, path_geo_model=path_geo_model) for det in det_list if type(det) == InfrasoundDetection]).prod(axis=0)
    result *= np.array([det.pdf(lat, lon, t) for det in det_list if type(det) == SeismicDetection]).prod(axis=0)
    prog_bar.increment(n=prog_step)
    return result

def joint_pdf_wrapper(args):
    return joint_pdf(*args)


def marginal_spatial_pdf(lat, lon, det_list, path_geo_model=None, prog_step=0, resol=100):
    # define seprate lists of infrasound and seismic only detections
    infr_det_list = [det for det in det_list if type(det) == InfrasoundDetection]
    infr_det_list2 = [det for det in infr_det_list if det.peakF_UTCtime < UTCDateTime("9999-01-01T00:00:00")]

    seis_det_list = [det for det in det_list if type(det) == SeismicDetection]

    infr2_cnt = len(infr_det_list2)
    seis_cnt = len(seis_det_list)

    if 3**infr2_cnt > resol:
        def temp(la, lo):
            infr_rngs = np.array([sph_proj.inv(det.longitude, det.latitude, lo, la)[2] / 1000.0 for det in infr_det_list2])
            infr_t1 = max(np.array([det.peakF_UTCtime - np.timedelta64(int(infr_rngs[n] / 0.2 * 1e3), 'ms') for n, det in enumerate(infr_det_list2)]))
            infr_t2 = min(np.array([det.peakF_UTCtime - np.timedelta64(int(infr_rngs[n] / 0.4 * 1e3), 'ms') for n, det in enumerate(infr_det_list2)]))

            if seis_cnt > 0:
                seis_rngs = np.array([sph_proj.inv(det.longitude, det.latitude, lo, la)[2] / 1000.0 for det in seis_det_list])
                seis_t1 = max(np.array([det.peakF_UTCtime - np.timedelta64(int(seis_rngs[n] / 2.0 * 1e3), 'ms') for n, det in enumerate(seis_det_list)]))
                seis_t2 = min(np.array([det.peakF_UTCtime - np.timedelta64(int(seis_rngs[n] / 8.0 * 1e3), 'ms') for n, det in enumerate(seis_det_list)]))
            else:
                seis_t1 = np.datetime64("0000-01-01T00:00:00")
                seis_t2 = np.datetime64("9999-01-01T00:00:00")

            t1 = max(infr_t1, seis_t1)
            t2 = min(infr_t2, seis_t2)

            t_vals = np.array([t1 + (t2 - t1) / (resol - 1) * m for m in range(resol)])
            pdf_vals = joint_pdf(np.array([la] * resol),np.array([lo] * resol), t_vals, det_list, path_geo_model=path_geo_model)

            return simps(pdf_vals, (t2 - t1).astype('m8[ms]').astype(float) * 1.0e-3 / (resol - 1) * range(resol))

        if len(np.atleast_1d(lat)) == 1:
            return temp(lat, lon)
        else:
            return np.array([temp(la, lon[n]) for n, la in enumerate(lat)])

    else:
        # compute azimuthal distribution for all infrasound detections        
        az_pdf = np.array([det.az_pdf(lat, lon, path_geo_model) if det.back_azimuth is not None else np.ones_like(lat) for det in infr_det_list]).prod(axis=0)

        if infr2_cnt > 0:
            # pull out the latitudes and longitudes for infrasonic and seismic detections
            infr_lats = np.array([det.latitude for det in infr_det_list2])
            infr_lons = np.array([det.longitude for det in infr_det_list2])
            infr_tms = np.array([(det.peakF_UTCtime - det_list[0].peakF_UTCtime).astype('m8[ms]').astype(float) * 1.0e-3 for det in infr_det_list2])

            seis_lats = np.array([det.latitude for det in seis_det_list])
            seis_lons = np.array([det.longitude for det in seis_det_list])
            seis_tms = np.array([(det.peakF_UTCtime - det_list[0].peakF_UTCtime).astype('m8[ms]').astype(float) * 1.0e-3 for det in seis_det_list])

            # Compute the index sequences for infrasound likelihoods
            sequences = []
            for seq in itertools.product(list(range(3)), repeat=infr2_cnt):
                sequences = sequences + [list(seq)]
            sequences = np.array(sequences, dtype=int)

            if len(np.atleast_1d(lat)) == 1:
                infr_rngs = sph_proj.inv(infr_lons, infr_lats, np.array([lon] * infr2_cnt), np.array([lat] * infr2_cnt), radians=False)[2] / 1000.0
                seis_rngs = sph_proj.inv(seis_lons, seis_lats, np.array([lon] * seis_cnt), np.array([lat] * seis_cnt), radians=False)[2] / 1000.0
                if path_geo_model:
                    mns = np.empty((infr2_cnt, 3))
                    vrs = np.empty((infr2_cnt, 3))
                    wts = np.empty((infr2_cnt, 3))

                    for n, (rng_n, az_n) in enumerate(zip(infr_rngs, np.array([det.__back_az for det in infr_det_list]))):
                        az_bin = infrasound.find_azimuth_bin(az_n - 180.0)

                        mns[n] = np.array([path_geo_model.rcel_mns[az_bin][m](min(rng_n, path_geo_model.rng_max)) for m in range(3)])
                        vrs[n] = np.array([path_geo_model.rcel_vrs[az_bin][m](min(rng_n, path_geo_model.rng_max)) for m in range(3)])
                        wts[n] = np.array([path_geo_model.rcel_wts[az_bin][m](min(rng_n, path_geo_model.rng_max)) for m in range(3)])
                else:
                    mns = np.array([infrasound.canon_rcel_mns] * infr2_cnt)
                    vrs = np.array([infrasound.canon_rcel_vrs] * infr2_cnt)
                    wts = np.array([infrasound.canon_rcel_wts] * infr2_cnt)
                temp1 = np.array([1.0 / det.sigma()**2 for det in seis_det_list])
                temp2 = np.array([seis_tms[n] - det.trvl_tm(seis_rngs[n]) for n, det in enumerate(seis_det_list)])

                a_seis = temp1.sum()
                b_seis = (temp2 * temp1).sum()
                c_seis = (temp2**2 * temp1).sum()

                rng_cel_pdf = 0.0
                for seq in sequences:
                    a, b, c, N = a_seis, b_seis, c_seis, 1.0

                    for n, ni in enumerate(zip(seq)):
                        dt = infr_tms[n] - infr_rngs[n] * mns[n][ni]
                        sig = infr_rngs[n] * vrs[n][ni]

                        a += 1.0 / sig**2
                        b += dt / sig**2
                        c += (dt / sig)**2
                        N *= wts[n][ni] / sig

                    rng_cel_pdf += N / np.sqrt(a) * np.exp(-1.0 / 2.0 * (c - b**2 / a))

                rng_cel_pdf /= np.power(2.0 * np.pi, (seis_cnt + infr2_cnt - 1.0) / 2.0)
                prog_bar.increment(n=prog_step)
            else:
                infr_rngs = sph_proj.inv(np.array([infr_lons] * len(lon)), np.array([infr_lats] * len(lat)), np.array([lon] * infr2_cnt).T, np.array([lat] * infr2_cnt).T, radians=False)[2].reshape(len(lat), infr2_cnt) / 1000.0
                seis_rngs = sph_proj.inv(np.array([seis_lons] * len(lon)), np.array([seis_lats] * len(lat)), np.array([lon] * seis_cnt).T, np.array([lat] * seis_cnt).T, radians=False)[2].reshape(len(lat), seis_cnt) / 1000.0

                if path_geo_model:
                    mns = np.empty((len(lat), infr2_cnt, 3))
                    vrs = np.empty((len(lat), infr2_cnt, 3))
                    wts = np.empty((len(lat), infr2_cnt, 3))

                    infr_azs = np.array([np.array([det.back_azimuth for det in det_list if type(det) == InfrasoundDetection])] * len(lat))
                    az_bins = infrasound.find_azimuth_bin(infr_azs.flatten() - 180.0, path_geo_model.az_bin_cnt)
                    rngs_eval = infr_rngs.flatten()
                    rngs_eval[rngs_eval > path_geo_model.rng_max] = path_geo_model.rng_max

                    for k in range(3):
                        mns_temp = np.empty_like(az_bins, dtype=float)
                        vrs_temp = np.empty_like(az_bins, dtype=float)
                        wts_temp = np.empty_like(az_bins, dtype=float)

                        for n_az in range(path_geo_model.az_bin_cnt):
                            if np.any(az_bins == n_az):
                                mns_temp[az_bins == n_az] = path_geo_model.rcel_mns[n_az][k](rngs_eval[az_bins==n_az])
                                vrs_temp[az_bins == n_az] = path_geo_model.rcel_vrs[n_az][k](rngs_eval[az_bins==n_az])
                                wts_temp[az_bins == n_az] = path_geo_model.rcel_wts[n_az][k](rngs_eval[az_bins==n_az])

                        mns[:, :, k] = mns_temp.reshape(len(lat), infr2_cnt)
                        vrs[:, :, k] = vrs_temp.reshape(len(lat), infr2_cnt)
                        wts[:, :, k] = wts_temp.reshape(len(lat), infr2_cnt)

                else:
                    mns = np.array([np.array([infrasound.canon_rcel_mns] * infr2_cnt)] * len(lat))
                    vrs = np.array([np.array([infrasound.canon_rcel_vrs] * infr2_cnt)] * len(lat))
                    wts = np.array([np.array([infrasound.canon_rcel_wts] * infr2_cnt)] * len(lat))

                temp1 = np.array([[1.0 / det.sigma()**2] * len(lat) for det in seis_det_list])
                temp2 = np.array([seis_tms[n] - det.trvl_tm(seis_rngs[:, n]) for n, det in enumerate(seis_det_list)])

                a_seis = temp1.sum(axis=0)
                b_seis = (temp2 * temp1).sum(axis=0)
                c_seis = (temp2**2 * temp1).sum(axis=0)

                rng_cel_pdf = 0.0
                for seq in sequences:
                    a, b, c, N = a_seis, b_seis, c_seis, 1.0

                    for n in range(infr2_cnt):
                        dt = infr_tms[n] - infr_rngs[:, n] * mns[:, n, seq[n]]
                        sig = infr_rngs[:, n] * vrs[:, n, seq[n]]

                        a = a + 1.0 / sig**2
                        b = b + dt / sig**2
                        c = c + (dt / sig)**2
                        N = N * wts[:, n, seq[n]] / sig

                    rng_cel_pdf += N / np.sqrt(a) * np.exp(-1.0 / 2.0 * (c - b**2 / a))
                rng_cel_pdf /= np.power(2.0 * np.pi, (seis_cnt + infr2_cnt - 1.0) / 2.0)
                prog_bar.increment(n=prog_step)
        else:
            rng_cel_pdf = 1.0

        return az_pdf * rng_cel_pdf

def marginal_spatial_pdf_wrapper(args):
    return marginal_spatial_pdf(*args)
