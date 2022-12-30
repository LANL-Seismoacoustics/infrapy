from scipy import signal
import numpy as np
import pyqtgraph as pg
from obspy import Trace


class IPSpectrogram(pg.ImageItem):
    def __init__():
        super().__init__()


    def setData(self, trace):
        # ingest obspy trace, get the data and calculate the spectrogram
        f, t, Sxx = signal.spectrogram(trace.data, trace.stats['sampling_rate'])

    
