from obspy.core.stream import Stream
from obspy.core import read as obsRead


def ip_read(pathname_or_url=None, format=None, headonly=False, starttime=None,
            endtime=None, nearest_sample=True, dtype=None, apply_calib=False,
            check_compression=True, **kwargs):

    stream = obsRead(pathname_or_url, format=format, headonly=headonly, starttime=starttime,
                     endtime=endtime, nearest_sample=nearest_sample, dtype=dtype, apply_calib=apply_calib,
                     check_compression=check_compression, **kwargs)

    stream.__class__ = IPStream
    stream.promote()

    return stream


class IPStream(Stream):

    __filtered = []

    def __init__(self, traces=None):
        super().__init__(traces)
        self.promote()

    def promote(self):
        self.resetFiteredTraces()

    def getFiltered(self):
        return self.__filtered

    def resetFiteredTraces(self):
        self.__filtered.clear()
        for trace in self.traces:
            self.__filtered.append(trace.copy())

        return
