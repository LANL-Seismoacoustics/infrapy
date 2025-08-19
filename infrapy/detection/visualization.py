# visualization.py
#
# Visualization methods for beamforming (fk) and detection (fd) results
#
# Philip Blom (pblom@lanl.gov)


import numpy as np

from obspy import UTCDateTime

from scipy.interpolate import interp1d

from scipy.signal import spectrogram, stft, cwt, morlet2

import matplotlib.pyplot as plt 
from matplotlib import cm

from . import beamforming_new



def plot_fk_json(det_dict, output_path=None, show_fig=True):

    fig, a = plt.subplots(2, figsize=(10, 4), sharex=True)
    a[1].set_ylabel("Tr. Vel. [m/s]")
    a[0].set_ylabel("Back Az. [deg.]")

    f_mean = np.mean(det_dict["fk"]["f-stat"])
    f_std = np.std(det_dict["fk"]["f-stat"])
    cmax = f_mean + 2.0 * f_std
    
    for k, det in enumerate(det_dict['det_info']):
        t1 = np.datetime64(det["peak f-stat time"]) + np.timedelta64(int(det["start/end"][0][0] * 1000.0), 'ms')
        t2 = np.datetime64(det["peak f-stat time"]) + np.timedelta64(int(det["start/end"][0][1] * 1000.0), 'ms')
        for ax in a:
            ax.axvspan(t1, t2, color="steelblue", alpha=0.25)
        
        a[0].annotate(str(k), (t1, 180.0), xytext=(2, 0), textcoords='offset pixels', horizontalalignment = "left", verticalalignment="top")
            

    t_vals = [np.datetime64(det_dict["wvfrm_info"][0][0]["starttime"]) + np.timedelta64(int(dt * 1000.0), 'ms') for dt in det_dict["fk"]["time"]]

    a[1].scatter(t_vals, det_dict["fk"]["tr vel"], c=det_dict["fk"]["f-stat"], cmap=cm.jet, s=5.0, vmax=cmax)
    sc_plt = a[0].scatter(t_vals, det_dict["fk"]["back az"], c=det_dict["fk"]["f-stat"], cmap=cm.jet, s=5.0, vmax=cmax)

    fig.colorbar(sc_plt, ax=a.ravel().tolist(), label="f-stat", location="right", aspect=20)

    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()


def plot_det_json(det_dict, single_det_index, param_set_index=None, output_path=None, show_fig=True):

    det_info = det_dict["det_info"][single_det_index]
        
    if param_set_index is None:
        param_set_index = np.argmax(np.array([np.max(fk_j["f-stat"]) for fk_j in det_info["fk"]]))

    for fk_j in det_info["fk"]:
        dt_vals = [np.datetime64(det_info["peak f-stat time"]) + np.timedelta64(int(dt * 1000.0), 'ms') for dt in fk_j["time"]]

    dt_min = min([(fk_j["time"][0]) for fk_j in det_info["fk"]])
    dt_max = max([(fk_j["time"][-1]) for fk_j in det_info["fk"]])

    t1 = np.datetime64(det_info["peak f-stat time"]) + np.timedelta64(int(dt_min * 1000.0), 'ms')
    t2 = np.datetime64(det_info["peak f-stat time"]) + np.timedelta64(int(dt_max * 1000.0), 'ms')

    t_lims = [t1, t2]

    fig = plt.figure(figsize=(12, 6), layout="constrained")
    spec = fig.add_gridspec(4, 5)

    ax1 = fig.add_subplot(spec[0, :3])
    ax2 = fig.add_subplot(spec[1, :3], sharex=ax1)
    ax3 = fig.add_subplot(spec[2, :3], sharex=ax1)
    ax4 = fig.add_subplot(spec[3, :3], sharex=ax1)

    ax1.set_xlim(t_lims)

    ax2.set_xlabel("")
    ax3.set_xlabel("")
    ax4.set_xlabel("")
    ax4.tick_params(axis='x', labelrotation=30)

    ax1.set_ylabel("Pressure [Pa]")
    ax2.set_ylabel("Back Azimuth [deg]")
    ax3.set_ylabel("Tr. Velocity [m/s]")
    ax4.set_ylabel("f-stat")


    t1 = np.datetime64(det_info["peak f-stat time"]) + np.timedelta64(int(det_info["start/end"][param_set_index][0] * 1000.0), 'ms')
    t2 = np.datetime64(det_info["peak f-stat time"]) + np.timedelta64(int(det_info["start/end"][param_set_index][1] * 1000.0), 'ms')

    for ax in [ax1, ax2, ax3, ax4]:
        ax.axvspan(t1, t2, color="lightsteelblue", alpha=0.5)

    t_vals = [np.datetime64(det_info["peak f-stat time"]) + np.timedelta64(int(dt * 1000.0), 'ms') for dt in det_info["beam"][param_set_index]["time"]]
    ax1.fill_between(t_vals, -np.array(det_info["beam"][param_set_index]['resid']), det_info["beam"][param_set_index]['resid'], color='r', alpha=0.5)
    ax1.plot(t_vals, det_info["beam"][param_set_index]['signal'], 'k', linewidth=0.5)

    for fk_j in det_info["fk"]:
        dt_vals = [np.datetime64(det_info["peak f-stat time"]) + np.timedelta64(int(dt * 1000.0), 'ms') for dt in fk_j["time"]]

        ax2.plot(dt_vals, fk_j["back az"], 'o', markersize=3)
        ax3.plot(dt_vals, fk_j["tr vel"], 'o', markersize=3)
        ax4.plot(dt_vals, fk_j["f-stat"], 'o', markersize=3)

    for ax in [ax1, ax2, ax3]:
        plt.setp(ax.get_xticklabels(), visible=False)
        
    ax5 = fig.add_subplot(spec[2:, 3:])
    
    if "fk_params" in det_dict:
        f_min = det_dict["fk_params"][param_set_index]["freq_min"]
        f_max = det_dict["fk_params"][param_set_index]["freq_max"]
        trace_id = det_dict["wvfrm_info"][param_set_index][0]["trace id"]
        trace_cnt = len(det_dict["wvfrm_info"][param_set_index])
    else:
        f_min = min([fk_j["freq_min"] for fk_j in det_info["fk_params"]])
        f_max = max([fk_j["freq_max"] for fk_j in det_info["fk_params"]])
        trace_id = det_info["wvfrm_info"][param_set_index][0]["trace id"]
        trace_cnt = len(det_info["wvfrm_info"][param_set_index])

    ax5.set_xlim(max(f_min / 5, det_info["spec"][param_set_index]["freq"][10]), min(f_max * 5, det_info["spec"][param_set_index]["freq"][-1]))

    ax5.set_xlabel("Frequency [Hz]")
    ax5.set_ylabel("Spectral Amplitude [Pa/Hz]")

    ax5.yaxis.set_label_position("right")
    ax5.yaxis.set_ticks_position("right")

    # combine spectral and residual curves across detection parameter sets
    if len(det_info["spec"]) > 1:
        f_lim = np.min([det_spec["freq"][-1] for det_spec in det_info["spec"]])       
        df = np.min([det_spec["freq"][1] for det_spec in det_info["spec"]])

        signal_interps = [interp1d(det_spec["freq"], det_spec["signal"]) for det_spec in det_info["spec"]]
        resid_interps = [interp1d(det_spec["freq"], det_spec["resid"]) for det_spec in det_info["spec"]]

        spec_freq = np.arange(0.0, f_lim, df)
        spec_sig = np.mean(np.array([spec_fit(spec_freq) for spec_fit in signal_interps]), axis=0)
        spec_resid = np.mean(np.array([resid_fit(spec_freq) for resid_fit in resid_interps]), axis=0)
    else:
        spec_freq = det_info["spec"][param_set_index]["freq"]
        spec_sig = det_info["spec"][param_set_index]["signal"]
        spec_resid = det_info["spec"][param_set_index]["resid"]

    ax5.loglog(spec_freq, spec_resid, '-r', linewidth=0.5, label="Residual")
    ax5.loglog(spec_freq, spec_sig, '-k', linewidth=0.5, label="Signal Estimate")
    ax5.axvspan(f_min, f_max, color='lightsteelblue', alpha=0.5, edgecolor=None, zorder=1)
    ax5.legend(fontsize=9)

    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.grid(True, zorder=10)

    ax6 = fig.add_subplot(spec[:2, 3:])
    ax6.axis('off')

    det_summary = 'Detection Summary\n' + '-' * 30
    det_summary = det_summary + '\nTrace ID: ' + trace_id + " (1st of " + str(trace_cnt) + " traces)"
    det_summary = det_summary + '\nDate/Time: ' + det_info["peak f-stat time"]
    det_summary = det_summary + '\nf-stat: ' + str(np.round(det_info["f-stat"],1))
    det_summary = det_summary + '\nBack Azimuth: ' + str(np.round(det_info["back az"],1)) + " deg (rel. N)"
    det_summary = det_summary + '\nTr. Velocity: ' + str(np.round(det_info["tr vel"],1)) + " m/s"
    det_summary = det_summary + '\nDuration: ' + str(det_info["start/end"][param_set_index][1] - det_info["start/end"][param_set_index][0]) + " sec"
    det_summary = det_summary + '\nFrequency Band: ' + str(f_min) + " - " + str(f_max) + " Hz"
   
    ax6.text(0.0, 1.0, det_summary, va="top", fontsize=12, bbox=dict(boxstyle="round, pad=0.25", fc="lightsteelblue", ec="black", lw=1))
    
    for ax in [ax1, ax2, ax3, ax4]:
        ax.axvspan(t1, t2, color="lightsteelblue", alpha=0.75)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=250) 

    if show_fig:
        plt.show()


def plot_wvfrms(det_dict, annotate_option="frequency", output_path=None, show_fig=True):

    param_set_cnt = len(det_dict["fk"])
        
    t0 = np.datetime64(det_dict["peak f-stat time"])

    dt_min = min([(beam_j["time"][0]) for beam_j in det_dict["beam"]])
    dt_max = max([(beam_j["time"][-1]) for beam_j in det_dict["beam"]])

    t_min = t0 + np.timedelta64(int(dt_min * 1000.0), 'ms')
    t_max = t0 + np.timedelta64(int(dt_max * 1000.0), 'ms')

    fig = plt.figure(figsize=(8, 1 + 1.5 * param_set_cnt), layout="constrained")
    spec = fig.add_gridspec(param_set_cnt, 1)

    ax0 = fig.add_subplot(spec[param_set_cnt - 1])

    ax0.set_xlabel("")
    ax0.tick_params(axis='x', labelrotation=30)
    ax0.set_ylabel("Pressure [Pa]")

    ax0.set_xlim((t_min, t_max))

    beam_0 = det_dict["beam"][0]
    t_vals = [t0 + np.timedelta64(int(dt * 1000.0), 'ms') for dt in beam_0["time"]]

    ax0.fill_between(t_vals, -np.array(beam_0['resid']), beam_0['resid'], color='r', alpha=0.5)
    ax0.plot(t_vals, beam_0['signal'], 'k', linewidth=0.5)

    t1 = t0 + np.timedelta64(int(det_dict["start/end"][0][0] * 1000.0), 'ms')
    t2 = t0 + np.timedelta64(int(det_dict["start/end"][0][1] * 1000.0), 'ms')
    ax0.axvspan(t1, t2, color="lightsteelblue", alpha=0.75)

    if annotate_option == "frequency":
        ann_text = str(det_dict["fk_params"][0]["freq_min"]) + " - " + str(det_dict["fk_params"][0]["freq_max"]) + " Hz"
    else:
        ann_text = "Parameter Index: 0"

    ax0.annotate(ann_text, (0.975, 0.95), xycoords='axes fraction', horizontalalignment = "right", verticalalignment="top")

    for k, beam_k in enumerate(det_dict["beam"][1:]):
        ax_k = fig.add_subplot(spec[param_set_cnt - (k + 2)], sharex=ax0)
        plt.setp(ax_k.get_xticklabels(), visible=False)

        t_vals = [t0 + np.timedelta64(int(dt * 1000.0), 'ms') for dt in beam_k["time"]]
        ax_k.fill_between(t_vals, -np.array(beam_k['resid']), beam_k['resid'], color='r', alpha=0.5)
        ax_k.plot(t_vals, beam_k['signal'], 'k', linewidth=0.5)

        t1 = t0 + np.timedelta64(int(det_dict["start/end"][k + 1][0] * 1000.0), 'ms')
        t2 = t0 + np.timedelta64(int(det_dict["start/end"][k + 1][1] * 1000.0), 'ms')
        ax_k.axvspan(t1, t2, color="lightsteelblue", alpha=0.75)

        if annotate_option == "frequency":
            ann_text = str(det_dict["fk_params"][k + 1]["freq_min"]) + " - " + str(det_dict["fk_params"][k + 1]["freq_max"]) + " Hz"
        else:
            ann_text = "Parameter Index: " + str(k + 1)

        ax_k.annotate(ann_text, (0.975, 0.95), xycoords='axes fraction', horizontalalignment = "right", verticalalignment="top")

    if output_path:
        plt.savefig(output_path, dpi=250) 

    if show_fig:
        plt.show()


def plot_sd_json(det_dict, log_scale_freq=False, output_path=None, show_fig=True):
    '''
    Visualize multiple spectral detection (sd) results

    '''  

    f, t, Sxx_log = det_dict['spectrogram']
    f = np.array(f)
    t = np.array(t)
    Sxx_log = np.array(Sxx_log)

    # Build normalized spectrogram
    peaks_history, thresh_history, times_history = det_dict['history']   
    freq_band_mask = np.logical_and(det_dict["sd_params"]["freq_min"] < np.array(f), np.array(f) < det_dict["sd_params"]["freq_max"])

    Sxx_norm = np.array([Sxx_log[:, tn] - peaks_history[np.argmin(abs(tn - np.array(times_history)))][::2] for tn, n in enumerate(t)]).T
    Sxx_norm[(~freq_band_mask)] = np.nan      

    _, a = plt.subplots(3, sharex=True, figsize=(8, 8))

    try:
        Sxx_mean = np.mean(Sxx_log)
        Sxx_std = np.std(Sxx_log)
        Sxx_min = np.min(Sxx_log)
        Sxx_max = np.max(Sxx_log)
    except:
        Sxx_mean = np.mean(Sxx_log[Sxx_log != -np.inf])
        Sxx_std = np.std(Sxx_log[Sxx_log != -np.inf])
        Sxx_min = np.min(Sxx_log[Sxx_log != -np.inf])
        Sxx_max = np.max(Sxx_log[Sxx_log != -np.inf])

    cmap_max = min(Sxx_mean + 2.0 * Sxx_std, Sxx_max)
    cmap_min = max(Sxx_mean - 2.0 * Sxx_std, Sxx_min)

    a[2].set_xlabel("Time (rel. " + det_dict["wvfrm_info"][0]["starttime"] + " [s]")

    a[0].set_ylabel("Frequency [Hz]")
    a[1].set_ylabel("Frequency [Hz]")
    a[2].set_ylabel("Frequency [Hz]")
  
    a[1].sharex(a[0])
    a[1].sharey(a[0])
    a[2].sharex(a[0])
    a[2].sharey(a[0])

    a[0].imshow(np.flipud(Sxx_log), extent=[t[0], t[-1], f[1], f[-1]], cmap=cm.jet, aspect='auto', vmin=cmap_min, vmax=cmap_max)
    a[1].imshow(np.flipud(Sxx_norm), extent=[t[0], t[-1], f[0], f[-1]], cmap=cm.gnuplot2, aspect='auto', vmin=0.0, vmax= min(np.nanmax(Sxx_norm), 2.0 * np.nanstd(Sxx_norm)))



    for k, det in enumerate(det_dict['det_info']):
        spec_pnts = np.array(det['spec pnts'])
        dt = UTCDateTime(det['peak f-stat time']) - UTCDateTime(det_dict["wvfrm_info"][0]["starttime"])
        a[2].plot(spec_pnts[:, 0] + dt, spec_pnts[:, 1], '.', markersize=1.5)
        a[2].annotate(str(k), (dt, max(spec_pnts[:, 1])), xytext=(0, 10), textcoords='offset pixels', horizontalalignment = "center")
        a[2].grid()

    if log_scale_freq:
        a[0].set_yscale('log')
        a[1].set_yscale('log')
        a[2].set_yscale('log')

    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()


def plot_sd_single_json(det_dict, single_det_index, log_scale_freq=False, output_path=None, show_fig=True):
    '''
    Visualize a single spectral detection (sd) result

    '''  

    det_info = det_dict["det_info"][single_det_index]

    f, t, Sxx_log = det_info['spectrogram']
    try:
        Sxx_mean = np.mean(Sxx_log)
        Sxx_std = np.std(Sxx_log)
        Sxx_min = np.min(Sxx_log)
        Sxx_max = np.max(Sxx_log)
    except:
        Sxx_mean = np.mean(Sxx_log[Sxx_log != -np.inf])
        Sxx_std = np.std(Sxx_log[Sxx_log != -np.inf])
        Sxx_min = np.min(Sxx_log[Sxx_log != -np.inf])
        Sxx_max = np.max(Sxx_log[Sxx_log != -np.inf])

    cmap_max = min(Sxx_mean + 2.0 * Sxx_std, Sxx_max)
    cmap_min = max(Sxx_mean - 2.0 * Sxx_std, Sxx_min)

    spec_pnts = np.array(det_info['spec pnts'])
    t1, t2 = min(spec_pnts[:, 0]), max(spec_pnts[:, 0])
    f1, f2 = min(spec_pnts[:, 1]), max(spec_pnts[:, 1])

    fig = plt.figure(figsize=(10, 5), layout="constrained")
    spec = fig.add_gridspec(3, 6)

    ax0 = fig.add_subplot(spec[0, :4])
    ax0.set_xlim([t[0], t[-1]])
    ax0.plot(det_info["waveform"][0], det_info["waveform"][1], '-k', linewidth=0.5)
    ax0.axvspan(t1, t2, color='lightsteelblue', alpha=0.5)
    ax0.set_xlabel("")
    ax0.set_ylabel("Press. [Pa]")
    ax0.grid()

    ax1 = fig.add_subplot(spec[1, :4], sharex=ax0)
    ax1.imshow(np.flipud(Sxx_log), extent=[t[0], t[-1], f[1], f[-1]], cmap=cm.jet, aspect='auto', vmin=cmap_min, vmax=cmap_max)
    ax1.set_ylabel("Frequency [Hz]")

    if log_scale_freq:
        ax1.set_yscale('log')

    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax2 = fig.add_subplot(spec[2, :4], sharex=ax1, sharey=ax1)
    ax2.plot(np.array(det_info['spec pnts'])[:, 0], np.array(det_info['spec pnts'])[:, 1], '.', markersize=2.0, color='black')
    ax2.set_xlabel("Time (rel. " + det_info['peak f-stat time'] + ") [s]")
    ax2.set_ylabel("Frequency [Hz]")
    ax2.grid()

    ax3 = fig.add_subplot(spec[1:, 4:])
    ax3.semilogx(det_info['bg spec'][0], det_info['bg spec'][1], '-k', linewidth=1.0, label="Background")
    ax3.semilogx(det_info['bg spec'][0], det_info['bg spec'][2], '--k', linewidth=1.0, label="Threshold")
    ax3.semilogx(det_info['spec'][0], det_info['spec'][1], '-r', linewidth=1.5, label="Detection (mean)")
    ax3.semilogx(det_info['spec'][0], det_info['spec'][2], '--r', linewidth=1.0, label="Detection (peak)")
    ax3.axvspan(f1, f2, color='lightsteelblue', alpha=0.5, edgecolor=None)
    ax3.grid()
    ax3.set_ylim(-90.0)

    ax3.set_xlabel("Frequency [Hz]")
    ax3.set_ylabel("Power Spectral Density [Pa^2/Hz]")
    ax3.yaxis.set_label_position("right")
    ax3.yaxis.set_ticks_position("right")

    ax3.legend(loc="lower left", prop={'size': 8})

    ax4 = fig.add_subplot(spec[:2, 4:])
    ax4.axis('off')

    det_summary = 'Detection Summary\n' + '-' * 30
    det_summary = det_summary + '\nTrace ID: ' + det_dict["wvfrm_info"][0]["trace id"]
    det_summary = det_summary + '\nDate/Time: ' + det_info["peak f-stat time"]
    det_summary = det_summary + '\nDuration: ' + str(t2 - t1) + " sec"
    det_summary = det_summary + '\nFrequency Band: ' + str(f1) + " - " + str(f2) + " Hz"
    
    ax4.text(0.0, 1.0, det_summary, va="top", fontsize=11, bbox=dict(boxstyle="round, pad=0.2", fc="lightsteelblue", ec="black", lw=1))

    if output_path:
        plt.savefig(output_path, dpi=250) 

    if show_fig:
        plt.show()






##########################################
## THE REST OF THESE ARE DEPRECATED AND ## 
##  WILL BE REMOVED IN A FUTURE UPDATE  ##
##########################################






def plot_fk1(stream, latlon, times, peaks, detections=None, title=None, output_path=None, show_fig=True, det_thresh=None):
    '''
    Visualize beamforming (fk) results with waveform data included

    '''   
    x, t, t0, _ = beamforming_new.stream_to_array_data(stream, latlon)

    f, a = plt.subplots(4, figsize=(10, 6), sharex=True)
    a[3].set_xlabel("Time")
    a[0].set_ylabel("F-stat")    
    a[1].set_ylabel("Tr. Vel. [m/s]")
    a[2].set_ylabel("Back Az. [deg.]")
    a[3].set_ylabel("Pr. [Pa]")

    a[3].plot(np.array([t0 + np.timedelta64(int(tn * 1000.0), 'ms') for tn in t]), x[0,:], '-k')
    a[2].plot(times, peaks[:, 0], '.k', markersize=4)
    a[1].plot(times, peaks[:, 1], '.k', markersize=4)
    a[0].plot(times, peaks[:, 2], '.k', markersize=4)

    if detections:
        for det in detections:
            t1 = det.peakF_UTCtime + np.timedelta64(int(det.start * 1000.0), 'ms')
            t2 = det.peakF_UTCtime + np.timedelta64(int(det.end * 1000.0), 'ms')
            for n in range(4):
                a[n].axvspan(t1, t2, color="steelblue")

    if det_thresh is not None:
        a[0].plot(det_thresh[0], det_thresh[1], '--k', linewidth=0.5)

    if title:
        a[0].set_title(title)

    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()

def plot_fk2(times, peaks, detections=None, title=None, output_path=None, show_fig=True):
    '''
    Visualize beamforming (fk) results without waveform data
    '''

    f, a = plt.subplots(3, figsize=(10, 6), sharex=True)
    a[2].set_xlabel("Time")
    a[0].set_ylabel("F-stat")    
    a[1].set_ylabel("Tr. Vel. [m/s]")
    a[2].set_ylabel("Back Az. [deg.]")

    a[2].plot(times, peaks[:, 0], '.k', markersize=4)
    a[1].plot(times, peaks[:, 1], '.k', markersize=4)
    a[0].plot(times, peaks[:, 2], '.k', markersize=4)

    if detections:
        for det in detections:
            t1 = det.peakF_UTCtime + np.timedelta64(int(det.start * 1000.0), 'ms')
            t2 = det.peakF_UTCtime + np.timedelta64(int(det.end * 1000.0), 'ms')
            for n in range(3):
                a[n].axvspan(t1, t2, color="steelblue")

    if title:
        a[0].set_title(title)
       
    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()


def plot_sd(trace, det_list, freq_band, spec_option="spectrogram", morlet_omega0=12.0, title=None, output_path=None, show_fig=False):
    '''
    Visualize multiple spectral detection (sd) results

    '''  
    trace.filter('bandpass', freqmin=freq_band[0], freqmax=freq_band[1])

    dt = trace.stats.delta
    nperseg = int((4.0 / freq_band[0]) / dt) 
    if spec_option == "spectrogram":
        f, t, Sxx = spectrogram(trace.data, 1.0 / dt, nperseg=nperseg, noverlap=int(nperseg * 0.8))
        Sxx_log = 10.0 * np.log10(Sxx)
    elif spec_option == "stft":
        f, t, Sxx = stft(trace.data, 1.0 / dt, nperseg=nperseg, noverlap=int(nperseg * 0.8))
        Sxx_log = 10.0 * np.log10(abs(Sxx))
    elif spec_option == "cwt":
        f, _, _ = spectrogram(trace.data, 1.0 / dt, nperseg=nperseg, noverlap=int(nperseg * 0.8))
        t = trace.times()       
        widths = morlet_omega0 / (2 * np.pi * f) * (1.0 / dt)
        Sxx_log = 10.0 * np.log10(abs(cwt(trace.data, morlet2, widths, w=morlet_omega0)))

    fig, a = plt.subplots(3, sharex=True, figsize=(9, 5))
    a[0].plot(trace.times(), trace.data, '-k')
    a[0].set_xlim(trace.times()[0], trace.times()[-1])

    cmap_max = np.mean(Sxx_log[Sxx_log != -np.inf]) + 2.0 * np.std(Sxx_log[Sxx_log != -np.inf])
    cmap_min = np.mean(Sxx_log[Sxx_log != -np.inf]) - 2.0 * np.std(Sxx_log[Sxx_log != -np.inf])
    a[1].imshow(np.flipud(Sxx_log), extent=[t[0], t[-1], f[1], f[-1]], cmap=cm.jet, aspect='auto', vmin=cmap_min, vmax=cmap_max)

    a[2].set_xlabel("Time [s]")
    a[2].set_ylabel("Frequency [Hz]")
    a[1].set_ylabel("Frequency [Hz]")
    a[0].set_ylabel("Amplitude")

    a[1].axhline(freq_band[0], color='0.5')
    a[1].axhline(freq_band[1], color='0.5')

    a[2].sharey(a[1])
    a[2].set_ylim(freq_band[0], freq_band[1])
    for det in det_list:
        Sxx_pnts = np.array(det['Sxx_points'])
        a[2].plot(Sxx_pnts[:, 0], Sxx_pnts[:, 1], '.', markersize=1.5)

    a[1].set_yscale('log')
    a[2].set_yscale('log')

    if title:
        a[0].set_title(title)
       
    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()


def plot_sd_single(trace, det_info, freq_band, spec_option="spectrogram", morlet_omega0=12.0, title=None, output_path=None, show_fig=False):
    '''
    Visualize a single spectral detection (sd) result

    '''  
    t_shift = UTCDateTime(det_info['Time (UTC)']) - UTCDateTime(trace.stats.starttime)

    dt = trace.stats.delta
    nperseg = int((4.0 / freq_band[0]) / dt) 
    if spec_option == "spectrogram":
        f, t, Sxx = spectrogram(trace.data, 1.0 / dt, nperseg=nperseg, noverlap=int(nperseg * 0.8))
        Sxx_log = 10.0 * np.log10(Sxx)
    elif spec_option == "stft":
        f, t, Sxx = stft(trace.data, 1.0 / dt, nperseg=nperseg, noverlap=int(nperseg * 0.8))
        Sxx_log = 10.0 * np.log10(abs(Sxx))
    elif spec_option == "cwt":
        f, _, _ = spectrogram(trace.data, 1.0 / dt, nperseg=nperseg, noverlap=int(nperseg * 0.8))
        t = trace.times()       
        widths = morlet_omega0 / (2 * np.pi * f) * (1.0 / dt)
        Sxx_log = 10.0 * np.log10(abs(cwt(trace.data, morlet2, widths, w=morlet_omega0)))

    Sxx_pnts = np.array(det_info['Sxx_points'])

    fig = plt.figure(figsize=(10, 5), layout="constrained")
    spec = fig.add_gridspec(3, 5)

    ax1 = fig.add_subplot(spec[1, :3])
    cmap_max = np.mean(Sxx_log[Sxx_log != -np.inf]) + 2.0 * np.std(Sxx_log[Sxx_log != -np.inf])
    cmap_min = np.mean(Sxx_log[Sxx_log != -np.inf]) - 2.0 * np.std(Sxx_log[Sxx_log != -np.inf])   
    ax1.imshow(np.flipud(Sxx_log), extent=[t[0] - t_shift, t[-1] - t_shift, f[1], f[-1]], cmap=cm.jet, aspect='auto', vmin=cmap_min, vmax=cmap_max)
    
    ax1.set_xlim([det_info['Start'] - 150.0, det_info['End'] + 150.0])
    ax1.axhline(freq_band[0], color='0.5')
    ax1.axhline(freq_band[1], color='0.5')
    ax1.set_yscale('log')
    ax1.set_ylabel("Frequency [Hz]")

    ax0 = fig.add_subplot(spec[2, :3], sharex=ax1, sharey=ax1)
    ax0.plot(Sxx_pnts[:, 0] - t_shift, Sxx_pnts[:, 1], '.', markersize=1.5, color='black')
    ax0.set_yscale('log')
    ax0.set_xlabel("Time (rel. " + det_info['Time (UTC)'] + ") [s]")
    ax0.set_ylabel("Frequency [Hz]")
    
    trace.filter('bandpass', freqmin=det_info['Freq Range'][0], freqmax=det_info['Freq Range'][1])

    ax2 = fig.add_subplot(spec[0, :3], sharex=ax1)
    ax2.plot(trace.times() - t_shift, trace.data, '-k')
    ax2.axvspan(det_info['Start'], det_info['End'], color='green', alpha=0.5)
    ax2.set_xlabel("")
    ax2.set_ylabel("Press. [Pa]")

    ax3 = fig.add_subplot(spec[1:, 3:])
    ax3.semilogx(det_info['Background Peaks'][0], det_info['Background Peaks'][1], '-k', linewidth=1.0, label="Background (Mean)")
    ax3.semilogx(det_info['Background Threshold'][0], det_info['Background Threshold'][1], '--k', linewidth=1.0, label="Detection Threshold")
    ax3.semilogx(det_info['Sxx_det_mean'][0], det_info['Sxx_det_mean'][1], '-b', linewidth=1.5, label="Detection (mean)")
    ax3.semilogx(det_info['Sxx_det_max'][0], det_info['Sxx_det_max'][1], '-r', linewidth=1.0, label="Detection (max)")
    ax3.axvspan(det_info['Freq Range'][0], det_info['Freq Range'][1], color='green', alpha=0.5, edgecolor=None)
    ax3.set_xlabel("Frequency [Hz]")
    ax3.set_ylabel("Power Spectral Density [Pa^2/Hz]")
    ax3.yaxis.set_label_position("right")
    ax3.yaxis.set_ticks_position("right")

    fig.legend(loc="upper right")

    if title:
        ax2.set_title(title)
       
    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()