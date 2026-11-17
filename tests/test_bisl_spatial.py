"""BISL regressions using synthetic detections and the real library/CLI paths."""
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from obspy import UTCDateTime
from scipy.integrate import simpson

from infrapy.event import bisl
from infrapy.propagation import infrasound, likelihoods
from infrapy.utils import data_io


class SpatialLocalizationTests(unittest.TestCase):
    def setUp(self):
        for name in ("canon_rcel_wts", "canon_rcel_mns", "canon_rcel_vrs"):
            self.addCleanup(setattr, infrasound, name, getattr(infrasound, name))
        infrasound._load_celerity_model("regional_lf")
        self.time_limits = tuple(np.datetime64(t) for t in (
            "2026-11-09T23:59:30", "2026-11-10T00:00:30"))

    def detection_info(self, east=False, timed=False):
        stations = [(-35.1, 117.0 if east else -117.0, 31.0),
                    (-34.9, 116.8 if east else -117.2, 90.0)]
        return [{
            "peak f-stat time": f"2026-11-10T00:01:{30 + 20 * n:02}" if timed else None,
            "start/end": [[-5.0, 5.0]], "f-stat": 8.0,
            "back az": azimuth, "tr vel": 340.0,
            "wvfrm_info": [[{"latitude": lat, "longitude": lon,
                              "trace id": f"XX.A{n}{element}..BDF"}
                             for element in range(3)]],
            "fk_params": [{"freq_min": 0.5, "freq_max": 2.0}],
        } for n, (lat, lon, azimuth) in enumerate(stations)]

    def options(self, east=False):
        return dict(bm_width=10.0, rng_max=1000.0, grid_resol=3,
                    ll_corner=(-35.0, 117.0 if east else -117.0),
                    ur_corner=(-34.8, 117.3 if east else -116.7),
                    latlon_resol=0.1, tm_lims=self.time_limits, tm_resol=10.0)

    def prepare(self, east=False, timed=False):
        detections = [data_io._det_dict_to_likelihood(info)
                      for info in self.detection_info(east, timed)]
        for detection in detections:
            self.assertIs(type(detection), likelihoods.InfrasoundDetection)
            self.assertEqual(bool(detection.peakF_UTCtime < UTCDateTime("9999-01-01")), timed)
        lat, lon, _, times = bisl.build_grid(
            detections, include_tms=timed, **self.options(east))
        self.assertEqual(lat.shape[:2], (3, 4))
        values = [det.pdf(lat, lon, times) if timed else det.az_pdf(lat, lon)
                  for det in detections]
        density = np.asarray(values).prod(axis=0)
        self.assertTrue(np.isfinite(density).all())
        self.assertTrue(np.any(density > 0))
        return detections, lat, lon, times, density

    def run_spatial(self, east=False):
        detections, lat, lon, _, density = self.prepare(east)
        result = bisl.run(detections, verbose=False, show_prog=False, **self.options(east))
        return result, lat, lon, density

    def assert_maximum(self, result, lat, lon, density):
        points = list(zip(lon.flat, lat.flat, density.flat))
        expected_lon, expected_lat, expected_density = max(points, key=lambda p: p[2])
        self.assertEqual(result["MaP_val"], expected_density)
        self.assertEqual(result["lat_MaP"], expected_lat)
        self.assertEqual(result["lon_MaP"], expected_lon)

    def test_direction_only_maximum_in_western_grid(self):
        self.assert_maximum(*self.run_spatial())

    def test_direction_only_maximum_in_eastern_grid(self):
        self.assert_maximum(*self.run_spatial(east=True))

    def test_export_coordinates_stay_paired_with_raw_density(self):
        result, lat, lon, density = self.run_spatial(east=True)
        expected = list(zip(lon.flat, lat.flat, density.flat))
        self.assertEqual(np.asarray(result["spatial_pdf"]).shape, (3, len(expected)))
        self.assertEqual(list(map(tuple, np.asarray(result["spatial_pdf"]).T)), expected)

    def test_spatial_result_exposes_existing_normalization(self):
        result, lat, lon, density = self.run_spatial(east=True)
        expected = simpson(simpson(density * np.cos(np.radians(lat)),
                                   x=np.unique(lon)), x=np.unique(lat))
        self.assertIn("norm", result)
        self.assertEqual(result["norm"], expected)
        self.assertGreater(result["norm"], 0)
        self.assertNotIn("t_MaP", result)
        self.assertNotIn("temporal_pdf", result)

    def test_existing_timed_maximum_still_uses_probability(self):
        detections, lat, lon, times, density = self.prepare(timed=True)
        result = bisl.run(detections, verbose=False, show_prog=False, **self.options())
        self.assert_maximum(result, lat, lon, density)
        expected = max(zip(times.flat, density.flat), key=lambda p: p[1])[0]
        self.assertEqual(result["t_MaP"], expected)
        self.assertGreater(result["norm"], 0)
        self.assertTrue(bisl.summarize(result))

    def test_CLI_persists_and_reuses_one_spatial_result(self):
        from click.testing import CliRunner
        from infrapy.cli.__main__ import main

        _, lat, lon, _, density = self.prepare(east=True)
        original = {"det_info": self.detection_info(east=True), "location": [],
                    "assoc_params": {"range_max": 1000.0}}
        config = """[LOC]
back_az_width = 10.0
range_max = 1000.0
grid_resol = 3
ll_corner = -35.0,117.0
ur_corner = -34.8,117.3
latlon_resol = 0.1
tm_min = 2026-11-09T23:59:30
tm_max = 2026-11-10T00:00:30
tm_resol = 10.0
celerity_model = regional_lf
pgm_file = None
atmo_data = None
alt_lims = None
alt_resol = None
grnd_snd_spd = None
c0_stdev = None
det_tm_stdev = None
az_limit = None
local_temp_dir = None
cpu_cnt = 1
"""
        with tempfile.TemporaryDirectory() as directory:
            event_path = Path(directory) / "direction-only.ev.json.gz"
            config_path = Path(directory) / "loc.cnfg"
            config_path.write_text(config)
            with gzip.open(event_path, "wt", encoding="UTF-8") as stream:
                json.dump(original, stream)
            args = ["event", "locate", "--ev-file", str(event_path),
                    "--cnfg-file", str(config_path), "--range-max", "1000"]
            first = CliRunner().invoke(main, args)
            self.assertEqual(first.exit_code, 0, first.output + repr(first.exception))
            first_bytes = event_path.read_bytes()
            with gzip.open(event_path, "rt", encoding="UTF-8") as stream:
                saved = json.load(stream)
            self.assertEqual({k: v for k, v in saved.items() if k != "location"},
                             {k: v for k, v in original.items() if k != "location"})
            self.assertEqual(len(saved["location"]), 1)
            result = saved["location"][0]["result"]
            self.assert_maximum(result, lat, lon, density)
            self.assertEqual(list(map(tuple, result["spatial_pdf"])),
                             [tuple(lon.flat), tuple(lat.flat), tuple(density.flat)])
            self.assertGreater(result["norm"], 0)
            self.assertTrue(bisl.summarize(result))
            second = CliRunner().invoke(main, args)
            self.assertEqual(second.exit_code, 0, second.output + repr(second.exception))
            self.assertEqual(event_path.read_bytes(), first_bytes)


if __name__ == "__main__":
    unittest.main()
