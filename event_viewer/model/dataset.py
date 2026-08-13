"""
A per-file dataset: the event table plus on-demand single-event access.

``EventDataset`` ties together an :class:`EventFileReader` (raw I/O) and a
:class:`DetectorModel` (geometry), exposing exactly what the controller and
visualisers need: the per-event DataFrame, the indices passing a cut, and a fully
built :class:`Event` (per-hit arrays with physical z + the row of metrics).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .event import Event


class EventDataset:
    """Per-file view combining the reader's table with single-event access."""

    def __init__(self, reader, detector):
        self.reader = reader
        self.detector = detector
        self._table = reader.event_table(0.0)

    # -------------------------------------------------------------- table --
    @property
    def table(self) -> pd.DataFrame:
        return self._table

    @property
    def n_events(self) -> int:
        return self.reader.n_events

    @property
    def has_metrics(self) -> bool:
        return self.reader.has_metrics

    def feature_columns(self):
        return self.reader.feature_columns()

    def passing(self, cut_model) -> np.ndarray:
        """Tree-entry indices passing ``cut_model`` (all events if empty)."""
        if cut_model is None or cut_model.is_empty:
            return np.arange(self.n_events)
        return cut_model.passing_indices(self._table)

    # ------------------------------------------------------------- event --
    def get_event(self, index: int, want_tracks: bool = False) -> Event:
        """Build the :class:`Event` at tree entry ``index``.

        ``want_tracks`` fetches the reconstructed ACTSTracks too (a podio
        Frame access, same cost class as the hits) -- skipped by default
        since most renders don't draw them (see ``show-overlays``'s
        "tracks" checkbox).
        """
        hits = self.reader.read_hits(index)
        slab = np.asarray(hits.get("hit_slab", []), dtype=np.int64)
        row = self._table.iloc[index].to_dict() if index < len(self._table) else {}
        tracks = []
        if want_tracks:
            read_tracks = getattr(self.reader, "read_tracks", None)
            if read_tracks:
                tracks = self._reframe_tracks(read_tracks(index))
        return Event(
            index=index,
            x=np.asarray(hits.get("hit_x", []), dtype=float),
            y=np.asarray(hits.get("hit_y", []), dtype=float),
            z=self.detector.slab_z_array(slab),
            slab=slab,
            chip=np.asarray(hits.get("hit_chip", []), dtype=np.int64),
            chan=np.asarray(hits.get("hit_chan", []), dtype=np.int64),
            energy=np.asarray(hits.get("hit_energy", []), dtype=float),
            hg=np.asarray(hits.get("hit_hg", []), dtype=float),
            lg=np.asarray(hits.get("hit_lg", []), dtype=float),
            metrics=row,
            tracks=tracks,
        )

    def _reframe_tracks(self, tracks):
        """Remap track-state z into the hits' own frame, when the two disagree.

        Simulation's ``ACTSTracks`` come out in the DD4hep frame (job4_tracking
        runs on the pre-flip ``SiPadHitsDigi``, matching the DD4hep-built ACTS
        surfaces) while this same file's hits are already in the test-beam
        frame (``SiPadHitsMapped``/``ECalHits``) -- two different z
        conventions in one file, so a track drawn as-is lands on the wrong
        side of (or outside) the detector entirely. Real test-beam tracks have
        no such mismatch (that repo's ACTSGeoSvc builds its surfaces straight
        from the test-beam z table already), so this is a no-op there.

        Detected from :meth:`PidFileReader.track_z_table`'s own sign: the
        test-beam frame is always <= 0 (``mappings/slab_z_positions.yml``),
        the DD4hep one always > 0 -- remap only when it's positive. x/y are
        untouched (confirmed against a real digitized.edm4hep.root: only z
        differs between the pre- and post-flip hit collections). Matched by
        NEAREST value, not an exact/rounded lookup: ``track_z_table()``'s
        entries come from a ``np.round`` over a float32 branch, which can
        disagree with plain ``round()`` right at the 0.05 mm boundary (e.g.
        49.349998 -- physically 49.35 -- rounds to 49.3 in float64 but 49.4 in
        float32), silently dropping a few of the 15 slabs if matched by key.
        A nearest-within-tolerance match is immune to that.
        """
        if not tracks:
            return tracks
        dd_table, tb_table = self._track_z_tables()
        if dd_table is None:
            return tracks
        out = []
        for tr in tracks:
            pts = [(x, y, self._remap_z(z, dd_table, tb_table))
                  for x, y, z in tr["points"]]
            out.append({**tr, "points": pts})
        return out

    # Half the 15 mm layer pitch would already be unambiguous; well under it
    # to stay clear of the one 30 mm gap (layer 10->11, the empty rail slot).
    _Z_MATCH_TOL_MM = 5.0

    @classmethod
    def _remap_z(cls, z, dd_table, tb_table):
        i = int(np.argmin(np.abs(dd_table - z)))
        return float(tb_table[i]) if abs(dd_table[i] - z) < cls._Z_MATCH_TOL_MM else z

    def _track_z_tables(self):
        """Cached ``(dd4hep_z_array, test_beam_z_array)``, index-aligned by
        slab, or ``(None, None)`` if no remap is needed/possible (see
        :meth:`_reframe_tracks`)."""
        if not hasattr(self, "_track_z_tables_cache"):
            self._track_z_tables_cache = (None, None)
            get_table = getattr(self.reader, "track_z_table", None)
            dd_table = get_table() if get_table else None
            if dd_table is not None and len(dd_table) and dd_table[0] > 0:
                tb = np.asarray(self.detector.slab_z_mm, dtype=float)
                n = min(len(dd_table), len(tb))
                self._track_z_tables_cache = (
                    np.asarray(dd_table[:n], dtype=float), tb[:n])
        return self._track_z_tables_cache

    def accumulate(self, indices) -> Event:
        """Aggregate the hits of many events into one pseudo-:class:`Event`.

        Hits from all ``indices`` are pooled and summed per pad (same ``slab`` and
        rounded ``x``/``y``), so the resulting "event" shows the *accumulated*
        energy deposited in each pad across the whole set -- the average shower
        footprint of, e.g., a cluster. ``chip``/``chan`` are not meaningful for an
        aggregate and are set to ``-1``.
        """
        from .._timing import timed

        hits = self.reader.all_hits()
        indices = list(indices)
        if not indices or "hit_x" not in hits or "hit_y" not in hits:
            empty = np.empty(0)
            return Event(index=-1, x=empty, y=empty, z=empty,
                         slab=empty.astype(np.int64), chip=empty.astype(np.int64),
                         chan=empty.astype(np.int64), energy=empty,
                         hg=empty, lg=empty, metrics={})
        with timed("accumulate (pool + groupby per pad)") as info:
            slab = np.concatenate([np.asarray(hits["hit_slab"][i]) for i in indices])
            x = np.concatenate([np.asarray(hits["hit_x"][i]) for i in indices])
            y = np.concatenate([np.asarray(hits["hit_y"][i]) for i in indices])
            energy = np.concatenate([np.asarray(hits["hit_energy"][i]) for i in indices])

            pooled = pd.DataFrame({"slab": slab.astype(int),
                                   "xr": np.round(x, 2), "yr": np.round(y, 2),
                                   "e": energy.astype(float)})
            agg = pooled.groupby(["slab", "xr", "yr"], as_index=False)["e"].sum()
            agg["e"] /= max(len(indices), 1)
            slab_g = agg["slab"].to_numpy(dtype=np.int64)
            info["events"] = len(indices)
            info["hits_pooled"] = int(slab.size)
            info["pads"] = int(slab_g.size)
        return Event(
            index=-1,
            x=agg["xr"].to_numpy(dtype=float),
            y=agg["yr"].to_numpy(dtype=float),
            z=self.detector.slab_z_array(slab_g),
            slab=slab_g,
            chip=np.full(slab_g.size, -1, dtype=np.int64),
            chan=np.full(slab_g.size, -1, dtype=np.int64),
            energy=agg["e"].to_numpy(dtype=float),
            hg=np.empty(0), lg=np.empty(0), metrics={},
        )
