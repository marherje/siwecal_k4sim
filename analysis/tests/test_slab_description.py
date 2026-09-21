"""
The per-slab technology description (mappings/slab_z_positions.yml) and its
loader, analysis/slab_description.py.

The description is what replaced the hard-coded "slab 12 is the chip-on-board"
in job3 / ChannelMapper / the tests, so these pin the facts the rest of the
chain relies on: slab 12 is the FEV11 COB with its own pad map and its own
threshold DAC, slab 14 carries the 650 um sensor, and a file without the block
still describes the old all-FEV10 detector.
"""

from __future__ import annotations

import os
import sys
import textwrap

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO_ROOT)

from analysis.slab_description import (  # noqa: E402
    DEFAULT_YAML, N_SLABS, load_slab_description)


@pytest.fixture(scope="module")
def slabs():
    return load_slab_description(DEFAULT_YAML)


def test_fifteen_slabs_everywhere(slabs):
    for name in ("technology", "sensor_thickness_um", "pad_map", "threshold_dac",
                 "slab_z_mm", "w_thickness_mm"):
        assert len(getattr(slabs, name)) == N_SLABS, name


def test_slab_12_is_the_chip_on_board(slabs):
    assert slabs.has_technology_block
    assert slabs.technology[12] == "FEV11_COB"
    assert all(t == "FEV10" for s, t in enumerate(slabs.technology) if s != 12)
    assert os.path.basename(slabs.pad_map[12]).startswith("fev11_cob")
    assert slabs.pad_map_overrides() == {12: slabs.pad_map[12]}
    assert os.path.isfile(slabs.pad_map[12])
    assert os.path.isfile(slabs.default_pad_map())


def test_slab_12_has_its_own_threshold_dac(slabs):
    # Run book: "Applied 220 to all slabs, all Asics but slab 12 (at 243)".
    assert slabs.threshold_dac[12] == 243
    assert slabs.slabs_with_threshold_override() == {12: 243}


def test_sensor_thickness(slabs):
    assert slabs.sensor_thickness_um[14] == 650
    assert all(t == 500 for s, t in enumerate(slabs.sensor_thickness_um) if s != 14)


def test_file_without_the_block_is_the_old_detector(tmp_path):
    p = tmp_path / "slabs.yml"
    p.write_text(textwrap.dedent("""
        slab_z_mm: [0, -15, -30, -45, -60, -75, -90, -105, -120, -135, -150, -180, -195, -210, -225]
        w_thickness_mm: [2.8, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6]
    """))
    d = load_slab_description(str(p), mappings_dir=os.path.join(REPO_ROOT, "mappings"))
    assert not d.has_technology_block
    assert d.pad_map_overrides() == {}
    assert d.slabs_with_threshold_override() == {}
    assert set(d.sensor_thickness_um) == {500}
    assert set(d.technology) == {"FEV10"}


def test_wrong_length_is_an_error(tmp_path):
    p = tmp_path / "slabs.yml"
    p.write_text("slab_technology: [FEV10, FEV10]\n")
    with pytest.raises(ValueError):
        load_slab_description(str(p))


def test_undeclared_technology_is_an_error(tmp_path):
    p = tmp_path / "slabs.yml"
    p.write_text("slab_technology: [" + ", ".join(["FEV10"] * 12 + ["FEV99", "FEV10", "FEV10"]) + "]\n")
    with pytest.raises(ValueError):
        load_slab_description(str(p))


def test_both_repos_carry_the_same_file():
    twin = os.path.join(REPO_ROOT, "..", "siwecal-tb2026", "mappings", "slab_z_positions.yml")
    if not os.path.isfile(twin):
        pytest.skip("siwecal-tb2026 checkout not next to this repo")
    assert open(DEFAULT_YAML).read() == open(twin).read()
