#!/usr/bin/env python3
"""Regression check for the Yuan18 wind dispatcher gate in run.c.

The run-loop should not enter the Yuan18 spawn dispatcher merely because a BH
has a non-empty reservoir. The dispatcher is useful only once the largest
global reservoir contains at least YUAN18_WIND_N_MIN target wind particles, and
particle rearrangement is required only after particles were actually spawned.
"""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
RUN_C = ROOT / "run.c"


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    sys.exit(1)


text = RUN_C.read_text()
match = re.search(
    r"#ifdef BH_YUAN18_WIND\b(?P<block>.*?)#endif",
    text,
    flags=re.S,
)
if not match:
    fail("could not locate BH_YUAN18_WIND block in run.c")

block = match.group("block")

if not re.search(
    r"if\s*\(\s*Max_Yuan18_WindReservoirMassUnits_fromSink_global\s*>=\s*YUAN18_WIND_N_MIN\s*\)",
    block,
):
    fail("Yuan18 wind dispatcher is not gated by YUAN18_WIND_N_MIN")

spawn_then_rearrange = re.search(
    r"yuan18_wind_particles_spawned\s*=\s*spawn_bh_yuan18_wind_feedback\s*\(&yuan18_wind_mass_spawned\)\s*;"
    r"(?P<after_spawn>.*?)"
    r"rearrange_particle_sequence\s*\(\s*\)\s*;",
    block,
    flags=re.S,
)
if not spawn_then_rearrange:
    fail("could not locate Yuan18 spawn call followed by particle rearrangement")

if not re.search(
    r"if\s*\(\s*yuan18_wind_particles_spawned\s*>\s*0\s*\)",
    spawn_then_rearrange.group("after_spawn"),
):
    fail("Yuan18 wind rearrangement is not guarded by actual spawned particles")

if 'printf("[Yuan18-wind]' in block and not re.search(
    r"if\s*\(\s*\(\s*ThisTask\s*==\s*0\s*\)\s*&&\s*\(\s*yuan18_wind_particles_spawned\s*>\s*0\s*\)\s*\)",
    block,
):
    fail("Yuan18 wind debug output is not limited to real spawn events")

print("PASS: Yuan18 wind dispatch waits for a full batch and rearranges only after spawning.")
