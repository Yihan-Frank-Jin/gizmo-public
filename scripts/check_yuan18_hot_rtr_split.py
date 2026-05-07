#!/usr/bin/env python3
"""Regression check for Yuan18 hot-mode radius bookkeeping."""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "galaxy_sf" / "blackholes" / "blackhole.c"


def fail(message):
    print(f"FAIL: {message}", file=sys.stderr)
    sys.exit(1)


text = SOURCE.read_text()
match = re.search(
    r"else if \(mdot_bondi > MIN_REAL_NUMBER\).*?else // mdot_bondi <= mdot_crit",
    text,
    flags=re.S,
)
if not match:
    fail("could not locate Yuan18 hot-mode block")

hot = match.group(0)

required_fragments = [
    "r_tr_physical",
    "r_tr_mdot",
    "r_tr_mdot = DMIN(r_tr_physical, x1min)",
    "r_inject = DMAX(r_tr_physical, x1min)",
    "mdot_bh = mdot_bondi * sqrt(3. * r_s / r_tr_mdot)",
    "v_wind = 0.2 * sqrt(All.G * BPP(n).BH_Mass / r_tr_physical)",
    "pow(r_inject / r_tr_physical",
]

for fragment in required_fragments:
    if fragment not in hot:
        fail(f"missing expected fragment: {fragment}")

print("Yuan18 hot-mode radius split is present")
