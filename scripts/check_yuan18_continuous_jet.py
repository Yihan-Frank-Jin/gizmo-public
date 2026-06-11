#!/usr/bin/env python3
"""Static regression checks for Yuan18 continuous jet feedback wiring."""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text()


def require(condition, message, failures):
    if not condition:
        failures.append(message)


def contains(path, pattern, failures, message=None):
    text = read(path)
    ok = re.search(pattern, text, re.MULTILINE | re.DOTALL) is not None
    require(ok, message or f"{path} missing pattern: {pattern}", failures)


def main():
    failures = []

    contains(
        "Template_Config.sh",
        r"BH_YUAN18_JET_CONTINUOUS",
        failures,
        "Template_Config.sh must document BH_YUAN18_JET_CONTINUOUS next to Yuan18 continuous feedback flags.",
    )
    contains(
        "allvars.h",
        r"defined\(BH_YUAN18_JET_CONTINUOUS\).*?!defined\(BH_YUAN18_ACCRETION\)",
        failures,
        "allvars.h must require BH_YUAN18_JET_CONTINUOUS to pair with BH_YUAN18_ACCRETION.",
    )
    contains("galaxy_sf/blackholes/blackhole.h", r"MyFloat\s+Yuan18_mdot_jet\s*;", failures)
    contains("galaxy_sf/blackholes/blackhole.h", r"MyFloat\s+Yuan18_v_jet\s*;", failures)
    contains("galaxy_sf/blackholes/blackhole.h", r"MyFloat\s+Yuan18_eps_jet\s*;", failures)
    contains("galaxy_sf/blackholes/blackhole.c", r"mdot_jet\s*=\s*0\.5\s*\*\s*mdot_bh\s*;", failures)
    contains("galaxy_sf/blackholes/blackhole.c", r"v_jet\s*=\s*0\.3\s*\*\s*C_LIGHT_CODE\s*;", failures)
    contains("galaxy_sf/blackholes/blackhole.c", r"eps_jet\s*=\s*0\s*;", failures)
    contains("galaxy_sf/blackholes/blackhole.c", r"mode_wind\s*==\s*4.*YUAN18_COS_ANG_JET", failures)
    contains(
        "galaxy_sf/blackholes/blackhole_swallow_and_kick.c",
        r"Yuan18_mdot_jet\s*\*\s*local\.Dt",
        failures,
        "blackhole_swallow_and_kick.c must couple continuous jet mass using Yuan18_mdot_jet * local.Dt.",
    )
    for suffix, dataset in (
        ("MASS", "Yuan18JetMass"),
        ("ENERGY", "Yuan18JetEnergy"),
        ("MOMENTUM", "Yuan18JetMomentum"),
        ("LASTMODE", "Yuan18JetLastMode"),
    ):
        contains("allvars.h", rf"IO_YUAN18_JET_{suffix}", failures)
        contains("io.c", rf"case\s+IO_YUAN18_JET_{suffix}:", failures)
        contains("io.c", rf"strcpy\(buf,\s*\"{dataset}\"\)", failures)
        contains("read_ic.c", rf"case\s+IO_YUAN18_JET_{suffix}:", failures)

    if failures:
        print("Yuan18 continuous jet check failed:")
        for failure in failures:
            print(f" - {failure}")
        return 1

    print("Yuan18 continuous jet check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
