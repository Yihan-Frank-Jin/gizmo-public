#!/usr/bin/env python3
"""Check Yuan18 wind mass accounting uses direct reservoir feeding.

Yuan18 solves mdot_bh and mdot_wind explicitly, so the physical BH mass should
grow only by mdot_bh. The wind reservoir should be fed directly from mdot_wind
without temporarily adding wind mass to BPP.BH_Mass and subtracting it later.
"""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
BLACKHOLE_C = ROOT / "galaxy_sf" / "blackholes" / "blackhole.c"


def require_absent(source: str, pattern: str, description: str) -> list[str]:
    if re.search(pattern, source, re.MULTILINE):
        return [f"unexpected {description}"]
    return []


def require_present(source: str, pattern: str, description: str) -> list[str]:
    if not re.search(pattern, source, re.MULTILINE | re.DOTALL):
        return [f"missing {description}"]
    return []


def main() -> int:
    source = BLACKHOLE_C.read_text()
    errors: list[str] = []

    errors += require_absent(
        source,
        r"BPP\(n\)\.BH_Mass\s*\+=\s*BlackholeTempInfo\[i\]\.Yuan18_mdot_wind\s*\*\s*dt",
        "temporary Yuan18 wind over-accretion into BPP(n).BH_Mass",
    )
    errors += require_absent(
        source,
        r"BPP\(n\)\.BH_Mass\s*-=\s*dm_wind_yuan18",
        "Yuan18 wind subtraction from BPP(n).BH_Mass during reservoir feed",
    )
    errors += require_present(
        source,
        r"double\s+dm_wind_yuan18\s*=\s*BlackholeTempInfo\[i\]\.Yuan18_mdot_wind\s*\*\s*dt;",
        "direct Yuan18 wind reservoir mass from mdot_wind * dt",
    )
    errors += require_present(
        source,
        r"BPP\(n\)\.BH_AccretionDeficit\s*\+=\s*BlackholeTempInfo\[i\]\.Yuan18_mdot_wind\s*\*\s*dt;",
        "Yuan18 wind contribution to BH_AccretionDeficit",
    )

    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1

    print("PASS: Yuan18 wind mass accounting feeds reservoir without BH_Mass over-accretion.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
