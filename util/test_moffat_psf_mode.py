#!/usr/bin/env python3
"""Minimal regression check for SIMLIB PSF_MODEL: MOFFAT.

Runs `snlc_sim.exe NOFILE ... SIMLIB_DUMP 2` on a tiny SIMLIB and verifies:
- NEA ratio (MOFFAT/GAUSS) matches the analytic expectation at fixed FWHM.
- Central-pixel fraction (FPIX) is smaller for MOFFAT than GAUSS.

This test is intentionally lightweight and does not require large inputs.
It does require a valid SNANA environment (notably $SNDATA_ROOT).
"""

from __future__ import annotations

import math
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict


def nea_ratio_moffat_over_gauss(beta: float) -> float:
    """Return NEA_M/NEA_G at fixed FWHM for a normalized circular Moffat PSF."""
    return (
        0.25
        * (2.0 * math.log(2.0))
        / (2.0 ** (1.0 / beta) - 1.0)
        * (2.0 * beta - 1.0)
        / ((beta - 1.0) ** 2)
    )


def parse_dump_obs(path: Path) -> Dict[str, str]:
    varnames = None
    row = None
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("VARNAMES:"):
            varnames = line.split()[1:]
        if line.startswith("ROW:") and row is None:
            row = line.split()[1:]
            break
    if not varnames or not row:
        raise RuntimeError(f"Could not parse VARNAMES/ROW from {path}")
    if len(row) < len(varnames):
        raise RuntimeError(f"ROW has fewer columns than VARNAMES in {path}")
    return dict(zip(varnames, row[: len(varnames)]))


def run_dump(snlc_sim: Path, simlib_path: Path, workdir: Path) -> Path:
    env = os.environ.copy()
    env.setdefault("SNANA_DIR", str(snlc_sim.parent.parent))
    cmd = [
        str(snlc_sim),
        "NOFILE",
        "SIMLIB_FILE",
        str(simlib_path),
        "SIMLIB_DUMP",
        "2",
        "GENFILTERS",
        "r",
        "SIMLIB_MINOBS",
        "0",
    ]
    subprocess.run(
        cmd,
        cwd=str(workdir),
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        env=env,
    )
    dumps = sorted(workdir.glob("SIMLIB_DUMP_OBS_*.TEXT"))
    if len(dumps) != 1:
        raise RuntimeError(f"Expected 1 OBS dump file, found {len(dumps)}: {dumps}")
    return dumps[0]


def make_simlib(psf_model: str, beta: float) -> str:
    lines = [
        "DOCUMENTATION:",
        "  PURPOSE: test_moffat_psf_mode",
        "DOCUMENTATION_END:",
        "",
        "SURVEY:   LSST",
        "FILTERS:  r",
        "PIXSIZE:  0.200",
        "PSF_UNIT: ARCSEC_FWHM",
    ]
    if psf_model.upper() == "MOFFAT":
        lines += [
            "PSF_MODEL:   MOFFAT",
            f"MOFFAT_BETA: {beta:.6f}",
        ]
    lines += [
        "",
        "BEGIN LIBGEN",
        "",
        "LIBID: 1",
        "NOBS:  1  RA:  0.0  DEC:  0.0",
        "",
        "#     MJD        ID   FLT GAIN NOISE SKYSIG  PSF1  PSF2  RATIO  ZPTAVG ZPTERR  MAG",
        "S: 60000.0 1 r 1.60 3.75 10.0 0.800 0.0 0.0 31.500 0.005 99.0",
        "END_LIBID: 1",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    snana_dir = Path(__file__).resolve().parents[1]
    snlc_sim = snana_dir / "bin" / "snlc_sim.exe"
    if not snlc_sim.exists():
        raise SystemExit(f"Missing {snlc_sim}; build SNANA first.")

    os.environ.setdefault("SNANA_DIR", str(snana_dir))

    sndata_root = os.environ.get("SNDATA_ROOT")
    if not sndata_root or not Path(sndata_root).exists():
        raise SystemExit("Missing or invalid $SNDATA_ROOT; source your SNANA setup first.")

    beta = 2.032994  # chosen so NEA_M/NEA_G ~ 2.451 at fixed FWHM
    expected_ratio = nea_ratio_moffat_over_gauss(beta)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        simlib_gauss = tmp_path / "TEST_GAUSS.SIMLIB"
        simlib_moffat = tmp_path / "TEST_MOFFAT.SIMLIB"
        simlib_gauss.write_text(make_simlib("GAUSS", beta), encoding="utf-8")
        simlib_moffat.write_text(make_simlib("MOFFAT", beta), encoding="utf-8")

        dump_gauss = run_dump(snlc_sim, simlib_gauss, tmp_path)
        row_gauss = parse_dump_obs(dump_gauss)
        nea_gauss = float(row_gauss["NEA"])
        fpix_gauss = float(row_gauss["FPIX"])

        # Remove the previous dump to avoid ambiguity.
        dump_gauss.unlink(missing_ok=True)

        dump_moffat = run_dump(snlc_sim, simlib_moffat, tmp_path)
        row_moffat = parse_dump_obs(dump_moffat)
        nea_moffat = float(row_moffat["NEA"])
        fpix_moffat = float(row_moffat["FPIX"])

    ratio = nea_moffat / nea_gauss
    if not (abs(ratio - expected_ratio) / expected_ratio < 0.02):
        raise SystemExit(
            f"FAIL: NEA ratio {ratio:.4f} differs from expected {expected_ratio:.4f} (beta={beta})"
        )

    if not (fpix_moffat < fpix_gauss):
        raise SystemExit(
            f"FAIL: expected FPIX(MOFFAT) < FPIX(GAUSS) but got {fpix_moffat:.6f} vs {fpix_gauss:.6f}"
        )

    print(
        "OK:",
        f"NEA ratio={ratio:.4f} (expected {expected_ratio:.4f}),",
        f"FPIX gauss={fpix_gauss:.6f} moffat={fpix_moffat:.6f}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

