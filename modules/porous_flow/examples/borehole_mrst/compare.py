#!/usr/bin/env python3
"""
Compare the finite-element Peaceman boreholes with the finite-volume results that were verified
against MRST.

gold/ holds the gold CSVs of FVPorousFlowPeacemanBorehole (idaholab/moose#34110), which
agree with MRST 2023b ad-blackoil to 0.4 Pa in bottom-hole pressure and 4e-6 kg/s in rate for the
single-phase cases.  Each finite-element input here reproduces one of those cases on the same
grid, with the porepressure at the element corners rather than the cell centres.

Every case is run with three values of re_constant: Peaceman's finite-difference value 0.28;
0.1594, the value Chen and Zhang derive for rectangular finite elements with the well at a node;
and 0.6, which is not derived but fitted to the single-phase rate case on this grid.  The
completions here are at element centres, where the borehole sees the average of the eight corner
pressures, so the equivalent radius is larger than Peaceman's rather than smaller.  For a
rate-controlled case the bottom-hole pressure is compared; for a pressure-controlled case, the
rates.

Usage, from this directory:
    python3 compare.py [path/to/porous_flow-opt]
"""

import csv
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# re_constant and its label in the table and plot
RE_CONSTANTS = {"0.28": "0.28", "0.1594": "0.1594", "0.6": "0.6 (fitted to this grid)"}

# input file, reference CSV, compared columns, label
CASES = [
    ("1ph_rate", "fv_borehole_rate_out", ["p_bot"], "1 phase, rate control"),
    ("1ph_bhp", "fv_borehole_bhp_out", ["bh_rate"], "1 phase, bhp control"),
    ("2ph_rate", "fv_borehole_2ph_rate_out", ["p_bot"], "gas injection, rate control"),
    ("2ph_bhp", "fv_borehole_2ph_bhp_out", ["gas_rate"], "gas injection, bhp control"),
    (
        "2ph_bhp_produce",
        "fv_borehole_2ph_bhp_produce_out",
        ["water_rate", "gas_rate"],
        "water and gas production, bhp control",
    ),
]


def read(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    # drop the t = 0 row, where every postprocessor is zero
    return [r for r in rows if float(r["time"]) > 0.0]


def run(exe, case, re_constant):
    base = f"{case}_re{re_constant}"
    subprocess.run(
        [
            exe,
            "-i",
            f"{case}.i",
            "--color",
            "off",
            f"GlobalParams/re_constant={re_constant}",
            f"Outputs/file_base={base}",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    return read(f"{base}.csv")


def main():
    exe = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "../../porous_flow-opt")

    results = {}
    for case, ref, _, _ in CASES:
        results[case] = {"fv": read(os.path.join("gold", f"{ref}.csv"))}
        for re_constant in RE_CONSTANTS:
            results[case][re_constant] = run(exe, case, re_constant)

    print(f"{'case':42s} {'quantity':11s} {'FV range':>22s}", end="")
    for re_label in RE_CONSTANTS.values():
        print(f" {'max |FE-FV| re=' + re_label:>40s}", end="")
    print()
    for case, _, columns, label in CASES:
        fv = results[case]["fv"]
        for col in columns:
            ref = [float(r[col]) for r in fv]
            print(f"{label:42s} {col:11s} {min(ref):10.5g} .. {max(ref):<10.5g}", end="")
            for re_constant in RE_CONSTANTS:
                fe = [float(r[col]) for r in results[case][re_constant]]
                diff = max(abs(a - b) for a, b in zip(fe, ref))
                print(f" {diff:40.4g}", end="")
            print()

    n_panels = sum(len(c[2]) for c in CASES)
    fig, axes = plt.subplots(n_panels, 1, figsize=(7, 2.6 * n_panels), sharex=True)
    panel = 0
    for case, _, columns, label in CASES:
        for col in columns:
            ax = axes[panel]
            panel += 1
            data = results[case]
            t = [float(r["time"]) for r in data["fv"]]
            ax.plot(t, [float(r[col]) for r in data["fv"]], "k-", label="FV (matches MRST)")
            for (re_constant, re_label), style in zip(RE_CONSTANTS.items(), ["C0o", "C1s", "C2^"]):
                ax.plot(
                    t,
                    [float(r[col]) for r in data[re_constant]],
                    style,
                    markersize=4,
                    label=f"FE re_constant={re_label}",
                )
            ax.set_title(f"{label}: {col}", fontsize=9)
            ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[-1].set_xlabel("time (s)")
    fig.tight_layout()
    fig.savefig("borehole_mrst_comparison.png", dpi=120)


if __name__ == "__main__":
    main()
