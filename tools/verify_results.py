#!/usr/bin/env python3
"""Check that the published numbers follow from the published data.

Run from the repository root:

    python3 tools/verify_results.py           # check only
    python3 tools/verify_results.py --write   # also write data/summary.json

Three checks, each printed as PASS or FAIL:

  A. The shipped step-5 scripts (analysis/step5_calc_Rk) are re-run on a copy of
     data/per_pressure and data/heat_flux_J.dat in a temporary folder. Their
     outputs must match the tables shipped in data/.
  B. R_K is recomputed with independent code (this file), from the same inputs
     and the same definition, and compared with data/Rk_pfh.dat. The heat flux
     is refitted from data/heat_flux_qInOut and compared with
     data/heat_flux_J.dat.
  C. Every number in the README results tables is compared with data/.
  D. The shipped spectral-overlap script (analysis/step6_sdhf) is re-run on
     data/spectral_heat_flux; the effect of the smoothing on the integrated
     spectrum, and the agreement of that integral with the thermostat heat
     flow, are checked.
  E. The numbers quoted in the README's running text (trend, ranges,
     temperatures, cross-section) are recomputed and looked up as written.

Needs Python 3.9+ with numpy, scipy and matplotlib. Exit code 0 means every
check passed.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from scipy.stats import linregress

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
STEP5 = os.path.join(ROOT, "analysis", "step5_calc_Rk")
PRESSURES = [50, 100, 200, 250, 300, 350, 400]

# Choices fixed in the shipped analysis (see docs/analysis-choices.md).
N_SHEETS = 9          # graphene sheets per wall
N_BLOCKS = 5          # time blocks for the error estimate
FIT_REGION = 100      # central PFH temperature bins used in the linear fit
FIT_FROM = 500        # thermostat-energy samples skipped as transient (x 10 ps)
SAMPLE_PS = 10.0      # thermostat energy written every 1e4 steps of 1 fs

results = []


def report(name, ok, detail=""):
    results.append((name, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))


def load_flux():
    """Steady heat flux per pressure (W/m^2) from data/heat_flux_J.dat."""
    d = np.loadtxt(os.path.join(DATA, "heat_flux_J.dat"))
    return {int(r[0]): r for r in d}


# ---------------------------------------------------------------- check A
def check_shipped_scripts():
    flux = load_flux()
    tmp = tempfile.mkdtemp(prefix="pfh_verify_")
    try:
        outs = os.path.join(tmp, "outs")
        os.makedirs(outs)
        os.makedirs(os.path.join(tmp, "plots"))
        for f in os.listdir(os.path.join(DATA, "per_pressure")):
            if f.endswith(("_cool_sheets.dat", "_hot_sheets.dat", "_pfh_temp.dat", "_pfh_dens.dat")):
                shutil.copy(os.path.join(DATA, "per_pressure", f), outs)
        for s in ["comp.Rk.block.py", "comp.Rk.py", "cmp.Rk_graphene_pfh.py", "comp.peak_bulk_dens.py"]:
            shutil.copy(os.path.join(STEP5, s), tmp)
        env = dict(os.environ, MPLBACKEND="Agg")

        def run(args):
            r = subprocess.run([sys.executable, "-I"] + args, cwd=tmp, env=env,
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"{' '.join(args)} failed:\n{r.stderr[-2000:]}")

        for p in PRESSURES:
            j = f"{flux[p][3]:.5e}"
            run(["comp.Rk.block.py", j, str(p)])
            run(["comp.Rk.py", j, str(p)])
        run(["cmp.Rk_graphene_pfh.py"])
        run(["comp.peak_bulk_dens.py"])

        for kind in ["kapitza_stats", "kapitza_average"]:
            same = []
            for p in PRESSURES:
                a = open(os.path.join(outs, f"{p}bar_{kind}.dat")).read()
                b = open(os.path.join(DATA, "per_pressure", f"{p}bar_{kind}.dat")).read()
                same.append(a == b)
            report(f"A. shipped scripts reproduce data/per_pressure/*_{kind}.dat",
                   all(same), f"{sum(same)}/{len(same)} files identical")

        a = np.loadtxt(os.path.join(tmp, "Rk_pfh.dat"))
        b = np.loadtxt(os.path.join(DATA, "Rk_pfh.dat"))
        report("A. shipped scripts reproduce data/Rk_pfh.dat", np.allclose(a, b, rtol=1e-12, atol=0))

        a = np.loadtxt(os.path.join(tmp, "density_stats.dat"))
        b = np.loadtxt(os.path.join(DATA, "density_stats.dat"))
        report("A. shipped scripts reproduce data/density_stats.dat", np.allclose(a, b, rtol=0, atol=1e-9))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------- check B
def rk_blocks(p, j, pfh_t_scale=1.0):
    """Independent R_K per block for one pressure.

    Definition (same as the shipped analysis): the PFH temperature profile,
    averaged over the two mirror-image slabs, is fitted with a straight line
    over its central FIT_REGION bins and extrapolated to the position of the
    graphene sheet in contact with the PFH (the innermost sheet of each wall).
    The jump between that sheet's temperature and the extrapolated PFH
    temperature, divided by the heat flux j, is R_K. The hot and cold
    interfaces are averaged. pfh_t_scale multiplies the PFH temperatures (used to
    test the effect of the carbon mass in the PFH temperature code).
    """
    pp = os.path.join(DATA, "per_pressure")
    hot = np.loadtxt(os.path.join(pp, f"{p}bar_hot_sheets.dat"))
    cool = np.loadtxt(os.path.join(pp, f"{p}bar_cool_sheets.dat"))
    temp = np.loadtxt(os.path.join(pp, f"{p}bar_pfh_temp.dat"))
    nb = temp.shape[0] // 2
    lo, hi = (nb - FIT_REGION) // 2, (nb + FIT_REGION) // 2
    z = temp[:nb, 0]
    z_hot = hot[N_SHEETS - 1, 1]     # innermost hot-wall sheet (right domain)
    z_cool = cool[N_SHEETS - 1, 1]   # innermost cold-wall sheet (right domain)
    out = {"hot": [], "cool": [], "av": [], "k": []}
    for b in range(N_BLOCKS):
        t_pfh = pfh_t_scale * 0.5 * (temp[:nb, 1 + b] + temp[nb:, 1 + b])
        t_hot = 0.5 * (hot[N_SHEETS - 1, 3 + b] + hot[2 * N_SHEETS - 1, 3 + b])
        t_cool = 0.5 * (cool[N_SHEETS - 1, 3 + b] + cool[2 * N_SHEETS - 1, 3 + b])
        slope, icpt = np.polyfit(z[lo:hi], t_pfh[lo:hi], 1)
        r_hot = (t_hot - (slope * z_hot + icpt)) / j
        r_cool = ((slope * z_cool + icpt) - t_cool) / j
        out["hot"].append(r_hot)
        out["cool"].append(r_cool)
        out["av"].append(0.5 * (r_hot + r_cool))
        out["k"].append(j / (1e10 * abs(slope)))
    return {k: np.array(v) for k, v in out.items()}


def check_independent():
    flux = load_flux()
    shipped = {int(r[0]): r for r in np.loadtxt(os.path.join(DATA, "Rk_pfh.dat"))}
    summary = []
    ok_rk = True
    for p in PRESSURES:
        j = flux[p][3]
        r = rk_blocks(p, j)
        mean, sd = r["av"].mean(), r["av"].std()   # population sd (ddof=0), as shipped
        # The shipped summary is built from block values printed to 5 significant figures,
        # so 1e-4 (mean) and 1e-3 (sd) relative, with no absolute slack, is agreement to rounding.
        ok_rk &= bool(np.isclose(mean, shipped[p][1], rtol=1e-4, atol=0)
                      and np.isclose(sd, shipped[p][2], rtol=1e-3, atol=0))
        summary.append({
            "pressure_bar": p,
            "heat_flux_W_m2": float(j),
            "R_K_m2K_W": float(mean),
            "R_K_sd_blocks_m2K_W": float(sd),
            "R_K_sem_m2K_W": float(r["av"].std(ddof=1) / np.sqrt(N_BLOCKS)),
            "R_K_hot_interface_m2K_W": float(r["hot"].mean()),
            "R_K_cold_interface_m2K_W": float(r["cool"].mean()),
            "PFH_thermal_conductivity_W_mK": float(r["k"].mean()),
        })
    report("B. independent R_K matches data/Rk_pfh.dat (mean and sd, all pressures)", ok_rk)

    # Heat flux: refit the thermostat energy and back out the cross-section it implies.
    areas = []
    ok_fit = True
    for p in PRESSURES:
        d = np.loadtxt(os.path.join(DATA, "heat_flux_qInOut", f"out.{p}_bar"))
        t = np.arange(d.shape[0]) * SAMPLE_PS
        e = 0.5 * (d[:, 2] - d[:, 1])          # mean of energy removed and energy added (eV)
        fit = linregress(t[FIT_FROM:], e[FIT_FROM:])
        q_per_interface_eV_ps = 0.5 * fit.slope  # heat splits over two paths from the hot wall
        ok_fit &= bool(np.isclose(q_per_interface_eV_ps, flux[p][2], rtol=1e-4))
        area_A2 = q_per_interface_eV_ps * 1.60218e-19 / 1e-12 / flux[p][3] / 1e-20
        areas.append(area_A2)
    report("B. refitted heat flow matches data/heat_flux_J.dat (eV/ps column)", ok_fit)
    # 100 bar started from its own configuration; the other six share one.
    a100 = areas[PRESSURES.index(100)]
    others = [a for p, a in zip(PRESSURES, areas) if p != 100]
    ok_area = (max(others) - min(others)) / np.mean(others) < 1e-4 and abs(a100 - np.mean(others)) / np.mean(others) < 2e-3
    report("B. cross-section implied by J: one value at six pressures, within 0.2 % at 100 bar",
           ok_area, f"{np.mean(others):.1f} A^2; 100 bar {a100:.1f} A^2")

    dens = {int(r[0]): r for r in np.loadtxt(os.path.join(DATA, "density_stats.dat"))}
    for s in summary:
        r = dens[s["pressure_bar"]]
        s.update({"bulk_density_kg_m3": float(r[1]),
                  "first_peak_density_hot_interface_kg_m3": float(r[2]),
                  "first_peak_distance_from_sheet_A": float(r[3])})
    return summary, areas


# ---------------------------------------------------------------- check D
def daniell(y, w):
    """Moving average over w points, made odd (same rule as analysis/step6_sdhf)."""
    if w % 2 == 0:
        w += 1
    return np.convolve(y, np.ones(w) / w, mode="same")


def check_spectral():
    """Spectral heat flux: re-run the shipped overlap script; check the smoothing
    claim and that the integrated spectrum equals the thermostat heat flow."""
    from scipy import integrate
    flux = load_flux()
    tmp = tempfile.mkdtemp(prefix="pfh_verify_sdhf_")
    try:
        run_dir = os.path.join(tmp, "step6_sdhf")
        os.makedirs(os.path.join(run_dir, "sdhf_consolidated"))
        for f in os.listdir(os.path.join(DATA, "spectral_heat_flux")):
            if f.endswith(("_fluid.dat", "_solid.dat")):
                shutil.copy(os.path.join(DATA, "spectral_heat_flux", f), os.path.join(run_dir, "sdhf_consolidated"))
        shutil.copy(os.path.join(ROOT, "analysis", "step6_sdhf", "plot.sdhf_overlap.py"), run_dir)
        r = subprocess.run([sys.executable, "-I", "plot.sdhf_overlap.py"], cwd=run_dir,
                           env=dict(os.environ, MPLBACKEND="Agg"), capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(r.stderr[-2000:])
        a = np.loadtxt(os.path.join(run_dir, "sdhf_overlaps.dat"))
        b = np.loadtxt(os.path.join(DATA, "spectral_heat_flux", "sdhf_overlaps.dat"))
        report("D. shipped script reproduces data/spectral_heat_flux/sdhf_overlaps.dat",
               np.allclose(a, b, rtol=1e-4, atol=0))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    worst_smooth, worst_flux = 0.0, 0.0
    for p in PRESSURES:
        for side in ["fluid", "solid"]:
            d = np.loadtxt(os.path.join(DATA, "spectral_heat_flux", f"{p}bar_{side}.dat"), skiprows=1)
            f, ptot = d[:, 0], d[:, 4]
            w = int(round(0.5 / ((f[1] - f[0]) / 1e12)))
            raw = integrate.trapezoid(ptot, f)
            smooth = integrate.trapezoid(daniell(ptot, w), f)
            worst_smooth = max(worst_smooth, abs(smooth - raw) / abs(raw))
            # The shipped spectra integrate over frequency (Hz) to eV/ps, the unit of flux[p][2].
            worst_flux = max(worst_flux, abs(abs(raw) - flux[p][2]) / flux[p][2])
    report("D. smoothing changes the integrated spectral flux by less than 0.7 %",
           worst_smooth < 0.007, f"largest change {100 * worst_smooth:.2f} %")
    report("D. integrated spectral flux equals the thermostat heat flow within 1.6 % (README claim)",
           worst_flux < 0.016, f"largest difference {100 * worst_flux:.2f} %")
    return worst_smooth, worst_flux


# ---------------------------------------------------------------- check C
def readme_rows(header_start):
    """Rows of the first README table whose header line starts with header_start."""
    lines = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read().splitlines()
    rows, inside = [], False
    for ln in lines:
        if ln.startswith("|") and ln.lstrip("| ").startswith(header_start):
            inside = True
            continue
        if inside:
            if not ln.startswith("|"):
                break
            if set(ln.replace("|", "").strip()) <= set("-: "):
                continue
            rows.append([c.strip() for c in ln.strip("|").split("|")])
    return rows


def num(cell):
    return [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", cell)]


def check_readme(summary):
    by_p = {s["pressure_bar"]: s for s in summary}
    rows = readme_rows("Pressure (bar) | R_K")
    ok = sorted(int(num(r[0])[0]) for r in rows) == PRESSURES and all(len(r) == 4 for r in rows)
    for r in rows if ok else []:
        p = int(num(r[0])[0])
        s = by_p[p]
        m, sd = num(r[1])[:2]
        hot, cold = num(r[2])[0], num(r[3])[0]
        ok &= abs(m - round(s["R_K_m2K_W"] * 1e8, 2)) < 1e-9 and abs(sd - round(s["R_K_sd_blocks_m2K_W"] * 1e8, 2)) < 1e-9
        ok &= abs(hot - round(s["R_K_hot_interface_m2K_W"] * 1e8, 2)) < 1e-9
        ok &= abs(cold - round(s["R_K_cold_interface_m2K_W"] * 1e8, 2)) < 1e-9
    report("C. README R_K table matches data (each pressure once, all four columns)", ok, f"{len(rows)} rows")

    rows = readme_rows("Pressure (bar) | Bulk density")
    ok = sorted(int(num(r[0])[0]) for r in rows) == PRESSURES and all(len(r) == 4 for r in rows)
    for r in rows if ok else []:
        p = int(num(r[0])[0])
        s = by_p[p]
        ok &= round(s["bulk_density_kg_m3"]) == num(r[1])[0]
        ok &= round(s["first_peak_density_hot_interface_kg_m3"]) == num(r[2])[0]
        ok &= abs(round(s["first_peak_distance_from_sheet_A"], 2) - num(r[3])[0]) < 1e-9
    report("C. README density table matches data (each pressure once)", ok, f"{len(rows)} rows")


def check_prose(summary, areas, worst_flux):
    """Numbers quoted in the README text, recomputed and looked up as written."""
    text = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()
    s = {d["pressure_bar"]: d for d in summary}
    lo, hi = s[PRESSURES[0]], s[PRESSURES[-1]]
    rk = np.array([s[p]["R_K_m2K_W"] for p in PRESSURES])
    sd = np.array([s[p]["R_K_sd_blocks_m2K_W"] for p in PRESSURES])
    x = np.array(PRESSURES, float)

    def slope_over_sigma(sig):
        # weighted least squares, errors taken as absolute
        a = np.vstack([x, np.ones_like(x)]).T
        w = 1 / sig ** 2
        cov = np.linalg.inv(a.T @ (a * w[:, None]))
        beta = cov @ (a.T @ (w * rk))
        return beta[0], abs(beta[0]) / np.sqrt(cov[0, 0])

    slope, z_sd = slope_over_sigma(sd)
    _, z_sem = slope_over_sigma(sd / 2)     # standard error = population sd / 2 for 5 blocks
    drop = [100 * (1 - s[p]["R_K_hot_interface_m2K_W"] / s[p]["R_K_cold_interface_m2K_W"]) for p in PRESSURES]
    t_hot, t_cold = [], []
    for p in PRESSURES:
        h = np.loadtxt(os.path.join(DATA, "per_pressure", f"{p}bar_hot_sheets.dat"))
        c = np.loadtxt(os.path.join(DATA, "per_pressure", f"{p}bar_cool_sheets.dat"))
        t_hot.append(np.mean(0.5 * (h[N_SHEETS - 1, 3:] + h[2 * N_SHEETS - 1, 3:])))
        t_cold.append(np.mean(0.5 * (c[N_SHEETS - 1, 3:] + c[2 * N_SHEETS - 1, 3:])))
    j = [s[p]["heat_flux_W_m2"] / 1e8 for p in PRESSURES]
    flux = load_flux()

    # Hot below cold, block by block.
    n_lower = sum(int(np.sum(rk_blocks(p, flux[p][3])["hot"] < rk_blocks(p, flux[p][3])["cool"])) for p in PRESSURES)

    # Known issue 1: time axis from recorded timesteps instead of row numbers.
    gain = []
    for p in PRESSURES:
        d = np.loadtxt(os.path.join(DATA, "heat_flux_qInOut", f"out.{p}_bar"))
        e = 0.5 * (d[:, 2] - d[:, 1])
        by_row = linregress(np.arange(d.shape[0])[FIT_FROM:] * SAMPLE_PS, e[FIT_FROM:]).slope
        by_step = linregress(d[FIT_FROM:, 0] * 1e-3, e[FIT_FROM:]).slope     # 1 fs steps -> ps
        gain.append(100 * (by_step / by_row - 1))

    # Known issue 2: carbon mass 12.017 in the PFH temperature code (6 of 20 atoms are carbon).
    eps = (6 / 20) * (12.017 / 12.011 - 1)
    d_hot, d_cold, d_av = [], [], []
    for p in PRESSURES:
        a, b = rk_blocks(p, flux[p][3]), rk_blocks(p, flux[p][3], 1 / (1 + eps))
        d_hot.append(100 * (b["hot"].mean() / a["hot"].mean() - 1))
        d_cold.append(100 * (b["cool"].mean() / a["cool"].mean() - 1))
        d_av.append(100 * (b["av"].mean() / a["av"].mean() - 1))

    a100 = areas[PRESSURES.index(100)]
    a_rest = np.mean([a for p, a in zip(PRESSURES, areas) if p != 100])
    expected = {
        "bulk density rise": f"bulk density rises by {100 * (hi['bulk_density_kg_m3'] / lo['bulk_density_kg_m3'] - 1):.0f} %",
        "first peak rise": f"first density peak by {100 * (hi['first_peak_density_hot_interface_kg_m3'] / lo['first_peak_density_hot_interface_kg_m3'] - 1):.0f} %",
        "R_K fall": f"R_K falls by {100 * (1 - hi['R_K_m2K_W'] / lo['R_K_m2K_W']):.1f} %",
        "trend slope": f"{abs(slope) * 100 * 1e8:.2f} × 10⁻⁸ m² K/W per 100 bar",
        "trend / sd": f"{z_sd:.1f} times its uncertainty",
        "trend / sem": f"{z_sem:.1f} times when",
        "hot vs cold": f"{min(drop):.0f} to {max(drop):.0f} % lower",
        "interface temperatures": f"about {np.mean(t_hot):.0f} K (hot) and {np.mean(t_cold):.0f} K (cold)",
        "heat flux range": f"{min(j):.1f} to {max(j):.1f} × 10⁸ W/m²",
        "spectral vs thermostat": f"within {100 * worst_flux:.1f} %",
        "cross-section": f"{a_rest:.0f} Å² ({a100:.0f} Å² at 100 bar)",
        "hot below cold, blocks": f"lower in {n_lower} of the {5 * len(PRESSURES)} time blocks",
        "time-axis issue": f"heat flow by {np.mean(gain):.1f} %",
        "carbon-mass issue, PFH temperature": f"PFH temperatures by {100 * eps:.3f} %",
        "carbon-mass issue, R_K": (f"raises R_K at the hot interface by {np.mean(d_hot):.1f} %, lowers it at the cold "
                                   f"interface by {-np.mean(d_cold):.1f} % and raises the averaged R_K by {np.mean(d_av):.2f} %"),
        "known-issue summary": f"R_K by {np.mean(gain):.1f} % and {np.mean(d_av):.2f} %",
    }
    flat = " ".join(text.split())   # ignore line wrapping
    missing = [f"{k}: '{v}'" for k, v in expected.items() if " ".join(v.split()) not in flat]
    report("E. numbers in the README text match the data", not missing,
           f"{len(expected) - len(missing)}/{len(expected)} found" + ("; missing " + "; ".join(missing) if missing else ""))


def main():
    write = "--write" in sys.argv[1:]
    check_shipped_scripts()
    summary, areas = check_independent()
    worst_smooth, worst_flux = check_spectral()
    check_readme(summary)
    check_prose(summary, areas, worst_flux)
    if write:
        out = {
            "description": "Graphene / n-perfluorohexane NEMD: results per pressure, recomputed by tools/verify_results.py from data/.",
            "units": {"R_K": "m^2 K / W", "heat_flux": "W / m^2", "density": "kg / m^3",
                      "distance": "Angstrom", "thermal_conductivity": "W / (m K)"},
            "error_definitions": {
                "R_K_sd_blocks": "population standard deviation (ddof=0) of 5 time-block values; the error bar used in the README and figures",
                "R_K_sem": "standard error of the mean of the 5 block values (ddof=1)"},
            "cross_section_implied_by_J_A2": [round(a, 2) for a in areas],
            "spectral_checks": {
                "largest_change_from_smoothing_fraction": round(worst_smooth, 5),
                "largest_difference_integrated_spectrum_vs_thermostat_flow_fraction": round(worst_flux, 5)},
            "results": summary,
        }
        with open(os.path.join(DATA, "summary.json"), "w") as f:
            json.dump(out, f, indent=2)
        print("wrote data/summary.json")
    n_fail = sum(1 for _, ok in results if not ok)
    print(f"\n{len(results) - n_fail}/{len(results)} checks passed")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
