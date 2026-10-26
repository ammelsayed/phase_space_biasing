#!/usr/bin/env python3
"""Analyze the phase-space biasing study.

For every (sample, module) run produced by study/run_study.sh, read the
unweighted LHE events and compare, for ptl1/ptl2/met/ht:
  * the physical (weighted) distributions -> closure check: with
    impact_xsec=.False. every biased sample must reproduce the reference;
  * the unweighted population per bin -> where the bias oversamples;
  * the per-bin relative statistical uncertainty sqrt(sum w^2)/sum w,
    and its ratio to the reference -> the statistics actually gained.

Outputs plots into <outdir> and a markdown summary on stdout.

Usage: python3 study/analyze_lhe.py [runs_dir] [outdir]
"""
import gzip
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LEPTONS = {11, 13, 15}
INVIS = {12, 14, 16}

SAMPLES = {
    "z": {
        "proc": "z012j_MLM",
        "tags": ["z012j_ref", "z012j_ptl1", "z012j_xsecflat"],
        "labels": {"z012j_ref": "no bias", "z012j_ptl1": "ptl1_bias",
                   "z012j_xsecflat": "xsec_flat_bias (exp 1)"},
        "obs": ["ptl1", "ptl2", "ht"],
        "bins": {"ptl1": np.arange(150, 1151, 20.0),
                 "ptl2": np.arange(50, 851, 20.0),
                 "met": np.arange(50, 1051, 25.0),
                 "ht": np.arange(200, 1601, 40.0)},
        "tails": {"ptl1": [300, 500, 800], "ptl2": [150, 250, 400],
                  "met": [], "ht": [400, 600, 1000]},
    },
    "w": {
        "proc": "w012j_MLM",
        "tags": ["w012j_ref", "w012j_ptl1", "w012j_met", "w012j_ptl1_met",
                 "w012j_ptl1_met_soft", "w012j_xsecflat"],
        "labels": {"w012j_ref": "no bias", "w012j_ptl1": "ptl1_bias",
                   "w012j_met": "met_bias",
                   "w012j_ptl1_met": "ptl1_met_bias (4+4)",
                   "w012j_ptl1_met_soft": "ptl1_met_bias (2+2)",
                   "w012j_xsecflat": "xsec_flat_bias (exp 1)"},
        "obs": ["ptl1", "met", "ht"],
        "bins": {"ptl1": np.arange(150, 1151, 20.0),
                 "ptl2": np.arange(50, 851, 20.0),
                 "met": np.arange(50, 1051, 25.0),
                 "ht": np.arange(150, 1151, 20.0)},
        "tails": {"ptl1": [300, 500, 800], "ptl2": [],
                  "met": [150, 250, 400], "ht": [300, 500, 800]},
    },
}

COLORS = ["k", "tab:blue", "tab:green", "tab:red", "tab:purple", "tab:orange"]


def read_lhe(path):
    """Yield per event: (pdgs, istups, pxs, pys, weight, bias_weight)."""
    events = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="ignore") as f:
        inevent = False
        npart = -1
        pdgs, istups, pxs, pys = [], [], [], []
        w, bias = 0.0, 1.0
        for line in f:
            s = line.strip()
            if s.startswith("<event"):
                inevent, npart = True, -1
                pdgs, istups, pxs, pys = [], [], [], []
                w, bias = 0.0, 1.0
                continue
            if s.startswith("</event>"):
                events.append((pdgs, istups, pxs, pys, w, bias))
                inevent = False
                continue
            if not inevent:
                continue
            if npart == -1:
                hdr = s.split()
                if len(hdr) >= 3:
                    npart, w = int(hdr[0]), float(hdr[2])
                continue
            if npart > 0:
                cols = s.split()
                pdgs.append(int(cols[0]))
                istups.append(int(cols[1]))
                pxs.append(float(cols[6]))
                pys.append(float(cols[7]))
                npart -= 1
                continue
            if "wgt id='bias'" in s or 'wgt id="bias"' in s:
                try:
                    bias = float(s.split(">")[1].split("<")[0])
                except (ValueError, IndexError):
                    pass
    return events


def observables(events):
    out = {k: [] for k in ("ptl1", "ptl2", "met", "ht")}
    for pdgs, istups, pxs, pys, w, bias in events:
        lpt = sorted(
            (math.hypot(px, py) for g, i, px, py in zip(pdgs, istups, pxs, pys)
             if i == 1 and abs(g) in LEPTONS),
            reverse=True,
        )
        metx = sum(-px for g, i, px in zip(pdgs, istups, pxs)
                   if i == 1 and abs(g) in INVIS)
        mety = sum(-py for g, i, py in zip(pdgs, istups, pys)
                   if i == 1 and abs(g) in INVIS)
        ht = sum(math.hypot(px, py) for g, i, px, py in zip(pdgs, istups, pxs, pys)
                 if i == 1 and abs(g) not in INVIS)
        out["ptl1"].append(lpt[0] if lpt else np.nan)
        out["ptl2"].append(lpt[1] if len(lpt) > 1 else np.nan)
        out["met"].append(math.hypot(metx, mety))
        out["ht"].append(ht)
    return {k: np.array(v) for k, v in out.items()}


def load_run(runs_dir, tag):
    path = None
    for cand in (f"events_{tag}.lhe", f"events_{tag}.lhe.gz"):
        p = os.path.join(runs_dir, cand)
        if os.path.exists(p):
            path = p
            break
    if path is None:
        return None
    events = read_lhe(path)
    obs = observables(events)
    w = np.array([e[4] for e in events])
    bias = np.array([e[5] for e in events])
    return {"obs": obs, "w": w, "bias": bias, "path": path, "n": len(events)}


def whist(x, w, bins):
    """Weighted histogram, its sqrt(sum w^2) error, population, and
    the per-bin relative statistical uncertainty sqrt(sum w^2)/sum w."""
    good = np.isfinite(x)
    x, w = x[good], w[good]
    h, edges = np.histogram(x, bins=bins, weights=w)
    h2, _ = np.histogram(x, bins=bins, weights=w * w)
    pop, _ = np.histogram(x, bins=bins)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(h > 0, np.sqrt(h2) / np.abs(h), np.nan)
    return h, pop, rel, edges


def tail_stats(x, w, thr):
    good = np.isfinite(x)
    x, w = x[good], w[good]
    m = x >= thr
    n = int(m.sum())
    sw2 = float((w[m] ** 2).sum())
    sw = float(w[m].sum())
    rel = math.sqrt(sw2) / sw if sw > 0 else float("nan")
    return n, rel


def load_sigma(runs_dir, tag):
    """Read the reported cross section from runs/xsec_<tag>.txt."""
    p = os.path.join(runs_dir, f"xsec_{tag}.txt")
    if not os.path.exists(p):
        return None, None
    for line in open(p):
        if "Cross-section" in line:
            tok = line.split(":")[1].split()
            try:
                return float(tok[0]), float(tok[2])
            except (IndexError, ValueError):
                return None, None
    return None, None


def main():
    runs_dir = sys.argv[1] if len(sys.argv) > 1 else "runs"
    outdir = sys.argv[2] if len(sys.argv) > 2 else "study/plots"
    os.makedirs(outdir, exist_ok=True)

    summary = {}
    for sname, scfg in SAMPLES.items():
        runs = {}
        for tag in scfg["tags"]:
            run = load_run(runs_dir, tag)
            if run is None:
                print(f"[WARN] missing events for {tag}, skipped")
                continue
            runs[tag] = run
        if not runs:
            continue
        reftag = scfg["tags"][0]
        summary[sname] = {}
        for obs in scfg["obs"]:
            bins = scfg["bins"][obs]
            centers = 0.5 * (bins[:-1] + bins[1:])
            hists = {}
            for i, (tag, run) in enumerate(runs.items()):
                hists[tag] = whist(run["obs"][obs], run["w"], bins)
            ref = hists[reftag]

            fig, (a1, a2, a3) = plt.subplots(
                3, 1, sharex=True, figsize=(7, 8),
                gridspec_kw={"height_ratios": [2, 1.4, 1.4]})
            for i, (tag, run) in enumerate(runs.items()):
                h, pop, rel, _ = hists[tag]
                norm = h.sum() if h.sum() > 0 else 1.0
                lbl = f"{scfg['labels'][tag]}  ({run['n']} evts)"
                a1.step(centers, h / norm, where="mid", color=COLORS[i], label=lbl)
                a2.step(centers, pop, where="mid", color=COLORS[i], label=lbl)
                if tag == reftag:
                    continue
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = rel / ref[2]
                a3.plot(centers, ratio, color=COLORS[i], drawstyle="steps-mid",
                        label=scfg["labels"][tag])
            a1.set_yscale("log")
            a1.set_ylabel(r"$\Sigma w$ normalized")
            a1.legend(fontsize=8)
            a1.set_title(f"{sname.upper()}-sample: {obs}")
            a2.set_yscale("log")
            a2.set_ylabel("events / bin")
            a3.axhline(1.0, color="gray", lw=0.8)
            a3.set_ylabel(r"$\delta_{\rm biased}/\delta_{\rm ref}$")
            a3.set_xlabel(f"{obs} [GeV]")
            a3.set_ylim(0, 5.0)
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, f"{sname}_{obs}.png"), dpi=130)
            plt.close(fig)

            # tail-yield table
            summary[sname][obs] = {}
            for tag, run in runs.items():
                row = {}
                for thr in scfg["tails"][obs]:
                    n, rel = tail_stats(run["obs"][obs], run["w"], thr)
                    row[str(thr)] = {"n": n, "rel_err": rel}
                summary[sname][obs][tag] = row

        # global info
        for tag, run in runs.items():
            sigma, sigma_err = load_sigma(runs_dir, tag)
            summary[sname].setdefault("_info", {})[tag] = {
                "n_events": run["n"],
                "sigma_pb": sigma,
                "sigma_err_pb": sigma_err,
                "sum_w": float(run["w"].sum()),
            }

    with open(os.path.join(outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)

    # markdown report
    print("# Biasing study summary\n")
    for sname, obs_map in summary.items():
        info = obs_map.get("_info", {})
        print(f"## {sname.upper()} sample\n")
        for tag, d in info.items():
            sig = (f"{d['sigma_pb']:.4g} +- {d['sigma_err_pb']:.2g} pb"
                   if d["sigma_pb"] is not None else "n/a")
            print(f"- `{tag}`: {d['n_events']} events, sigma = {sig}, "
                  f"sum(w) = {d['sum_w']:.4g}")
        print()
        for obs, tags in obs_map.items():
            if obs.startswith("_"):
                continue
            reftag = [t for t in tags if "ref" in t][0]
            print(f"### {obs}\n")
            print("| threshold [GeV] | " + " | ".join(tags) + " |")
            print("|---" * (len(tags) + 1) + "|")
            for thr in next(iter(tags.values())):
                ref_rel = tags[reftag][str(thr)]["rel_err"]
                cells = []
                for tag in tags:
                    d = tags[tag][str(thr)]
                    if tag == reftag or ref_rel == 0:
                        cells.append(f"{d['n']} / {d['rel_err']:.3f} / -")
                    else:
                        imp = ref_rel / d["rel_err"] if d["rel_err"] else float("nan")
                        cells.append(f"{d['n']} / {d['rel_err']:.3f} / {imp:.1f}x")
                print(f"| {thr} | " + " | ".join(cells) + " |")
            print()


if __name__ == "__main__":
    main()
