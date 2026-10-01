"""Prospective adaptive paired multiprecision timing experiments.

Design:
- 15 paired repetitions for every group.
- Extension to 31 repetitions only when the prospectively frozen ambiguity
  rule is triggered.
- The ambiguity rule is evaluated mechanically from the first 15 pairs:
  extend if the paired percentile bootstrap interval for the ratio of medians
  contains 1 OR the two-IQR separation rule is not met. A direction-consistency
  guard is also required for early stopping.

No pilot timing is mixed into the repeated samples. No method or tolerance is
selected after seeing which method wins.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import platform
import random
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import mpmath
import numpy as np
from mpmath import mp
from mpmath.libmp import BACKEND

from .core import (
    DenseProblem,
    FIELDS,
    coefficients,
    iterate,
    method,
    norm,
    pack_vec,
    text,
    traces,
    unpack_vec,
)
from .costs import eta, work

ROOT = Path(__file__).resolve().parents[1]


def stamp():
    return datetime.now(timezone.utc).isoformat()


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def fhash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def code_hash():
    paths = [ROOT / "run_external_adaptive.py"] + sorted((ROOT / "pn_ext_adaptive").glob("*.py"))
    return digest({str(p.relative_to(ROOT)): fhash(p) for p in paths})


def js(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def csvout(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


def gzout(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False)


def env():
    try:
        import gmpy2
        g = {"version": gmpy2.version(), "gmp": gmpy2.mp_version()}
    except ImportError:
        g = None
    return {
        "utc": stamp(),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "machine": platform.machine(),
        "python": sys.version,
        "executable": sys.executable,
        "numpy": np.__version__,
        "mpmath": mpmath.__version__,
        "mpmath_backend": BACKEND,
        "gmpy2": g,
        "reported_logical_cpus": os.cpu_count(),
        "command": sys.argv,
        "timer": vars(time.get_clock_info("perf_counter")),
        "code_sha256": code_hash(),
        "thread_environment": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "arithmetic": "Serial mpmath mpf / NumPy object arrays; no BLAS in the MP solver",
        "timer_scope": "Full solve including copies, counters and terminal residuals; excluding input setup, error norms, trace serialization and output I/O",
    }


def gid(g):
    return f"{g['suite']}_{g['field']}_n{g['n']}_d{g['delta'].replace('.', 'p')}_e{g['digits']}"


def experiment_plan():
    groups = []
    tolerances = [50, 200, 800, 1024]
    for field in ("H1", "H5"):
        n = FIELDS[field]
        for suite, methods in (("external6", ["P5", "M6"]), ("external8", ["P7", "M8"])):
            for digits in tolerances:
                groups.append({
                    "suite": suite,
                    "field": field,
                    "n": n,
                    "delta": "0.30",
                    "digits": digits,
                    "methods": methods,
                })
    return {
        "schema": 2,
        "status": "PROSPECTIVE",
        "working_dps": 1300,
        "verification_dps": 2600,
        "max_cycles": 20,
        "warmups": 2,
        "initial_repetitions": 15,
        "extension_target_repetitions": 31,
        "bootstrap_resamples": 10000,
        "seed": 20260927,
        "initialization": "D",
        "reference_kappa": 1,
        "cost_scope": "Canonical reference work only; NOT calibrated multiprecision time weights",
        "input_rule": "Frozen float64 A/B/C lifted as exact dyadics; analytic constants/alpha in MP; one working-precision x0 exactly lifted for validation",
        "gate_relative_limit": "1e-10",
        "gate_guard_digits": 40,
        "adaptive_rule": {
            "description": "Stop at 15 only if BOTH diagnostics are clear and point in the same direction; otherwise extend mechanically to 31 paired repetitions.",
            "extend_if_bootstrap_interval_contains_one": True,
            "extend_if_two_iqr_filter_not_met": True,
            "require_direction_consistency_for_early_stop": True,
            "decision_is_mechanical_after_initial_stage": True,
        },
        "groups": groups,
        "design_scope": "H1 and H5 only; P5/M6 and P7/M8; delta=0.30; all four pre-proposed tolerances. No group is added or removed based on the observed winner.",
    }


def validate_plan(plan):
    groups = plan["groups"]
    if not groups:
        raise ValueError("Empty experiment plan")
    if plan["working_dps"] < max(g["digits"] for g in groups) + 150:
        raise ValueError("At least 150 guard digits required")
    if plan["verification_dps"] < 2 * plan["working_dps"]:
        raise ValueError("Verification dps must be at least twice working dps")
    r0 = int(plan["initial_repetitions"])
    r1 = int(plan["extension_target_repetitions"])
    if plan["warmups"] < 1 or r0 < 7 or r1 <= r0:
        raise ValueError("Require >=1 warmup, >=7 initial repetitions, and a larger extension target")
    if len({gid(g) for g in groups}) != len(groups):
        raise ValueError("Duplicate group identifiers")
    allowed = {"external6": ["P5", "M6"], "external8": ["P7", "M8"]}
    for g in groups:
        if g["field"] not in ("H1", "H5"):
            raise ValueError("This frozen tranche is restricted to H1/H5")
        if g["n"] != FIELDS[g["field"]]:
            raise ValueError("Field/dimension mismatch")
        if g["suite"] not in allowed or g["methods"] != allowed[g["suite"]]:
            raise ValueError("External-comparator pair mismatch")
        if g["delta"] != "0.30" or g["digits"] not in (50, 200, 800, 1024):
            raise ValueError("Unexpected distance or tolerance")
        for key in g["methods"]:
            method(key)
    rule = plan.get("adaptive_rule", {})
    if not all(rule.get(k) is True for k in (
        "extend_if_bootstrap_interval_contains_one",
        "extend_if_two_iqr_filter_not_met",
        "require_direction_consistency_for_early_stop",
        "decision_is_mechanical_after_initial_stage",
    )):
        raise ValueError("Adaptive rule must match the prospectively declared design")


def input_name(n, d):
    return f"x0_n{n}_d{d.replace('.', 'p')}.json"


def prepare(out, plan):
    directory = Path(out) / "inputs"
    directory.mkdir(parents=True, exist_ok=True)
    for n in sorted({g["n"] for g in plan["groups"]}):
        path = directory / f"coefficients_n{n}.npz"
        if not path.exists():
            A, B, C = coefficients(n)
            np.savez_compressed(path, A=A, B=B, C=C)
    identities = {}
    with mp.workdps(plan["working_dps"]):
        for n, d in sorted({(g["n"], g["delta"]) for g in plan["groups"]}):
            p = DenseProblem("H1", n, directory)
            identities[str(n)] = p.input_sha256
            path = directory / input_name(n, d)
            value = {"dps": mp.dps, "n": n, "delta": d, "x": pack_vec(p.initial(d))}
            if path.exists() and json.loads(path.read_text()) != value:
                raise ValueError("An incompatible x0 already exists; use a new output directory")
            js(path, value)
    js(directory / "identity.json", identities)
    return input_hash(out)


def input_hash(out):
    return digest({p.name: fhash(p) for p in sorted((Path(out) / "inputs").glob("*")) if p.is_file()})


def x0_from(out, g):
    return unpack_vec(json.loads((Path(out) / "inputs" / input_name(g["n"], g["delta"])).read_text())["x"])


def serial(p, r):
    return {
        **{k: r[k] for k in ("method", "status", "cycles", "active_outer_cycles", "active_predictors", "counts")},
        "residual": text(r["residual"]),
        "error": text(norm(r["x"] - p.alpha)),
        "history": traces(p, r),
    }


def precision_gate(out, g, key, plan):
    """Empirical precision/stability check; not interval certification."""
    D, V = plan["working_dps"], plan["verification_dps"]
    timing = {}
    with mp.workdps(D):
        p = DenseProblem(g["field"], g["n"], Path(out) / "inputs")
        x0 = x0_from(out, g)
        start = time.perf_counter()
        low = iterate(p, x0, method(key), mp.mpf(f"1e-{g['digits']}"), plan["max_cycles"], trace=True)
        timing["low_run_s"] = time.perf_counter() - start
        low_doc = serial(p, low)
    with mp.workdps(V):
        p = DenseProblem(g["field"], g["n"], Path(out) / "inputs")
        x0 = x0_from(out, g)
        tol = mp.mpf(f"1e-{g['digits']}")
        start = time.perf_counter()
        high = iterate(p, x0, method(key), tol, plan["max_cycles"], trace=True)
        timing["high_run_s"] = time.perf_counter() - start
        high_doc = serial(p, high)
        rr = norm(p.F(low["x"]))
        ee = norm(low["x"] - p.alpha)
        floor = mp.power(10, -D + plan["gate_guard_digits"]) * (1 + norm(p.alpha))
        diffs = []
        for k, (xl, xh) in enumerate(zip(low["history"], high["history"])):
            e = norm(xh - p.alpha)
            d = norm(xl - xh) / max(e, floor)
            diffs.append({
                "iteration": k,
                "relative_difference_or_floor_scaled": text(d),
                "floor_censored": bool(e <= floor),
                "pass": bool(d < mp.mpf(plan["gate_relative_limit"])),
            })
        passed = (
            low["status"] == "converged"
            and high["status"] == "converged"
            and low["cycles"] == high["cycles"]
            and rr <= tol
            and all(t["pass"] for t in diffs)
        )
        return {
            "group": gid(g),
            "method": key,
            "passed": bool(passed),
            "working_dps": D,
            "verification_dps": V,
            "checked_residual_of_low_iterate": text(rr),
            "checked_error_of_low_iterate": text(ee),
            "low_cycles": low["cycles"],
            "high_cycles": high["cycles"],
            "active_outer_cycles": low["active_outer_cycles"],
            "iterate_comparison": diffs,
            "feasibility_durations_not_comparative_timings": timing,
            "low": low_doc,
            "high": high_doc,
        }


def freeze(out, plan_path=None, neutral_note=""):
    out = Path(out)
    path = out / "protocol_frozen.json"
    if path.exists():
        raise ValueError("Protocol already frozen; do not overwrite")
    if any((out / d).exists() for d in ("sessions_initial", "sessions_extension")):
        raise ValueError("Cannot freeze over prior observations")
    plan = experiment_plan() if plan_path is None else json.loads(Path(plan_path).read_text())
    validate_plan(plan)
    if plan_path is not None and digest(plan["groups"]) != digest(experiment_plan()["groups"]) and not neutral_note.strip():
        raise ValueError("A modified grid requires a prospectively stated neutral design/resource reason")
    plan["status"] = "FROZEN_BEFORE_REPEATED_TIMINGS"
    plan["freeze_utc"] = stamp()
    plan["neutral_design_note"] = neutral_note or "All 16 pre-proposed external-comparator groups retained"
    plan["code_sha256"] = code_hash()
    plan["input_sha256"] = prepare(out, plan)
    plan["protocol_sha256"] = digest(plan)
    js(path, plan)
    js(out / "freeze_environment.json", env())
    return plan


def frozen(out):
    p = json.loads((Path(out) / "protocol_frozen.json").read_text())
    wanted = p.pop("protocol_sha256")
    if digest(p) != wanted:
        raise ValueError("Frozen protocol was modified")
    p["protocol_sha256"] = wanted
    if p["code_sha256"] != code_hash():
        raise ValueError("Code changed since freeze; use a new protocol/output")
    if p["input_sha256"] != input_hash(out):
        raise ValueError("Frozen inputs changed")
    return p


def describe(plan):
    initial = plan["initial_repetitions"] * sum(len(g["methods"]) for g in plan["groups"])
    extra_per_flagged = (plan["extension_target_repetitions"] - plan["initial_repetitions"]) * 2
    return {
        "groups": len(plan["groups"]),
        "method_configurations": sum(len(g["methods"]) for g in plan["groups"]),
        "initial_timed_full_solves": initial,
        "maximum_total_timed_full_solves_if_all_groups_extend": plan["extension_target_repetitions"] * sum(len(g["methods"]) for g in plan["groups"]),
        "additional_timed_full_solves_per_extended_group": extra_per_flagged,
        "working_dps": plan["working_dps"],
        "verification_dps": plan["verification_dps"],
    }


def _bootstrap_summary(data, plan, g, baseline, count, seed_tag):
    b = np.asarray(data[baseline], dtype=float)
    rng = np.random.default_rng(int(digest([plan["seed"], gid(g), seed_tag, count])[:16], 16))
    inds = rng.integers(0, count, (plan["bootstrap_resamples"], count))
    result = {}
    for key, values in data.items():
        a = np.asarray(values, dtype=float)
        ma, mb = float(np.median(a)), float(np.median(b))
        ia = float(np.percentile(a, 75) - np.percentile(a, 25))
        ib = float(np.percentile(b, 75) - np.percentile(b, 25))
        boot = np.median(a[inds], axis=1) / np.median(b[inds], axis=1)
        low, high = np.percentile(boot, [2.5, 97.5])
        result[key] = {
            "median_s": ma,
            "iqr_s": ia,
            "median_ratio": ma / mb,
            "paired_bootstrap_ratio_low": float(low),
            "paired_bootstrap_ratio_high": float(high),
            "two_iqr_filter": bool(abs(ma - mb) > 2 * (ia + ib)),
        }
    return result


def summary_rows(raw, plan, g, expected_count, stage_label):
    data = {key: [] for key in g["methods"]}
    for row in raw:
        data[row["method"]].append(row["elapsed_ns"] * 1e-9)
    if any(len(v) != expected_count for v in data.values()):
        raise ValueError(f"Expected {expected_count} observations per method")
    baseline = g["methods"][-1]
    stats = _bootstrap_summary(data, plan, g, baseline, expected_count, stage_label)
    rows = []
    for key in g["methods"]:
        example = next(r for r in raw if r["method"] == key)
        c = example["counts"]
        rows.append({
            "group": gid(g),
            "stage": stage_label,
            "suite": g["suite"],
            "field": g["field"],
            "n": g["n"],
            "delta": g["delta"],
            "tolerance_digits": g["digits"],
            "method": key,
            "baseline": baseline,
            "repetitions": expected_count,
            **stats[key],
            "cycles": example["cycles"],
            "active_outer_cycles": example["active_outer_cycles"],
            "reference_eta": eta(key, g["field"], g["n"]),
            "reference_work": work(c, g["field"], g["n"]),
            **c,
        })
    return rows


def _load_jsonl(path):
    rows = []
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def _initial_marker(out, id):
    return Path(out) / "initial_completed" / f"{id}.json"


def _extension_marker(out, id):
    return Path(out) / "extension_completed" / f"{id}.json"


def _run_group(out, plan, g, phase, repetitions, repeat_offset, initial_reference=None):
    out = Path(out)
    session = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    root = out / ("sessions_initial" if phase == "initial" else "sessions_extension") / session
    js(root / "environment.json", env())
    id = gid(g)
    gd = root / id
    gd.mkdir(parents=True)

    print(f"Precision gate ({phase}): {id}", flush=True)
    gates = {}
    failed = False
    for key in g["methods"]:
        try:
            r = precision_gate(out, g, key, plan)
            gates[key] = r
            gzout(gd / f"gate_{key}.json.gz", r)
            if not r["passed"]:
                failed = True
        except Exception as e:
            failed = True
            js(gd / f"gate_{key}_error.json", {"error": repr(e)})
    if failed:
        js(gd / "status.json", {"status": "GATE_FAILED_NO_TIMINGS", "group": g, "phase": phase, "session": session})
        raise ArithmeticError(f"Precision gate failed for {id}; inspect {gd}")

    if initial_reference is not None:
        for key in g["methods"]:
            ref = initial_reference[key]
            if gates[key]["low_cycles"] != ref["low_cycles"] or gates[key]["low"]["counts"] != ref["low"]["counts"]:
                js(gd / "status.json", {"status": "EXTENSION_GATE_MISMATCH", "group": g, "phase": phase, "session": session})
                raise ArithmeticError(f"Extension precision gate disagrees with the initial gated path for {id}/{key}")

    with mp.workdps(plan["working_dps"]):
        p = DenseProblem(g["field"], g["n"], out / "inputs")
        x0 = x0_from(out, g)
        tol = mp.mpf(f"1e-{g['digits']}")
        for _ in range(plan["warmups"]):
            for key in g["methods"]:
                iterate(p, x0, method(key), tol, plan["max_cycles"])
        raw = []
        rng = random.Random(int(digest([plan["seed"], id, phase])[:16], 16))
        logfile = gd / "observations.jsonl"
        for j in range(repetitions):
            rep = repeat_offset + j
            keys = list(g["methods"])
            rng.shuffle(keys)
            for position, key in enumerate(keys):
                start = time.perf_counter_ns()
                r = iterate(p, x0, method(key), tol, plan["max_cycles"])
                elapsed = time.perf_counter_ns() - start
                valid = r["status"] == "converged" and r["cycles"] == gates[key]["low_cycles"] and r["counts"] == gates[key]["low"]["counts"]
                row = {
                    "session": session,
                    "phase": phase,
                    "group": id,
                    "repeat": rep,
                    "position": position,
                    "method": key,
                    "elapsed_ns": elapsed,
                    "valid": valid,
                    "status": r["status"],
                    "cycles": r["cycles"],
                    "active_outer_cycles": r["active_outer_cycles"],
                    "active_predictors": r["active_predictors"],
                    "counts": r["counts"],
                    "residual": text(r["residual"]),
                    "error": text(norm(r["x"] - p.alpha)),
                }
                with logfile.open("a", encoding="utf-8") as f:
                    f.write(json.dumps(row) + "\n")
                    f.flush()
                raw.append(row)
                if not valid:
                    raise ArithmeticError("Timing result diverged from the precision-gated trajectory; group incomplete")
            print(f"  {id}: {phase} repetition {j+1}/{repetitions}", flush=True)

    rows = summary_rows(raw, plan, g, repetitions, phase)
    csvout(gd / "summary.csv", rows)
    return session, gd, gates, raw, rows


def run_initial(out, max_groups=None, suite=None):
    out = Path(out)
    plan = frozen(out)
    completed = 0
    groups = [g for g in plan["groups"] if suite is None or g["suite"] == suite]
    (out / "initial_completed").mkdir(exist_ok=True)
    for g in groups:
        id = gid(g)
        marker = _initial_marker(out, id)
        if marker.exists():
            continue
        session, gd, gates, raw, rows = _run_group(
            out,
            plan,
            g,
            "initial",
            plan["initial_repetitions"],
            0,
        )
        js(marker, {
            "status": "INITIAL_COMPLETE",
            "session": session,
            "directory": str(gd.relative_to(out)),
            "protocol_sha256": plan["protocol_sha256"],
            "group": g,
        })
        completed += 1
        summarize(out)
        if max_groups is not None and completed >= max_groups:
            break
    summarize(out)


def _read_initial_gates(out, marker, g):
    root = Path(out) / marker["directory"]
    gates = {}
    for key in g["methods"]:
        with gzip.open(root / f"gate_{key}.json.gz", "rt", encoding="utf-8") as f:
            gates[key] = json.load(f)
    return gates


def _group_raw_from_marker(out, marker):
    return _load_jsonl(Path(out) / marker["directory"] / "observations.jsonl")


def _decision_from_initial_rows(rows, g):
    primary = g["methods"][0]
    r = next(x for x in rows if x["method"] == primary)
    ratio = float(r["median_ratio"])
    low = float(r["paired_bootstrap_ratio_low"])
    high = float(r["paired_bootstrap_ratio_high"])
    two_iqr = bool(r["two_iqr_filter"])
    bootstrap_clear = high < 1.0 or low > 1.0
    direction_consistent = (ratio < 1.0 and high < 1.0) or (ratio > 1.0 and low > 1.0)
    early_stop = bool(two_iqr and bootstrap_clear and direction_consistent)
    reasons = []
    if not two_iqr:
        reasons.append("two_iqr_filter_not_met")
    if low <= 1.0 <= high:
        reasons.append("bootstrap_interval_contains_one")
    if bootstrap_clear and not direction_consistent:
        reasons.append("direction_mismatch")
    if early_stop:
        reasons = ["both_diagnostics_clear_same_direction_after_15"]
    return {
        "group": gid(g),
        "primary_method": primary,
        "baseline": g["methods"][1],
        "median_ratio_after_15": ratio,
        "bootstrap_low_after_15": low,
        "bootstrap_high_after_15": high,
        "two_iqr_filter_after_15": two_iqr,
        "bootstrap_clear_after_15": bootstrap_clear,
        "direction_consistent_after_15": direction_consistent,
        "extend_to_31": not early_stop,
        "decision_reasons": reasons,
    }


def decide_extensions(out):
    out = Path(out)
    plan = frozen(out)
    decision_path = out / "extension_decisions.json"
    initial_markers = sorted((out / "initial_completed").glob("*.json")) if (out / "initial_completed").exists() else []
    if len(initial_markers) != len(plan["groups"]):
        raise ValueError("Complete all 16 initial groups before freezing extension decisions")

    decisions = []
    for g in plan["groups"]:
        marker = json.loads(_initial_marker(out, gid(g)).read_text())
        if marker["protocol_sha256"] != plan["protocol_sha256"]:
            raise ValueError("Mixed protocol hashes")
        raw = _group_raw_from_marker(out, marker)
        rows = summary_rows(raw, plan, g, plan["initial_repetitions"], "initial")
        decisions.append(_decision_from_initial_rows(rows, g))

    doc = {
        "schema": 1,
        "created_utc": stamp(),
        "protocol_sha256": plan["protocol_sha256"],
        "rule": plan["adaptive_rule"],
        "initial_repetitions": plan["initial_repetitions"],
        "extension_target_repetitions": plan["extension_target_repetitions"],
        "decisions": decisions,
    }
    if decision_path.exists():
        old = json.loads(decision_path.read_text())
        # Ignore timestamp for deterministic re-check.
        a = {k: v for k, v in old.items() if k != "created_utc"}
        b = {k: v for k, v in doc.items() if k != "created_utc"}
        if a != b:
            raise ValueError("Existing extension decisions do not match the mechanically recomputed rule")
        return old
    js(decision_path, doc)
    csvout(out / "extension_decisions.csv", decisions)
    summarize(out)
    return doc


def _decision_map(out):
    path = Path(out) / "extension_decisions.json"
    if not path.exists():
        raise ValueError("Run the decision command after all initial groups are complete")
    doc = json.loads(path.read_text())
    plan = frozen(out)
    if doc["protocol_sha256"] != plan["protocol_sha256"]:
        raise ValueError("Decision/protocol mismatch")
    return {d["group"]: d for d in doc["decisions"]}


def run_extensions(out, max_groups=None, suite=None):
    out = Path(out)
    plan = frozen(out)
    decisions = _decision_map(out)
    (out / "extension_completed").mkdir(exist_ok=True)
    completed = 0
    extra = plan["extension_target_repetitions"] - plan["initial_repetitions"]
    groups = [g for g in plan["groups"] if suite is None or g["suite"] == suite]
    for g in groups:
        id = gid(g)
        if not decisions[id]["extend_to_31"]:
            continue
        marker = _extension_marker(out, id)
        if marker.exists():
            continue
        initial_marker = json.loads(_initial_marker(out, id).read_text())
        initial_gates = _read_initial_gates(out, initial_marker, g)
        session, gd, gates, raw, rows = _run_group(
            out,
            plan,
            g,
            "extension",
            extra,
            plan["initial_repetitions"],
            initial_reference=initial_gates,
        )
        js(marker, {
            "status": "EXTENSION_COMPLETE",
            "session": session,
            "directory": str(gd.relative_to(out)),
            "protocol_sha256": plan["protocol_sha256"],
            "group": g,
            "added_repetitions": extra,
        })
        completed += 1
        summarize(out)
        if max_groups is not None and completed >= max_groups:
            break
    summarize(out)


def _final_raw_for_group(out, plan, g, decisions):
    id = gid(g)
    im = json.loads(_initial_marker(out, id).read_text())
    raw = _group_raw_from_marker(out, im)
    status = "EARLY_STOP_15"
    if decisions and decisions[id]["extend_to_31"]:
        em = _extension_marker(out, id)
        if not em.exists():
            return raw, "EXTENSION_PENDING"
        ext = json.loads(em.read_text())
        raw += _group_raw_from_marker(out, ext)
        status = "EXTENDED_TO_31"
    return raw, status


def summarize(out):
    out = Path(out)
    plan = frozen(out)
    initial_rows = []
    for g in plan["groups"]:
        marker_path = _initial_marker(out, gid(g))
        if not marker_path.exists():
            continue
        marker = json.loads(marker_path.read_text())
        raw = _group_raw_from_marker(out, marker)
        initial_rows.extend(summary_rows(raw, plan, g, plan["initial_repetitions"], "initial"))
    csvout(out / "combined_initial_summary.csv", initial_rows)

    decisions = None
    if (out / "extension_decisions.json").exists():
        decisions = _decision_map(out)

    final_rows = []
    group_status = []
    for g in plan["groups"]:
        id = gid(g)
        if not _initial_marker(out, id).exists():
            group_status.append({"group": id, "status": "INITIAL_PENDING"})
            continue
        raw, status = _final_raw_for_group(out, plan, g, decisions)
        expected = len(raw) // len(g["methods"])
        if status != "EXTENSION_PENDING":
            rows = summary_rows(raw, plan, g, expected, "final")
            for r in rows:
                r["adaptive_status"] = status
            final_rows.extend(rows)
        group_status.append({"group": id, "status": status, "repetitions_available": expected})
    csvout(out / "combined_final_summary.csv", final_rows)
    csvout(out / "group_status.csv", group_status)

    initial_complete = sum(_initial_marker(out, gid(g)).exists() for g in plan["groups"])
    flagged = 0
    extension_complete = 0
    if decisions:
        flagged = sum(bool(d["extend_to_31"]) for d in decisions.values())
        extension_complete = sum(
            _extension_marker(out, gid(g)).exists()
            for g in plan["groups"]
            if decisions[gid(g)]["extend_to_31"]
        )
    final_complete = initial_complete == len(plan["groups"]) and decisions is not None and extension_complete == flagged
    js(out / "coverage.json", {
        "planned_groups": len(plan["groups"]),
        "initial_complete_groups": initial_complete,
        "extension_decisions_frozen": decisions is not None,
        "groups_flagged_for_extension": flagged,
        "extensions_complete": extension_complete,
        "final_complete": final_complete,
        "rule": plan["adaptive_rule"],
    })
    return final_rows


def status(out):
    out = Path(out)
    plan = frozen(out)
    summarize(out)
    coverage = json.loads((out / "coverage.json").read_text())
    print(json.dumps(coverage, indent=2))
    if (out / "extension_decisions.json").exists():
        doc = json.loads((out / "extension_decisions.json").read_text())
        print("\nExtension decisions:")
        for d in doc["decisions"]:
            print(f"  {d['group']}: {('EXTEND TO '+str(plan['extension_target_repetitions'])) if d['extend_to_31'] else ('STOP AT '+str(plan['initial_repetitions']))}; {', '.join(d['decision_reasons'])}")


def self_test(out="validation/self_test.json"):
    checks = []
    plan = experiment_plan()
    validate_plan(plan)
    checks.append({"name": "plan_validates", "pass": True})
    checks.append({"name": "groups_16", "pass": len(plan["groups"]) == 16})
    checks.append({"name": "initial_480", "pass": describe(plan)["initial_timed_full_solves"] == 480})
    checks.append({"name": "max_992", "pass": describe(plan)["maximum_total_timed_full_solves_if_all_groups_extend"] == 992})

    # Synthetic decision checks, independent of machine timings.
    g = plan["groups"][0]
    clear = [
        {"method": g["methods"][0], "median_ratio": 0.90, "paired_bootstrap_ratio_low": 0.88, "paired_bootstrap_ratio_high": 0.93, "two_iqr_filter": True},
        {"method": g["methods"][1], "median_ratio": 1.0, "paired_bootstrap_ratio_low": 1.0, "paired_bootstrap_ratio_high": 1.0, "two_iqr_filter": False},
    ]
    amb_boot = [
        {"method": g["methods"][0], "median_ratio": 0.98, "paired_bootstrap_ratio_low": 0.95, "paired_bootstrap_ratio_high": 1.02, "two_iqr_filter": True},
        clear[1],
    ]
    amb_iqr = [
        {"method": g["methods"][0], "median_ratio": 1.08, "paired_bootstrap_ratio_low": 1.04, "paired_bootstrap_ratio_high": 1.12, "two_iqr_filter": False},
        clear[1],
    ]
    checks.append({"name": "clear_stops_15", "pass": _decision_from_initial_rows(clear, g)["extend_to_31"] is False})
    checks.append({"name": "bootstrap_ambiguity_extends", "pass": _decision_from_initial_rows(amb_boot, g)["extend_to_31"] is True})
    checks.append({"name": "iqr_ambiguity_extends", "pass": _decision_from_initial_rows(amb_iqr, g)["extend_to_31"] is True})
    passed = all(c["pass"] for c in checks)
    js(out, {"passed": passed, "checks": checks, "code_sha256": code_hash()})
    if not passed:
        raise AssertionError(checks)
    return checks
