#!/usr/bin/env python3
"""Inventory of the Biopython API surface.

Walks every Bio.* submodule, imports it, and dumps:
- classes / functions / methods / constants with signatures and docstrings;
- deprecated flags (from import warnings and docstring text);
- module size in lines;
- import errors (optional dependencies).

Output files (next to this script):
- biopython_inventory.csv  -- one public API object per row;
- biopython_inventory.json -- full dump, including private members;
- INDEX.md                 -- per-package summary.
"""
import csv
import inspect
import json
import sys
import time
import warnings
from datetime import date, datetime
from pathlib import Path

import Bio
import importlib
import pkgutil

OUT_DIR = Path(__file__).resolve().parent

# Dunders that carry an API contract -- keep those, skip the rest.
DUnder_OK = {"__init__", "__iter__", "__next__", "__call__", "__getitem__",
             "__setitem__", "__len__", "__enter__", "__exit__", "__eq__"}


def first_line(text, limit=120):
    if not text:
        return ""
    return text.strip().splitlines()[0].strip()[:limit]


def module_loc(mod):
    try:
        path = inspect.getsourcefile(mod)
        if path and path.endswith(".py"):
            with open(path, encoding="utf-8", errors="replace") as f:
                return sum(1 for _ in f)
    except (TypeError, OSError):
        pass
    return 0  # C extension or namespace package


def classify(val):
    if inspect.ismodule(val):
        return None
    if inspect.isclass(val):
        return "exception" if issubclass(val, Exception) else "class"
    if inspect.isroutine(val) or callable(val):
        return "function"
    return "constant"


def safe_signature(val):
    try:
        return str(inspect.signature(val))
    except (ValueError, TypeError):
        return ""


def method_rows(cls):
    """Own methods of a class (from __dict__), one row each."""
    rows = []
    for name, val in vars(cls).items():
        if name.startswith("__") and name.endswith("__") and name not in DUnder_OK:
            continue
        func = val
        kind = "method"
        if isinstance(val, (classmethod, staticmethod)):
            func = val.__func__
            kind = "classmethod" if isinstance(val, classmethod) else "staticmethod"
        elif isinstance(val, property):
            kind = "property"
            func = val.fget
        elif not (inspect.isfunction(val) or inspect.isroutine(val)):
            kind = "attribute"
        sig = "(property)" if kind == "property" else safe_signature(func)
        doc = inspect.getdoc(func) if not inspect.isdatadescriptor(val) else inspect.getdoc(val)
        rows.append({
            "name": name,
            "qualname": f"{cls.__name__}.{name}",
            "kind": kind,
            # Imported from another module (e.g. numpy helpers) -> not this API.
            "origin": "defined",
            "visibility": "private" if name.startswith("_") else "public",
            "signature": sig,
            "deprecated": "deprecat" in (doc or "").lower(),
            "doc": first_line(doc),
        })
    return rows


def main():
    t0 = time.perf_counter()
    mods = sorted(
        (m.name, m.ispkg) for m in pkgutil.walk_packages(Bio.__path__, "Bio.")
        if not m.name.startswith("Bio._")
    )
    print(f"biopython {Bio.__version__}, python {sys.version.split()[0]}, "
          f"modules to scan: {len(mods)}")
    print(f"scan started {datetime.now():%Y-%m-%d %H:%M:%S}")

    index = {}
    for name, ispkg in mods:
        entry = {"module": name, "package": name.split(".")[1] if "." in name else "Bio",
                 "ispkg": ispkg, "loc": 0, "all_exports": [], "import_error": "",
                 "deprecation_at_import": False, "members": []}
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                mod = importlib.import_module(name)
            except Exception as exc:  # optional deps: reportlab, DB drivers, etc.
                entry["import_error"] = f"{type(exc).__name__}: {exc}"
                index[name] = entry
                continue
        entry["deprecation_at_import"] = any(
            w.category in (DeprecationWarning, PendingDeprecationWarning, FutureWarning)
            for w in caught)
        entry["loc"] = module_loc(mod)
        entry["all_exports"] = sorted(getattr(mod, "__all__", []) or [])
        all_names = set(entry["all_exports"])

        for mname, val in vars(mod).items():
            if mname.startswith("__") and mname.endswith("__"):
                continue
            # Keep only what this module defines or re-exports via __all__.
            defined = getattr(val, "__module__", None) == name
            exported = mname in all_names
            if not (defined or exported) or inspect.ismodule(val):
                continue
            kind = classify(val)
            doc = inspect.getdoc(val)
            row = {"name": mname, "qualname": mname, "kind": kind or "unknown",
                   "visibility": "private" if mname.startswith("_") else "public",
                   "origin": "defined" if defined else "re-exported",
                   "signature": safe_signature(val) if kind == "function" else "",
                   "deprecated": "deprecat" in (doc or "").lower(),
                   "doc": first_line(doc)}
            entry["members"].append(row)
            if kind == "class":
                entry["members"].extend(method_rows(val))
        index[name] = entry

    # --- CSV: public objects, one per row ---
    csv_path = OUT_DIR / "biopython_inventory.csv"
    cols = ["package", "module", "name", "qualname", "kind", "visibility",
            "origin", "signature", "deprecated", "doc"]
    n_csv = 0
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for entry in index.values():
            for r in entry["members"]:
                if r["visibility"] != "public" or r["kind"] == "attribute":
                    continue
                w.writerow({"package": entry["package"], "module": entry["module"],
                            **{k: r[k] for k in cols if k not in ("package", "module")}})
                n_csv += 1

    # --- JSON: full dump ---
    json_path = OUT_DIR / "biopython_inventory.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump({"biopython_version": Bio.__version__, "generated": date.today().isoformat(),
                   "modules": index}, f, ensure_ascii=False, indent=1)

    # --- summaries ---
    totals = {"class": 0, "exception": 0, "function": 0, "method": 0,
              "classmethod": 0, "staticmethod": 0, "property": 0, "constant": 0}
    failed, deprecated = [], []
    for entry in index.values():
        if entry["import_error"]:
            failed.append((entry["module"], entry["import_error"]))
        if entry["deprecation_at_import"]:
            deprecated.append(entry["module"])
        for r in entry["members"]:
            if r["kind"] in totals:
                totals[r["kind"]] += 1
    by_pkg = {}
    for entry in index.values():
        p = by_pkg.setdefault(entry["package"], {"modules": 0, "classes": 0,
                              "functions": 0, "methods": 0, "constants": 0, "loc": 0})
        p["modules"] += 1
        p["loc"] += entry["loc"]
        for r in entry["members"]:
            if r["kind"] in ("class", "exception"):
                p["classes"] += 1
            elif r["kind"] == "function":
                p["functions"] += 1
            elif r["kind"] in ("method", "classmethod", "staticmethod", "property"):
                p["methods"] += 1
            else:
                p["constants"] += 1

    md = [f"# Biopython {Bio.__version__} API inventory",
          f"Generated {date.today().isoformat()} by `scan_biopython.py`.",
          "",
          f"**Totals:** {len(index)} modules ({len(failed)} failed to import), "
          f"{sum(e['loc'] for e in index.values()):,} lines of code, "
          f"{totals['class'] + totals['exception']} classes, {totals['function']} functions, "
          f"{totals['method'] + totals['classmethod'] + totals['staticmethod']} methods, "
          f"{totals['property']} properties, {totals['constant']} constants.",
          f"Public objects in CSV: {n_csv:,}.", "",
          "| Package | Modules | Classes | Functions | Methods | Constants | LOC |",
          "|---|---:|---:|---:|---:|---:|---:|"]
    for pkg, p in sorted(by_pkg.items(),
                         key=lambda kv: -(kv[1]["functions"] + kv[1]["methods"] + kv[1]["classes"])):
        md.append(f"| {pkg} | {p['modules']} | {p['classes']} | {p['functions']} | "
                  f"{p['methods']} | {p['constants']} | {p['loc']:,} |")
    if deprecated:
        md += ["", "## Deprecation warnings at import", ""]
        md += [f"- `{m}`" for m in sorted(deprecated)]
    if failed:
        md += ["", "## Failed imports (optional dependencies)", ""]
        md += [f"- `{m}` — {err}" for m, err in sorted(failed)]
    (OUT_DIR / "INDEX.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print(f"done in {time.perf_counter() - t0:.1f}s -> {csv_path.name}, "
          f"{json_path.name}, INDEX.md | public objects: {n_csv:,}")
    print("totals:", totals)


if __name__ == "__main__":
    main()