#!/usr/bin/env pvpython
"""merge_and_purge_3d.py — per-timestep merge + IMMEDIATE purge for the 3D EVALUATOR.

Why a separate script: the calcium/2D tool (merge_and_purge.py) is hardcoded to the
`Output_t*.pvtu` / `Output_p*_t<N>.vtu` naming. The 3D growth app writes MANY per-field
series instead — `eikanol_*`, `LSF_*`, `Eikonal_Inside_*`, `Results_*`, `Curvatura_Kapla_*`,
`Flujo_Calcio_*`, `Debug_*`, ... — so the 2D tool matches ZERO files here and frees nothing.

This generalizes it: for EVERY `<field>_t<N>.pvtu` in --root, merge its partitions
-> merged/<field>_t<N>.vtu, verify the merged file is non-empty, THEN delete that
pvtu and all its `<field>_p*_t<N>.vtu` partitions. Inodes drop continuously.

Safe: a timestep's fragments are removed only after its merged .vtu exists and is >0 bytes.
Idempotent: re-running skips already-purged timesteps (their pvtu is gone).

Usage:
  source /scratch/flore0a/Modules.sh mesa
  pvbatch merge_and_purge_3d.py --root /scratch/flore0a/3D_evaluator/<run>/growth --dry_run
  pvbatch merge_and_purge_3d.py --root /scratch/flore0a/3D_evaluator/<run>/growth
"""
import os, glob, re, argparse
import paraview.simple as pv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True,
                    help="dir containing <field>_t<N>.pvtu (e.g. .../<run>/growth)")
    ap.add_argument("--dry_run", action="store_true",
                    help="report what WOULD be merged/deleted; touch nothing")
    args = ap.parse_args()
    root = args.root.rstrip("/")
    merged = os.path.join(root, "merged")
    if not args.dry_run:
        os.makedirs(merged, exist_ok=True)

    pvtus = sorted(glob.glob(os.path.join(root, "*_t*.pvtu")))
    fields = {}
    for p in pvtus:
        m = re.match(r"(.+)_t\d+$", os.path.basename(p)[:-5])
        if m:
            fields[m.group(1)] = fields.get(m.group(1), 0) + 1
    print(f"[3dpurge] {root}: {len(pvtus)} pvtu across {len(fields)} fields", flush=True)
    for f, c in sorted(fields.items(), key=lambda kv: -kv[1]):
        print(f"           {c:>5}  {f}", flush=True)
    if args.dry_run:
        frags = len(glob.glob(os.path.join(root, "*_p[0-9]*_t[0-9]*.vtu")))
        print(f"[3dpurge] DRY RUN: would merge {len(pvtus)} timesteps and free "
              f"~{frags + len(pvtus)} inodes (partitions + pvtu)", flush=True)
        return

    pv.ResetSession()
    done = 0
    for pvtu in pvtus:
        stem = os.path.basename(pvtu)[:-5]              # <field>_tXXXX
        m = re.match(r"(.+)_t(\d+)$", stem)
        if not m:
            continue
        field, tnum = m.group(1), m.group(2)
        out_vtu = os.path.join(merged, stem + ".vtu")
        try:
            reader = pv.XMLPartitionedUnstructuredGridReader(FileName=[pvtu])
            pv.UpdatePipeline()
            mb = pv.MergeBlocks(Input=reader)
            pv.UpdatePipeline()
            pv.SaveData(out_vtu, proxy=mb)             # overwrite
            pv.Delete(mb); pv.Delete(reader)
        except Exception as e:                          # noqa: BLE001
            print(f"[3dpurge] MERGE FAIL {pvtu}: {e} -> keep fragments", flush=True)
            continue
        if os.path.exists(out_vtu) and os.path.getsize(out_vtu) > 0:
            try:
                os.remove(pvtu)
            except OSError:
                pass
            for part in glob.glob(os.path.join(root, f"{field}_p*_t{tnum}.vtu")):
                try:
                    os.remove(part)
                except OSError:
                    pass
            done += 1
            if done % 25 == 0:
                print(f"[3dpurge] {os.path.basename(root)} {done}/{len(pvtus)} compacted",
                      flush=True)
        else:
            print(f"[3dpurge] verify FAIL (empty merged) {stem} -> keep fragments", flush=True)
    print(f"[3dpurge] DONE {root}: {done}/{len(pvtus)} timesteps compacted", flush=True)


if __name__ == "__main__":
    main()
