#!/usr/bin/env pvpython
import os, glob, re, sys
import argparse
import paraview.simple as pv


def extract_timestep(path: str):
    m = re.search(r"Output_t(\d+)\.pvtu$", os.path.basename(path))
    return int(m.group(1)) if m else None


def extract_pvtu_files(root: str, t0=None, t1=None, snapshot_interval=None):
    """
    Select pvtu files in [t0, t1] and consistent with snapshot_interval.
    Assumes naming: Output_t<digits>.pvtu
    """
    pvtu_pattern = os.path.join(root, "Output_t*.pvtu")
    pvtu_files_all = sorted(glob.glob(pvtu_pattern))

    pvtu_files = []
    for pvtu in pvtu_files_all:
        t = extract_timestep(pvtu)
        if t is None:
            continue

        if t0 is not None and t < t0:
            continue
        if t1 is not None and t > t1:
            continue

        # UG4 prints only when: step % snapshot_interval == 0
        if snapshot_interval is not None and snapshot_interval > 0:
            if (t % snapshot_interval) != 0:
                continue

        pvtu_files.append(pvtu)

    print(f"[select] Found {len(pvtu_files_all)} pvtu files total")
    print(f"[select] Selected {len(pvtu_files)} pvtu files (t0={t0}, t1={t1}, snap={snapshot_interval})")
    return pvtu_files
 
 
def merge_pvtu_to_vtu(root=None, list_to_delete=None,t0=None, t1=None, snapshot_interval=None):
    merged_dir = os.path.join(root, "merged")
    os.makedirs(merged_dir, exist_ok=True) #create merge directory

    pvtu_files = extract_pvtu_files(root=root, t0=t0, t1=t1, snapshot_interval=snapshot_interval)
    name_list = open(list_to_delete, "w") if list_to_delete else None # open list to register files to delete after merge

    pv.ResetSession()
    for pvtu in pvtu_files:
        stem = os.path.basename(pvtu).replace(".pvtu", "")     #remove the .pvtu extension
        base = stem + ".vtu" 
        out_vtu = os.path.join(merged_dir, base)
        if os.path.exists(out_vtu) and os.path.getsize(out_vtu) > 0:
            if name_list: # merged already exists -> still register pvtu for deletion
                name_list.write(pvtu + "\n")
            continue
        print("[merge]", pvtu, "->", out_vtu)

        reader = pv.XMLPartitionedUnstructuredGridReader(FileName=[pvtu])
        pv.UpdatePipeline()

        m = pv.MergeBlocks(Input=reader)
        pv.UpdatePipeline()
        pv.SaveData(out_vtu, proxy=m)
        
        # si se creó el output, registramos qué borrar
        if os.path.exists(out_vtu) and os.path.getsize(out_vtu) > 0:
            if name_list:
                name_list.write(pvtu + "\n")
                 
    if name_list:
        name_list.close()
    print("[1] Done merging partitions")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Folder containing Output_t*.pvtu")
    ap.add_argument("--list_to_delete", default="delete_after_merge.txt")
    ap.add_argument("--t0", type=int, default=None, help="First timestep index (inclusive)")
    ap.add_argument("--t1", type=int, default=None, help="Last timestep index (inclusive)")
    ap.add_argument("--snapshot_interval", type=int, default=None,
                    help="Keep only timesteps where t %% snapshot_interval == 0")

    args = ap.parse_args()
    
    merge_pvtu_to_vtu(
         root=args.root,
         list_to_delete=args.list_to_delete,
         t0=args.t0,
         t1=args.t1,
         snapshot_interval=args.snapshot_interval
    )
    
if __name__ == "__main__":
    main()











