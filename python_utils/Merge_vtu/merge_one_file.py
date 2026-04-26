#!/usr/bin/env pvpython
import os, glob, sys, argparse
import paraview.simple as pv

def pack_vtu_series(root: str):
    root = os.path.abspath(root)
    merged_dir = os.path.join(root, "merged")

    root_vtus = sorted(glob.glob(os.path.join(root, "Output_t*.vtu")))
    merged_vtus = sorted(glob.glob(os.path.join(merged_dir, "Output_t*.vtu")))

    if merged_vtus:
        directory = merged_dir
        vtu_files = merged_vtus
        source = "merged"
    elif root_vtus:
        directory = root
        vtu_files = root_vtus
        source = "root"
    else:
        print("ERROR: no VTU files found.")
        print("Checked:")
        print(f"  root:   {root}")
        print(f"  merged: {merged_dir}")
        sys.exit(1)

    print(f"[2] Using VTU source: {source}")
    print(f"[2] VTU directory: {directory}")
    print(f"[2] Found {len(vtu_files)} VTU files")

    # Important: save PVD in root so outputs are centralized.
    pvd_path = os.path.join(root, "Simulation.pvd")

    with open(pvd_path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <Collection>\n')

        for i, fp in enumerate(vtu_files):
            rel_name = os.path.relpath(fp, start=root)
            rel_name = rel_name.replace(os.sep, "/")

            f.write(
                f'    <DataSet timestep="{i}" group="" part="0" file="{rel_name}"/>\n'
            )

        f.write('  </Collection>\n')
        f.write('</VTKFile>\n')

    print(f"[2] PVD created: {pvd_path}")

    # Load time series in ParaView
    pv.ResetSession()

    print(f"[3] Loading PVD time series: {pvd_path}")
    reader = pv.OpenDataFile(pvd_path)
    pv.UpdatePipeline(proxy=reader)

    print("[3] Merging blocks")
    merged = pv.MergeBlocks(Input=reader)
    pv.UpdatePipeline(proxy=merged)

    out_vtkhdf = os.path.join(root, "Simulation.vtkhdf")
    out_xmf = os.path.join(root, "Simulation.xmf")

    try:
        print(f"[4] Writing VTKHDF: {out_vtkhdf}")
        pv.SaveData(out_vtkhdf, proxy=merged, WriteAllTimeSteps=1)
        print("[4] SUCCESS VTKHDF")
        final_output = out_vtkhdf

    except Exception as e:
        print(f"[4] VTKHDF failed: {e}")
        print(f"[4] Writing fallback XDMF: {out_xmf}")
        pv.SaveData(out_xmf, proxy=merged, WriteAllTimeSteps=1)
        print("[4] SUCCESS XDMF")
        final_output = out_xmf

    # Write deletion manifest only after successful compact export.
    delete_list_path = os.path.join(root, "vtu_files_to_delete.txt")

    with open(delete_list_path, "w") as f:
        for fp in vtu_files:
            f.write(fp + "\n")

    print(f"[5] Deletion list written for {source} VTU files: {delete_list_path}")

    return final_output, pvd_path, delete_list_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True) # root is outside the merged directory.
    args = ap.parse_args()
    pack_vtu_series(root=args.root)


if __name__ == "__main__":
    main()