#!/usr/bin/env pvpython
import os, glob, sys, argparse
import paraview.simple as pv

def pack_vtu_series(root: str):
    directory = os.path.join(root, "merged")
    vtu_files = sorted(glob.glob(os.path.join(directory, "Output_t*.vtu")))
    print(f"[2] Found {len(vtu_files)} merged vtu files")

    if not vtu_files:
        print("ERROR: no VTU files found")
        sys.exit(1)

    # create PVD
    pvd_path = os.path.join(directory, "Simulation.pvd")

    with open(pvd_path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection">\n')
        f.write('<Collection>\n')
        for i, fp in enumerate(vtu_files):
            name = os.path.basename(fp)
            f.write(
                f'<DataSet timestep="{i}" group="" part="0" file="{name}"/>\n'
            )
        f.write('</Collection>\n')
        f.write('</VTKFile>\n')
    print("[2] PVD created:", pvd_path)

    # load time series
    pv.ResetSession()
    reader = pv.OpenDataFile(pvd_path)
    pv.UpdatePipeline()
    
    m = pv.MergeBlocks(Input=reader)
    pv.UpdatePipeline()

    out_vtkhdf = os.path.join(root, "Simulation.vtkhdf")
    out_xmf    = os.path.join(root, "Simulation.xmf")
    try:
        print("[3] Writing VTKHDF:", out_vtkhdf)
        pv.SaveData(out_vtkhdf, proxy=m, WriteAllTimeSteps=1)
        print("[3] SUCCESS VTKHDF")
    except Exception as e:
        print("[3] VTKHDF failed:", e)
        print("[3] Writing XDMF:", out_xmf)
        pv.SaveData(out_xmf, proxy=m, WriteAllTimeSteps=1)
        print("[3] SUCCESS XDMF")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True) # root is outside the merged directory.
    args = ap.parse_args()
    pack_vtu_series(root=args.root)


if __name__ == "__main__":
    main()