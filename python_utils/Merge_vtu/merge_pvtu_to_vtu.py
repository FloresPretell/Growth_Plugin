#!/usr/bin/env pvpython
import os, glob, re, sys
import argparse
import paraview.simple as pv


def merge_pvtu_to_vtu(root=None, list_to_delete=None):
    merged_dir = os.path.join(root, "merged")
    os.makedirs(merged_dir, exist_ok=True) #create merge directory

    # ====== STEP 1: merge particiones por timestep ======
    pvtu_pattern = os.path.join(root, "Output_t*.pvtu")
    pvtu_files = sorted(glob.glob(pvtu_pattern))
    print(f"[1] Found {len(pvtu_files)} pvtu files")

    name_list = open(list_to_delete, "w") if list_to_delete else None # open list to register files to delete after merge

    pv.ResetSession()
    for pvtu in pvtu_files:
        
        stem = os.path.basename(pvtu).replace(".pvtu", "")     #remove the .pvtu extension
        base = stem + ".vtu" 
        out_vtu = os.path.join(merged_dir, base)
        if os.path.exists(out_vtu):
            continue
        print("[merge]", pvtu, "->", out_vtu)

        reader = pv.XMLPartitionedUnstructuredGridReader(FileName=[pvtu])
        pv.UpdatePipeline()

        m = pv.MergeBlocks(Input=reader)
        pv.UpdatePipeline()
        pv.SaveData(out_vtu, proxy=m)
        
        # si se creó el output, registramos qué borrar
        merged_count = 0
        if os.path.exists(out_vtu) and os.path.getsize(out_vtu) > 0:
            merged_count += 1
            if name_list:
                name_list.write(pvtu + "\n")
                parts = sorted(glob.glob(os.path.join(root, stem + "_*.vtu")))
                for p in parts:
                    name_list.write(p + "\n")
                    
    if name_list:
        name_list.close()
    print("[1] Done merging partitions")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Folder containing Output_t*.pvtu")
    ap.add_argument("--list_to_delete", default="delete_after_merge.txt")
    args = ap.parse_args()
    
    merge_pvtu_to_vtu(root=args.root, list_to_delete=args.list_to_delete)
    
if __name__ == "__main__":
    main()












# # ====== STEP 2: crear PVD con todos los tiempos ======
# vtu_files = sorted(glob.glob(os.path.join(merged_dir, "ResultadosUnidos_t*.vtu")))
# if not vtu_files:
#     print("ERROR: no merged vtu files found")
#     sys.exit(1)

# def t_of(name):
#     m = re.search(r"_t0*([0-9]+)\.vtu$", os.path.basename(name))
#     return int(m.group(1)) if m else 0

# pvd_path = os.path.join(merged_dir, "simulation.pvd")
# with open(pvd_path, "w") as f:
#     f.write('<?xml version="1.0"?>\n')
#     f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
#     f.write('  <Collection>\n')
#     for fp in vtu_files:
#         f.write(f'    <DataSet timestep="{t_of(fp)}" group="" part="0" file="{os.path.basename(fp)}"/>\n')
#     f.write('  </Collection>\n')
#     f.write('</VTKFile>\n')

# print("[2] PVD created:", pvd_path)

# # ====== STEP 3: pack tiempos en 1 (o 2) archivos ======
# pv.ResetSession()
# sim = OpenDataFile(pvd_path)
# UpdatePipeline()

# # intenta 1 archivo VTKHDF
# try:
#     print("[3] Trying VTKHDF ->", out_vtkhdf)
#     SaveData(out_vtkhdf, proxy=sim)
#     print("[3] SUCCESS VTKHDF")
# except Exception as e:
#     print("[3] VTKHDF failed:", e)
#     print("[3] Trying XDMF ->", out_xmf)
#     SaveData(out_xmf, proxy=sim)
#     print("[3] SUCCESS XDMF (xmf+h5)")

# print("DONE") 