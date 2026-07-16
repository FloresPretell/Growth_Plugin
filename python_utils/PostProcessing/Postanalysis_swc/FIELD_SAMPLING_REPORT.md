# Feasibility Report: SWC-Based Field Sampling
*PostProcessing Stage 5 — Quantitative concentration dataset on neuronal morphology*

---

## 1. What files exist and where

### SWC skeleton files
- **Location pattern**: `<case_dir>/Simulation2D/merged/swc/<frame_idx>/clean_skeleton.swc`
- **Example**: `/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/merged/swc/033/clean_skeleton.swc`
- **File type**: Standard SWC format — but **coordinates are pixel coordinates** from the 2048×2048 rendered PNG (see Section 3 for the transform)
- **One SWC per simulation frame**: D7.5_V3 has 44 frames (000–043), each with `clean_skeleton.swc` and `pixel_skeleton.swc`
- **Two SWC variants**:
  - `clean_skeleton.swc` — branchpoint clusters merged, degree-2 chains compressed, spurs pruned. Contains only topological junction nodes (soma root, branchpoints, tips). **Recommended for sampling.**
  - `pixel_skeleton.swc` — every skeleton pixel retained (~hundreds to thousands of nodes per frame). Better spatial resolution but much larger.

### Simulation field files
- **Location pattern**: `<case_dir>/Simulation2D/Simulation.vtkhdf`
- **Example**: `/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/Simulation.vtkhdf`
- **Format**: VTKHDF (VTK HDF5 transient time series) — single file containing all time steps
- **Fallback**: `<case_dir>/Simulation2D/merged/Output_t*.vtu` (VTU series, one file per step)
- **Verified on D7.5_V3**: 44 time steps, t = 0, 1, 2, …, 43

### Available point-data arrays (confirmed in VTKHDF)

| Field name | Biological meaning | Relevant region |
|------------|-------------------|-----------------|
| `u_t` | Tubulin concentration | Intracellular (lsf < 0) |
| `u_u` | MAP2 free | Intracellular |
| `u_b` | MAP2 bound | Intracellular |
| `u_p` | MAP2 phosphorylated | Intracellular |
| `u_ca_cyt` | Cytosolic calcium | Intracellular |
| `inhibitor` | Inhibitor signal | Extracellular (lsf > 0) |
| `lsf` | Level-set function (< 0 inside cell, > 0 outside) | Everywhere |
| `curvature` | Interface curvature | Near interface |
| `Vel_Tubuline` | Tubulin advection velocity (vector, 3 components) | Intracellular |
| `Flujo_Calcio` | Calcium flux vector | Intracellular |
| `Flujo_inhibicion` | Inhibitor flux vector | Extracellular |
| `GrowthVel` *(cell data)* | Growth velocity at interface | Interface |

All primary biological fields (`u_t`, `u_u`, `u_b`, `u_p`, `u_ca_cyt`, `inhibitor`, `lsf`) are stored as **POINT data**, making them directly interpolatable at arbitrary coordinates using `vtkProbeFilter`.

---

## 2. Time-step alignment

Both sources use the same frame index (integer 0–43):
- VTK time step index `i` → time value `t = float(i)` (confirmed: TimestepValues = [0.0, 1.0, …, 43.0])
- SWC frame folder `NNN` → integer `int("NNN") = i`

Mapping is direct: SWC folder `"033"` → VTK time step index 33, time value 33.0.

---

## 3. Coordinate system — THE CRITICAL ISSUE

This is the most important technical concern. The two data sources use **different coordinate systems**.

### SWC coordinates (pixel space)
The SWC extraction pipeline works entirely in image pixel space.  
Looking at `util_postanalisis_export_swc.py` lines 296–300 and 319–322:

```python
y, x = node           # (row, col) in skeleton image
f.write(f"... {x:.3f} {-y:.3f} 0.000 ...")
```

Therefore:
- `x_swc` = **pixel column** (0 = left edge, 2047 = right edge)
- `y_swc` = **negative pixel row** (0 = top edge, −2047 = bottom edge)
- `z_swc` = 0.000 (always)
- `radius` = 1.000 (always — no actual radius was estimated)

### Simulation coordinates (physical space)
The VTK mesh uses physical coordinates. For D7.5_V3 at t=0:
- x ∈ [−3.0, 3.0]
- y ∈ [−1.550, 3.850]
- z = 0.0 (2D simulation embedded in 3D)

Units are simulation-specific (likely micrometers or scaled units based on the UG4 model parameters).

### The coordinate transform (analytically reconstructed)

The PNG rendering pipeline (`util_postanalisis_export_png.py`) uses a **fixed parallel projection camera** set at t=0:
1. Camera center: `(center_x, center_y) = ((xmin+xmax)/2, (ymin+ymax)/2)`
2. Parallel scale: `S = max(domain_width, domain_height) / 2`
3. For a square W×H = 2048×2048 image, the view spans `[center_x−S, center_x+S] × [center_y−S, center_y+S]`

**Pixel to physical transform:**
```
x_phys = center_x + (x_swc / 2048 − 0.5) × 2S
y_phys = center_y + (0.5 + y_swc / 2048) × 2S
```

**Uniform pixel scale factor:**
```
scale = 2S / 2048  [physical units per pixel, isotropic]
```

For D7.5_V3: S = 3.0, scale ≈ 0.002930 units/pixel.

### Verification (D7.5_V3, frame 033)

| Quantity | Value |
|---------|-------|
| SWC root (soma): `(x_swc, y_swc)` | (1032, −1907) |
| Computed physical: `(x_phys, y_phys)` | (0.023, −1.437) |
| Domain y-min | −1.550 |
| lsf at (0.023, −1.437, 0) via vtkProbeFilter | **−0.018** (inside cell ✓) |
| u_t at soma, t=0 | 0.249 |
| u_ca_cyt at soma, t=0 | 0.251 |

The lsf < 0 confirms the soma node is correctly placed inside the cell boundary. ✓

### Key assumptions for the transform
1. **The image was rendered with `margin_factor=0`** (confirmed by inspecting the export code: `margin_factor=0.0` is hardcoded).
2. **The camera uses t=0 bounds** (confirmed: `Camara_position` always calls `UpdatePipeline(time=t0)`).
3. **The mesh is Eulerian (fixed)** — bounds do not change across time steps. The same transform applies to all frames.
4. **The image size is 2048×2048** (hardcoded in the export pipeline and SLURM scripts).

All four assumptions are verified from the code. The transform is safe to use.

---

## 4. Sampling method

**Recommended approach**: `vtkProbeFilter` (VTK Python, available via `pvpython`)

```
vtkProbeFilter:
  SourceData  = unstructured mesh from VTKHDF reader (at requested time step)
  InputData   = vtkPolyData with SWC node positions (physical coords)
  Output      = same points with interpolated field values + vtkValidPointMask
```

- `vtkValidPointMask = 1` → point is inside the mesh domain, values are valid
- `vtkValidPointMask = 0` → point is outside the mesh — flagged in CSV, not silently dropped
- All primary fields are point data — bilinear interpolation via cell shape functions

**Why not ParaView pipeline API?**  
The `paraview.simple` API adds substantial overhead (pipeline objects, rendering state). For batch processing of 44 frames × 150 cases, using the raw VTK Python API (`import vtk`) via `pvpython` is faster and simpler.

---

## 5. Expected output dataset — CSV schema

```
case_id           — e.g. "D7p5_V3"
time_step         — integer frame index (0, 1, …, N-1)
time              — float time value from VTK (equals time_step for this dataset)
swc_frame         — zero-padded string, e.g. "033"
swc_node_id       — SWC node id (1-indexed, matches SWC file)
parent_id         — parent node id (−1 for soma root)
node_type         — SWC type: 1=soma, 3=dendrite
x_px              — pixel column (SWC x_swc)
y_px              — pixel row, negated (SWC y_swc)
x_phys            — physical x coordinate
y_phys            — physical y coordinate
z_phys            — physical z (always 0.0)
distance_from_soma — Euclidean distance from soma in physical units
path_length_from_soma — tree path length from soma in physical units
branch_order      — number of branchpoints on path from soma
is_tip            — 1 if terminal node (growth cone / dendrite tip)
u_t               — tubulin concentration
u_u               — MAP2_free concentration
u_b               — MAP2_bound concentration
u_p               — MAP2_phosphorylated concentration
u_ca_cyt          — cytosolic calcium concentration
inhibitor         — inhibitor concentration (meaningful for extracellular nodes)
lsf               — level-set value (< 0 inside cell; use as sanity check)
valid_point       — 1 if VTK probe returned a valid value, 0 if outside domain
```

---

## 6. Feasibility verdict

| Criterion | Status |
|-----------|--------|
| SWC files exist on disk | ✅ |
| VTKHDF files exist on disk | ✅ |
| Frame count matches (44 each for D7.5_V3) | ✅ |
| Field names identified (u_t, u_u, u_b, u_p, u_ca_cyt, inhibitor, lsf) | ✅ |
| Coordinate transform analytically recoverable from t=0 bounds | ✅ |
| vtkProbeFilter successful (lsf=−0.018 at soma, valid=True) | ✅ |
| pvpython available via `module load paraview/6.1.0-mesa` | ✅ |
| Scale of dataset: 151 cases × ~44 frames × ~25–100 nodes = manageable | ✅ |

**Sampling is fully feasible.** The implementation plan proceeds to `sample_fields_on_swc.py`.

---

## 7. Remaining limitations and caveats

1. **`clean_skeleton.swc` contains only topological junctions** (soma + branchpoints + tips). Branches between junctions are represented by a single edge with no intermediate nodes. If fine-grained spatial profiles are needed, use `--swc-type pixel` to include all skeleton pixels.

2. **Radius is always 1.0 pixel** — no neurite radius was estimated from the PNG. The physical radius column is therefore `1.0 × scale ≈ 0.003 units` everywhere and is not biologically meaningful.

3. **Inhibitor field is extracellular** — its value at intracellular SWC nodes is an extrapolation artifact. Use `lsf > 0` to filter if needed.

4. **Camera is set at t=0 bounds** — if the simulation domain was resized after t=0 in any case, the transform breaks. In practice, UG4 level-set simulations use a fixed Eulerian mesh, so this is not a risk.

5. **Cross-case coordinate units** — the physical coordinates differ between cases (D and V parameters affect domain extents). Path lengths and distances are in the same physical units across cases only if the simulation domain is identically sized. Verify with `x_phys` ranges in the CSV.

---

## 8. Implementation plan

| Step | Script | Status |
|------|--------|--------|
| Field sampling | `sample_fields_on_swc.py` | **Implemented** |
| Plotting | `plot_swc_field_samples.py` | **Implemented** |
| Batch processing (all cases) | `batch_sample_all_cases.sh` | Skeleton provided in script header |

**Test command (dry run, D7.5_V3, 10 nodes):**
```bash
module load paraview/6.1.0-mesa
pvpython sample_fields_on_swc.py \
    --sim-file /scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/Simulation.vtkhdf \
    --swc-dir  /scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/merged/swc \
    --out-csv  /scratch/flore0a/AnalysisResults/D7p5_V3_fields.csv \
    --case-id  D7p5_V3 \
    --dry-run
```
