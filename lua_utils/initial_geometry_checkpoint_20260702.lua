-- =====================================================================
-- initial_geometry_checkpoint_20260702.lua
--
-- Adds two things the model can call as functions (see -swcinit / -loadInitLSF
-- flags wired into Model_snew_cone_selection_3D_20260702.lua):
--
--  (1) InitialgeometrySWC(...) — build the initial level-set interface from an
--      SWC reference tree (a variable/constant-radius tube around the skeleton),
--      then redistance it to a clean signed-distance function. Mirrors the
--      contract of the built-in Initialgeometry(): returns the finished lsf.
--      This lets the SAME model be seeded from a SIM neuron OR a REAL neuron
--      (just point -swcEdges at the corresponding CSV + use the matching box).
--
--  (2) Save/LoadInitialInterface(...) — checkpoint the finished initial LSF to a
--      ConnectionViewer .vec file and load it back, so the (slow) Eikonal
--      redistancing is done ONCE and reused. The redistance is the expensive
--      part of init (150-300 pseudo-time steps); loading is ~instant.
--
--      HARD CONSTRAINT: a checkpoint is valid ONLY for the identical mesh —
--      same grid + numRefs + numPreRefs + process count. ReadVector requires
--      gridsize == num_indices(); a mismatch is undefined. Encode the mesh
--      signature in the filename (the model does this) and never mix meshes.
--
-- Requires (loaded/defined by the model BEFORE these are called):
--   * InitialgeometryStructure3d  (initial_geometry_structure_3d_20260630.lua)
--   * globals: OutflowSubsets, initNumEikonalSteps
-- =====================================================================

--- Build the initial interface from an SWC reference tree and redistance it.
-- @param ApproxSpace_lsf approximation space carrying the scalar lsf ("ca_cyt")
-- @param gridfunction_ls the LSF grid function to fill (modified in place)
-- @param reference_edges_csv path to the 3D reference-tree edges CSV
-- @return the finished (redistanced) lsf grid function
function InitialgeometrySWC(ApproxSpace_lsf, gridfunction_ls, reference_edges_csv)
    local lsf = gridfunction_ls

    -- 1. raw signed-distance tube: phi = dist_to_nearest_segment - radius
    InitialgeometryStructure3d(ApproxSpace_lsf, lsf, reference_edges_csv, "direct_phi_only")

    -- 2. redistance to a clean SDF (same Eikonal idiom the model uses at init).
    --    This is the slow step the checkpoint is designed to skip on reload.
    local steps = initNumEikonalSteps or 300
    local sdf_old = GridFunction(ApproxSpace_lsf)
    local sdf_new = GridFunction(ApproxSpace_lsf)
    local sdf_cfl = GridFunction(ApproxSpace_lsf)
    sdf_old:set(0.0); sdf_new:set(0.0); sdf_cfl:set(0.0)

    local Eik = HiResFluxBasedLSM()
    Eik:prepare_for_SDF(sdf_old, sdf_new)
    Eik:set_LSF(lsf)
    Eik:set_dt(1e-5)
    Eik:set_time_control(0.85, 0.95)
    Eik:save_CourantNumber_to(sdf_cfl)
    Eik:set_dirichlet_data(0.0)
    Eik:set_outflow_boundary(OutflowSubsets)
    Eik:set_nr_of_steps(steps)
    Eik:set_verbose(false)
    Eik:advect()
    lsf:assign(sdf_new)

    print("[initSWC] SWC interface built + redistanced (" .. steps ..
          " eikonal steps) from " .. tostring(reference_edges_csv))
    return lsf
end

--- Save the finished initial interface to a reusable checkpoint.
-- Uses UG4's restart primitive SaveToFile (one combined file). It serializes the
-- parallel storage mask, so ReadFromFile can restore a fully-usable vector --
-- unlike SaveVectorForConnectionViewer/LoadVector, which leave the loaded vector
-- as PST_UNDEFINED and break the first downstream op (VTK/Eikonal/transport).
function SaveInitialInterface(gridfunction_ls, path)
    SaveToFile(gridfunction_ls, path)
    print("[ckpt] saved initial interface -> " .. tostring(path))
end

--- Load a previously-saved initial interface, skipping recompute.
-- ReadFromFile restores BOTH the values and the parallel storage type, so the
-- loaded lsf is immediately usable by VTKOutput / Eikonal / transport.
-- Valid ONLY for the identical mesh + process count the checkpoint was saved with.
function LoadInitialInterface(gridfunction_ls, path)
    ReadFromFile(gridfunction_ls, path)
    print("[ckpt] loaded initial interface <- " .. tostring(path) ..
          "  (skipped Eikonal redistancing)")
end
