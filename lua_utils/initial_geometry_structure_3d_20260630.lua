-- =====================================================================
-- initial_geometry_structure_3d_20260630.lua
--
-- 3D port of Growth_Plugin/lua_utils/initial_geometry_structure.lua.
-- Builds the level-set phi as a signed distance to a 3D reference tree
-- (an SWC skeleton given as edges), i.e. a variable-radius tube:
--
--     phi(x,y,z) = dist_to_nearest_segment(x,y,z) - radius_of_that_segment
--
-- phi < 0 inside the tube (the neuron), phi > 0 outside. Interpolated onto
-- the LSF grid function via Interpolate(LuaUserNumber(...)).
--
-- CSV columns (from swc_to_edges3d.py):
--   edge_id, x0, y0, z0, x1, y1, z1, radius
-- =====================================================================

local function file_exists(path)
    local f = io.open(path, "r")
    if f then f:close(); return true end
    return false
end

local function split_csv_line(line)
    local fields = {}
    for field in string.gmatch(line, "([^,]+)") do
        fields[#fields + 1] = field:gsub("^%s+", ""):gsub("%s+$", "")
    end
    return fields
end

local function load_reference_tree_edges_3d(csv_path)
    local f, err = io.open(csv_path, "r")
    if not f then
        error("Cannot open 3D reference edge CSV: " .. tostring(csv_path) .. ": " .. tostring(err), 0)
    end
    local header = f:read("*l")
    if header == nil then
        error("Empty 3D reference edge CSV: " .. tostring(csv_path), 0)
    end
    -- column index lookup
    local cols = split_csv_line(header)
    local idx = {}
    for i, name in ipairs(cols) do idx[name] = i end
    for _, req in ipairs({"x0","y0","z0","x1","y1","z1","radius"}) do
        if idx[req] == nil then
            error("3D reference CSV missing column '" .. req .. "': " .. tostring(csv_path), 0)
        end
    end

    local edges = {}
    for line in f:lines() do
        if line ~= "" then
            local fld = split_csv_line(line)
            edges[#edges + 1] = {
                x0 = tonumber(fld[idx.x0]), y0 = tonumber(fld[idx.y0]), z0 = tonumber(fld[idx.z0]),
                x1 = tonumber(fld[idx.x1]), y1 = tonumber(fld[idx.y1]), z1 = tonumber(fld[idx.z1]),
                radius = tonumber(fld[idx.radius]),
            }
        end
    end
    f:close()
    if #edges == 0 then
        error("No edges loaded from 3D reference CSV: " .. tostring(csv_path), 0)
    end
    return edges
end

-- squared distance from point p to segment a-b, plus the parametric t
local function dist2_point_segment_3d(px,py,pz, ax,ay,az, bx,by,bz)
    local dx, dy, dz = bx-ax, by-ay, bz-az
    local len2 = dx*dx + dy*dy + dz*dz
    local t = 0.0
    if len2 > 1e-30 then
        t = ((px-ax)*dx + (py-ay)*dy + (pz-az)*dz) / len2
        if t < 0.0 then t = 0.0 elseif t > 1.0 then t = 1.0 end
    end
    local cx, cy, cz = ax + t*dx, ay + t*dy, az + t*dz
    local ex, ey, ez = px-cx, py-cy, pz-cz
    return ex*ex + ey*ey + ez*ez
end

-- module-level state (set by InitialgeometryStructure3d)
local STRUCTURE_REFERENCE_TREE_EDGES_3D = {}

-- returns phi = dist_to_nearest_segment - radius_of_that_segment
function StructureReferenceTreePhi3d(x, y, z, t)
    local edges = STRUCTURE_REFERENCE_TREE_EDGES_3D
    local best_d2 = math.huge
    local best_radius = 0.0
    for i = 1, #edges do
        local e = edges[i]
        local d2 = dist2_point_segment_3d(x, y, z, e.x0, e.y0, e.z0, e.x1, e.y1, e.z1)
        if d2 < best_d2 then
            best_d2 = d2
            best_radius = e.radius
        end
    end
    return math.sqrt(best_d2) - best_radius
end

-- Main entry: initialise the LSF grid function from a 3D reference-edge CSV.
--   ApproxSpace_lsf : approximation space carrying the lsf component "ca_cyt"
--   gridfunction_ls : the LSF grid function to fill
--   reference_edges_csv : path to the 3D edges CSV
--   mode : "direct_phi_only" (default)
function InitialgeometryStructure3d(ApproxSpace_lsf, gridfunction_ls, reference_edges_csv, mode)
    if reference_edges_csv == nil or reference_edges_csv == "" then
        error("InitialgeometryStructure3d requires -referenceTreeEdgesCSV", 0)
    end
    if not file_exists(reference_edges_csv) then
        error("InitialgeometryStructure3d cannot find reference edge CSV: " .. tostring(reference_edges_csv), 0)
    end
    mode = mode or "direct_phi_only"

    STRUCTURE_REFERENCE_TREE_EDGES_3D = load_reference_tree_edges_3d(reference_edges_csv)
    print("[init3d] loaded " .. #STRUCTURE_REFERENCE_TREE_EDGES_3D ..
          " reference edges from " .. reference_edges_csv)

    local lsf = gridfunction_ls
    lsf:set(0.0)
    Interpolate(LuaUserNumber("StructureReferenceTreePhi3d"), lsf, "ca_cyt", 0.0)
    print("[init3d] phi interpolated (signed-distance tube around the SWC tree)")

    if mode ~= "direct_phi_only" then
        print("[init3d] mode '" .. tostring(mode) .. "' not enabled; using direct_phi_only.")
    end
    return lsf
end
