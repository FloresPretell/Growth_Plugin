-- Analytical reference-tree level-set initialization for production use.

local function file_exists(path)
    local f = io.open(path, "r")
    if f ~= nil then
        f:close()
        return true
    end
    return false
end

local function split_csv_line(line)
    local fields = {}
    for field in string.gmatch(line .. ",", "([^,]*),") do
        fields[#fields + 1] = (field:gsub("^%s+", ""):gsub("%s+$", ""))
    end
    return fields
end

local function load_reference_tree_edges(csv_path)
    local f, err = io.open(csv_path, "r")
    if f == nil then
        error("Cannot open reference edge CSV: " .. tostring(csv_path) .. ": " .. tostring(err), 0)
    end

    local header_line = f:read("*l")
    if header_line == nil then
        f:close()
        error("Empty reference edge CSV: " .. tostring(csv_path), 0)
    end
    local header = split_csv_line(header_line)
    local idx = {}
    for i = 1, #header do
        idx[header[i]] = i
    end

    local edge_id_col = idx.edge_id or idx.ref_edge_id
    if edge_id_col == nil then
        f:close()
        error("Reference edge CSV is missing edge_id/ref_edge_id: " .. tostring(csv_path), 0)
    end

    local edges = {}
    for line in f:lines() do
        if line ~= nil and line ~= "" then
            local fields = split_csv_line(line)
            edges[#edges + 1] = {
                edge_id = tonumber(fields[edge_id_col]),
                ref_vertex_id_0 = tonumber(fields[idx.ref_vertex_id_0]),
                ref_vertex_id_1 = tonumber(fields[idx.ref_vertex_id_1]),
                x0 = tonumber(fields[idx.x0]),
                y0 = tonumber(fields[idx.y0]),
                z0 = tonumber(fields[idx.z0]),
                x1 = tonumber(fields[idx.x1]),
                y1 = tonumber(fields[idx.y1]),
                z1 = tonumber(fields[idx.z1]),
                length = tonumber(fields[idx.length])
            }
        end
    end
    f:close()
    return edges
end

local function closest_point_on_segment_2d(x, y, ax, ay, bx, by)
    local vx = bx - ax
    local vy = by - ay
    local wx = x - ax
    local wy = y - ay
    local denom = vx * vx + vy * vy
    local t = 0.0
    if denom > 0.0 then
        t = (wx * vx + wy * vy) / denom
        if t < 0.0 then
            t = 0.0
        elseif t > 1.0 then
            t = 1.0
        end
    end

    local closest_x = ax + t * vx
    local closest_y = ay + t * vy
    local dx = closest_x - x
    local dy = closest_y - y
    local distance = math.sqrt(dx * dx + dy * dy)
    return closest_x, closest_y, t, distance
end

local function distance_to_reference_tree_2d(x, y, edges)
    local best_distance = math.huge
    local best_edge_id = -1
    local best_x = 0.0
    local best_y = 0.0

    for i = 1, #edges do
        local edge = edges[i]
        local closest_x, closest_y, _, distance = closest_point_on_segment_2d(
            x,
            y,
            edge.x0,
            edge.y0,
            edge.x1,
            edge.y1
        )
        if distance < best_distance then
            best_distance = distance
            best_edge_id = edge.edge_id
            best_x = closest_x
            best_y = closest_y
        end
    end

    return best_distance, best_edge_id, best_x, best_y, best_x - x, best_y - y
end

local STRUCTURE_REFERENCE_TREE_EDGES = {}
local STRUCTURE_REFERENCE_TREE_RADIUS = 0.0

function StructureReferenceTreePhi2d(x, y, t)
    local d = distance_to_reference_tree_2d(x, y, STRUCTURE_REFERENCE_TREE_EDGES)
    return d - STRUCTURE_REFERENCE_TREE_RADIUS
end

function InitialgeometryStructure(ApproxSpace_lsf, gridfunction_ls, reference_edges_csv, radius, mode)
    if reference_edges_csv == nil or reference_edges_csv == "" then
        error("InitialgeometryStructure requires -referenceTreeEdgesCSV", 0)
    end
    if not file_exists(reference_edges_csv) then
        error("InitialgeometryStructure cannot find reference edge CSV: " .. tostring(reference_edges_csv), 0)
    end
    if radius == nil or radius <= 0.0 then
        error("InitialgeometryStructure requires a positive radius", 0)
    end

    mode = mode or "direct_phi_only"
    STRUCTURE_REFERENCE_TREE_EDGES = load_reference_tree_edges(reference_edges_csv)
    STRUCTURE_REFERENCE_TREE_RADIUS = radius

    local lsf = gridfunction_ls
    lsf:set(0.0)
    Interpolate(LuaUserNumber("StructureReferenceTreePhi2d"), lsf, "ca_cyt", 0.0)

    if mode == "direct_phi_only" then
        return lsf
    end

    print("InitialgeometryStructure: mode '" .. tostring(mode) .. "' is not enabled for production yet; returning direct_phi_only result.")
    return lsf
end
