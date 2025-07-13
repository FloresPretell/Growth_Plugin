

--- Gathers and returns all simulation parameters
-- Retrieves each parameter from the command line via util.GetParamNumber,
-- falling back to predefined defaults when no flag is provided.
-- Parameters are structured into three sections:
--   * Intracellular: tubulin, MAP proteins, and calcium
--   * Extracellular: inhibitor production rate, capacity, diffusion, and degradation
--   * LevelSet: elongation rate for the level-set method
--
-- @usage
-- ug_load_script("plugins/NeuroGrowth/lua_utils/test_parameters.lua")
-- local params = GetParams()
-- print(params.Intracellular.tubulin.diffusion) (obtain the number value)
--
-- @return table params
--   A nested table containing all parameters ready for use in the simulation  
function GetParams()

    -- #region Define all the parameters in the code
    ---@class TubulinParams
    ---@field diffusion    number "Diffusion coefficient of tubulin: 5 × (10 μm)^2/(500 s) = 1 μm²/s"
    ---@field velocity     number "Transport coefficient (direction × intensity):  5 × (10 μm)/(500 s) = 0.1 μm/s"
    ---@field somaFlux     number "Tubulin flux at soma boundary: 0.25 × 40 μM = 10 μM"
    ---@field assemblyRate number "Tubulin assembly rate: 50 x 1 /(500 s · 6.25 μM) = 0.016 μM⁻¹·s⁻¹ " 
    ---@field assemblyRateCurvatureminimun   table<number, table<number, number>>   "Minimum curvature threshold for assembly rate, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field assemblyRateCurvaturemaximun   table<number, table<number, number>>   "Maximum curvature threshold for assembly rate, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"


    ---@class MapParams
    ---@field diffusion number "Diffusion coefficient of MAP proteins: 1 × (10 μm)^2/(500 s) = 0.2 μm²/s"
    ---@field source    number "Source coefficient (direction × intensity): -0.5 × (6.25 μM)/(500 s) = -0.00625 μM/s"
    ---@field K1 number "Reaction constant K1: 1 x 1/500 s = 0.002 s⁻¹"
    ---@field K2 number "Reaction constant K2: 0.001 x 1/500 s = 2e-6 s⁻¹ "
    ---@field K3 number "Reaction constant K3: 1 x 1/500 s = 0.002 s⁻¹"
    ---@field K4 number "Reaction constant K4: 0.001x 1/500 s = 2e-6 s⁻¹"

    ---@class CalciumParams
    ---@field initial number "Initial calcium concentration: 0.25 x 0.5uM = 0.125 uM"
    ---@field difussion number "Calcium diffusion coefficient: 1 × (10 μm)^2/(500 s) = 0.2 μm²/s"

    ---@class IntracellularParams
    ---@field tubulin  TubulinParams
    ---@field map      MapParams
    ---@field calcium  CalciumParams


    ---@class InhibitionParams
    ---@field rateSource  number "Inhibitor production rate coefficient: 5 x 1/500 s = 0.01 s⁻¹"
    ---@field capacity    number "Inhibitor capacity: 0.01 × 1.35 μM = 0.0135 μM"
    ---@field diffusion   number "Inhibitor diffusion coefficient: 1 × (10 μm)^2/(500 s) = 0.2 μm²/s"
    ---@field degradation number "Inhibitor degradation rate: 10 x 1/500 s = 0.02 s⁻¹"

    ---@class ExtracellularParams
    ---@field inhibition InhibitionParams


    ---@class elongationParams
    ---@field rate number "Elongation rate: 225 × (10 μm)/(500 s · 40 μM · 6.25 μM) = 0.018 μm/(s·μM²)"

    ---@class LevelSetParams
    ---@field elongation  elongationParams

    ---@class Params
    ---@field Intracellular IntracellularParams
    ---@field Extracellular ExtracellularParams
    ---@field LevelSet LevelSetParams

    ---@type Params
    -- #endregion
    local params = {
    -- Intracellular parameters
    Intracellular = {
    tubulin = {
        diffusion    = util.GetParamNumber("-tubulinDiffusion",      5,   "Diffusion coefficient of tubulin"),   -- 5 × (10 μm)^2/(500 s) = 1 μm²/s
        velocity     = util.GetParamNumber("-tubulinVelocity",       5,   "Transport coefficient (direction × intensity)"), -- 5 × (10 μm)/(500 s) = 0.1 μm/s
        somaFlux     = util.GetParamNumber("-somaTubulinBoundaryFlux", 0.25, "Tubulin flux at soma boundary"), -- 0.25 × 40 μM = 10 μM
        assemblyRate = util.GetParamNumber("-tubulinAssemblyRate",   50,  "Tubulin assembly rate"), -- 50 x 1 /(500 s · 6.25 μM) = 0.016 μM⁻¹·s⁻¹
        assemblyRateCurvatureminimun = {
            [2] = {  -- dim = 2
                [5] = util.GetParamNumber("-assamblyCurvaturemin",  5.5, "Minimun Curvature for assambly rate in 2D with refinemnte 5 or less"),
                [6] = util.GetParamNumber("-assamblyCurvaturemin",  10, "Minimun Curvature for assambly rate in 2D with refinemnte 6 or more"),
            },
            [3] = {  -- dim = 3
                [5] = util.GetParamNumber("-assamblyCurvaturemin",  20, "Minimun Curvature for assambly rate in 3D without defining refinment"),
            },
        },
        assemblyRateCurvaturemaximun = {
            [2] = {  -- dim = 2
                [5] = util.GetParamNumber("-assamblyCurvaturemin",  18.5, "Minimun Curvature for assambly rate in 2D with refinemnte 5 or less"),
                [6] = util.GetParamNumber("-assamblyCurvaturemin",  20, "Minimun Curvature for assambly rate in 2D with refinemnte 6 or more"),
            },
            [3] = {  -- dim = 3
                [5] = util.GetParamNumber("-assamblyCurvaturemin",  27, "Minimun Curvature for assambly rate in 3D without defining refinment"),
            },
        },  
    },
    map = {
        diffusion = util.GetParamNumber("-mapDiffusion", 1,    "Diffusion coefficient of MAP proteins"), -- 1 × (10 μm)^2/(500 s) = 0.2 μm²/s
        source    = util.GetParamNumber("-mapSource",    -0.5, "Source coefficient (direction × intensity)"), -- -0.5 × (6.25 μM)/(500 s) = -0.00625 μM/s
        K1 = util.GetParamNumber("-reactionRateK1", 1,     "Reaction constant K1"), -- 1 x 1/500 s = 0.002 s⁻¹
        K2 = util.GetParamNumber("-reactionRateK2", 0.001, "Reaction constant K2"), -- 0.001 x 1/500 s = 2e-6 s⁻¹
        K3 = util.GetParamNumber("-reactionRateK3", 1,     "Reaction constant K3"),  -- 1 x 1/500 s = 0.002 s⁻¹
        K4 = util.GetParamNumber("-reactionRateK4", 0.001, "Reaction constant K4"), -- 0.001x 1/500 s = 2e-6 s⁻¹
    },
    calcium = {
        initial = util.GetParamNumber("-initialCalcium", 0.25, "Initial calcium concentration"), -- 0.25 x 0.5uM = 0.125 uM
        difussion = util.GetParamNumber("-calciumDiffusion", 1, "Calcium diffusion coefficient"), -- 1 × (10 μm)^2/(500 s) = 0.2 μm²/s
    },
    },

    -- Extracellular parameters
    Extracellular = {
    inhibition = {
        rateSource  = util.GetParamNumber("-inhibitionRateSource",  5,   "Inhibitor production rate coefficient"), -- 5 x 1/500 s = 0.01 s⁻¹
        capacity    = util.GetParamNumber("-inhibitionCapacity",    0.01,"Inhibitor capacity"), -- 0.01 × 1.35 μM = 0.0135 μM
        diffusion   = util.GetParamNumber("-inhibitionDiffusion",   1,   "Inhibitor diffusion coefficient"),  -- 1 × (10 μm)^2/(500 s) = 0.2 μm²/s
        degradation = util.GetParamNumber("-inhibitionDegradation",10,   "Inhibitor degradation rate"), -- 10 x 1/500 s = 0.02 s⁻¹
    },
    },

    -- Interface parameters for level set method
    LevelSet = {
    -- Velocity parameters
    elongation = {
        rate = util.GetParamNumber("-elongationRate", 225, "Elongation rate"), -- 225 × (10 μm)/(500 s · 40 μM · 6.25 μM) = 0.018 μm/(s·μM²)
    },
    },
    }

    return params
end

