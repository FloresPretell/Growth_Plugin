

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
    ---@field initial number "Initial tubulin concentration: 0.3 adim "
    ---@field diffusion    number "Diffusion coefficient of tubulin: 5 × (10 μm)^2/(500 s) = 1 μm²/s"
    ---@field velocity     number "Transport coefficient (direction × intensity):  5 × (10 μm)/(500 s) = 0.1 μm/s"
    ---@field somaFlux     number "Tubulin flux at soma boundary: 0.25 × 40 μM = 10 μM"
    ---@field assemblyRate number "Tubulin assembly rate: 50 x 1 /(500 s · 6.25 μM) = 0.016 μM⁻¹·s⁻¹ " 
    ---@field assemblyRateCurvatureminimun   table<number, table<number, number>>   "Minimum curvature threshold for assembly rate, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field assemblyRateCurvaturemaximun   table<number, table<number, number>>   "Maximum curvature threshold for assembly rate, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"


    ---@class MapParams
    ---@field initial   number "Initial MAP protein concentration: 0 adim"
    ---@field diffusion number "Diffusion coefficient of MAP proteins: 1 × (10 μm)^2/(500 s) = 0.2 μm²/s"
    ---@field source    number "Source coefficient (direction × intensity): -0.5 × (6.25 μM)/(500 s) = -0.00625 μM/s"
    ---@field K1 number "Reaction constant K1: 1 x 1/500 s = 0.002 s⁻¹"
    ---@field K2 number "Reaction constant K2: 0.001 x 1/500 s = 2e-6 s⁻¹ "
    ---@field K3 number "Reaction constant K3: 1 x 1/500 s = 0.002 s⁻¹"
    ---@field K4 number "Reaction constant K4: 0.001x 1/500 s = 2e-6 s⁻¹"
    ---@field phosphorylationBCurvatureMinimun table<number, table<number, number>>   "Minimum curvature threshold for Phosphorilation of B (remove B), keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field phosphorylationBCurvatureMaximun table<number, table<number, number>>   "Maximum curvature threshold for Phosphorilation of B (remove B), keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field despholymerizationBCurvatureMinimun table<number, table<number, number>>   "Minimum curvature threshold for Phosphorilation of B (remove B), keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field despholymerizationBCurvatureMaximun table<number, table<number, number>>   "Maximum curvature threshold for Phosphorilation of B (remove B), keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field synthesisUCurvatureMinimun table<number, table<number, number>>   "Minimum curvature threshold for syntesis MAP2 free (U) on growth cones, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field synthesisUCurvatureMaximun table<number, table<number, number>>   "Maximum curvature threshold for syntesis MAP2 free (U) on growth cones, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field phosphorylationPCurvatureMinimun table<number, table<number, number>>   "Minimum curvature threshold for Phosphorilation of B (receive P), keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field phosphorylationPCurvatureMaximun table<number, table<number, number>>   "Maximum curvature threshold for Phosphorilation of B (receive P), keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
   

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
    ---@field synthesisInhCurvatureMinimun table<number, table<number, number>>   "Minimum curvature threshold for syntesis MAP2 free (U) on growth cones, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field synthesisInhCurvatureMaximun table<number, table<number, number>>   "Maximum curvature threshold for syntesis MAP2 free (U) on growth cones, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
  

    ---@class ExtracellularParams
    ---@field inhibition InhibitionParams


    ---@class elongationParams
    ---@field rate number "Elongation rate: 225 × (10 μm)/(500 s · 40 μM · 6.25 μM) = 0.018 μm/(s·μM²)"
    ---@field VelocityCurvatureMinimun table<number, table<number, number>>   "Minimum curvature threshold for cone area where it will deform , keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field VelocityCurvatureMaximun table<number, table<number, number>>   "Maximum curvature threshold for cone area where it will deform, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field inhibitionMinConcentrationAvoidance number "Concentration of Inhibitor to start the avoidance process, value less than this number is elongation"
    ---@field inhibitionMinConcentrationRetraction number "Concentration of Inhibitor to start the retraction process"
    ---@field tubulineConcentrationThreshold number "Concentration of Tubulin to start the growth process; and if the value is less, the interface wont deform"

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
                initial = util.GetParamNumber("-initialTubulin", 0.3, "Initial tubulin concentration"),
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
                initial = util.GetParamNumber("-initialMAP", 0, "Initial MAP protein concentration"),
                diffusion = util.GetParamNumber("-mapDiffusion", 1,    "Diffusion coefficient of MAP proteins"), -- 1 × (10 μm)^2/(500 s) = 0.2 μm²/s
                source    = util.GetParamNumber("-mapSource",    -0.5, "Source coefficient (direction × intensity)"), -- -0.5 × (6.25 μM)/(500 s) = -0.00625 μM/s
                K1 = util.GetParamNumber("-reactionRateK1", 1,     "Reaction constant K1"), -- 1 x 1/500 s = 0.002 s⁻¹
                K2 = util.GetParamNumber("-reactionRateK2", 0.001, "Reaction constant K2"), -- 0.001 x 1/500 s = 2e-6 s⁻¹
                K3 = util.GetParamNumber("-reactionRateK3", 1,     "Reaction constant K3"),  -- 1 x 1/500 s = 0.002 s⁻¹
                K4 = util.GetParamNumber("-reactionRateK4", 0.001, "Reaction constant K4"), -- 0.001x 1/500 s = 2e-6 s⁻¹

                phosphorylationBCurvatureMinimun = {
                    [2] = {
                        [5] = util.GetParamNumber("-PhosmapBCurvaturemin",  3,  "Min interval (2D, ref≤5)"),
                        [6] = util.GetParamNumber("-PhosmapBCurvaturemin", 10,  "Min interval (2D, ref≥6)"),
                    },
                    [3] = {
                        [5] = util.GetParamNumber("-PhosmapBCurvaturemin", 20,  "Min interval (3D, use ref≤5)"),
                    },
                    },
                phosphorylationBCurvatureMaximun = {
                    [2] = {
                        [5] = util.GetParamNumber("-PhosmapBCurvaturmax", 18.5, "Max interval (2D, ref≤5)"),
                        [6] = util.GetParamNumber("-PhosmapBCurvaturmax", 20,   "Max interval (2D, ref≥6)"),
                    },
                    [3] = {
                        [5] = util.GetParamNumber("-PhosmapBCurvaturmax", 27,   "Max interval (3D, use ref≤5)"),
                    },
                },


                despholymerizationBCurvatureMinimun = {
                    [2] = {
                        [5] = util.GetParamNumber("-DesPolymapBCurvaturemin",  5.5,  "Min interval (2D, ref≤5)"),
                        [6] = util.GetParamNumber("-DesPolymapBCurvaturemin", 10,  "Min interval (2D, ref≥6)"),
                    },
                    [3] = {
                        [5] = util.GetParamNumber("-DesPolymapBCurvaturemin", 20,  "Min interval (3D, use ref≤5)"),
                    },
                    },
                despholymerizationBCurvatureMaximun = {
                    [2] = {
                        [5] = util.GetParamNumber("-DesPolymapBCurvaturemax", 18.5, "Max interval (2D, ref≤5)"),
                        [6] = util.GetParamNumber("-DesPolymapBCurvaturemax", 20,   "Max interval (2D, ref≥6)"),
                    },
                    [3] = {
                        [5] = util.GetParamNumber("-DesPolymapBCurvaturemax", 27,   "Max interval (3D, use ref≤5)"),
                    },
                },



                synthesisUCurvatureMinimun = {
                    [2] = {  -- dim=2
                        [5] = util.GetParamNumber("-synthUCurvaturemin",  4,"Min interval for U Synthesis (2D, numRefs ≤ 5)"),
                        [6] = util.GetParamNumber("-synthUCurvaturemin", 10,"Min interval for U Synthesis (2D, numRefs ≥ 6)"),
                    },
                    [3] = {  -- dim=3
                        [5] = util.GetParamNumber("-synthUCurvaturemin", 20,"Min interval for U Synthesis (3D, use numRefs ≤ 5)"),
                    },
                },
                synthesisUCurvatureMaximun = {
                    [2] = {
                        [5] = util.GetParamNumber("-synthUCurvaturemax", 18.5,"Max interval for U Synthesis (2D, numRefs ≤ 5)"),
                        [6] = util.GetParamNumber("-synthUCurvaturemax", 20,"Max interval for U Synthesis (2D, numRefs ≥ 6)"),
                    },
                    [3] = {
                        [5] = util.GetParamNumber("-synthUCurvaturemax", 27,"Max interval for U Synthesis (3D, use numRefs ≤ 5)"),
                    },
                },

                phosphorylationPCurvatureMinimun = {
                    [2] = {
                        [5] = util.GetParamNumber("-PhosmapPCurvaturemin",  5.5,  "Min interval (2D, ref≤5)"),
                        [6] = util.GetParamNumber("-PhosmapPCurvaturemin", 10,  "Min interval (2D, ref≥6)"),
                    },
                    [3] = {
                        [5] = util.GetParamNumber("-PhosmapPCurvaturemin", 20,  "Min interval (3D, use ref≤5)"),
                    },
                    },
                phosphorylationPCurvatureMaximun = {
                    [2] = {
                        [5] = util.GetParamNumber("-PhosmapPCurvaturmax", 18.5, "Max interval (2D, ref≤5)"),
                        [6] = util.GetParamNumber("-PhosmapPCurvaturmax", 20,   "Max interval (2D, ref≥6)"),
                    },
                    [3] = {
                        [5] = util.GetParamNumber("-PhosmapPCurvaturmax", 27,   "Max interval (3D, use ref≤5)"),
                    },
                },


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

                synthesisInhCurvatureMinimun = {
                    [2] = {  -- dim=2
                        [5] = util.GetParamNumber("-synthICurvaturemin",  4,"Min interval for Inhibitor synthesis  (2D, numRefs ≤ 5)"),
                        [6] = util.GetParamNumber("-synthICurvaturemin", 3,"Min interval for Inhibitor synthesis (2D, numRefs ≥ 6)"),
                    },
                    [3] = {  -- dim=3
                        [5] = util.GetParamNumber("-synthICurvaturemin", 20,"Min interval for Inhibitor synthesis (3D, use numRefs ≤ 5)"),
                    },
                },
                synthesisInhCurvatureMaximun = {
                    [2] = {
                        [5] = util.GetParamNumber("-synthICurvaturemax", 25,"Max interval for Inhibitor synthesis (2D, numRefs ≤ 5)"),
                        [6] = util.GetParamNumber("-synthICurvaturemax", 20,"Max interval for Inhibitor synthesis (2D, numRefs ≥ 6)"),
                    },
                    [3] = {
                        [5] = util.GetParamNumber("-synthICurvaturemax", 27,"Max interval for Inhibitor synthesis (3D, use numRefs ≤ 5)"),
                    },
                },

            },
        },

    -- Interface parameters for level set method
        LevelSet = {

        -- Velocity parameters
            elongation = {
                rate = util.GetParamNumber("-elongationRate", 225, "Elongation rate"), -- 225 × (10 μm)/(500 s · 40 μM · 6.25 μM) = 0.018 μm/(s·μM²)

                VelocityCurvatureMinimun = {
                    [2] = {
                        [5] = util.GetParamNumber("-VelCurvaturemin",  7,  "Min interval (2D, ref≤5)"),
                        [6] = util.GetParamNumber("-VelCurvaturemin", 10.5,  "Min interval (2D, ref≥6)"),
                    },
                    [3] = {
                        [5] = util.GetParamNumber("-VelCurvaturemin", 20,  "Min interval (3D, use ref≤5)"),
                    },
                    },
                VelocityCurvatureMaximun = {
                    [2] = {
                        [5] = util.GetParamNumber("-VelCurvaturemax", 18, "Max interval (2D, ref≤5)"),
                        [6] = util.GetParamNumber("-VelCurvaturemax", 20,   "Max interval (2D, ref≥6)"),
                    },
                    [3] = {
                        [5] = util.GetParamNumber("-VelCurvaturemax", 27,   "Max interval (3D, use ref≤5)"),
                    },
                },

                inhibitionMinConcentrationAvoidance = util.GetParamNumber("-inhibitionAvoidance", 0.0062, "Concentration of Inhibitor to start the avoidance process, value less than this number is elongation"),
                inhibitionMinConcentrationRetraction = util.GetParamNumber("-inhibitionRetraction", 0.0074, "Concentration of Inhibitor to start the retraction process"),

                tubulineConcentrationThreshold = util.GetParamNumber("-tubulineGrowthThreshold", 0.4, "Concentration of Tubulin to start the growth process, value less than this number is retraction"),
            },
        },
    }

    return params
end

