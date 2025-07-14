ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")


--- Function that Gathers and returns all simulation parameters
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


    ---@class curvature2DParams
    ---@field isotropic number "Isotropic curvature for 2D level set method: 0.007"
    ---@field anisotropic number "Anisotropic curvature for 2
    ---@class curvature3DParams
    ---@field isotropicFirst number "Isotropic curvature for 3D level set method: 0.00012"
    ---@field anisotropic number "Anisotropic curvature for 3D level set method: 0.11"
    ---@field isotropicSecond number "Isotropic curvature for 3D level set method (second time): 0.0001"

    ---@class smoothParams
    ---@field Transporteq table<string, number> "Diffusion term of the Implicit discretization of the transport equation, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field Curvature2D curvature2DParams 
    ---@field Curvature3D curvature3DParams
    ---@field InterfaceVelocity table<number, number> "Interface velocity for level set method, keyed by [ref_def],  where ref_def is 4 (numRefs ≤ 4) or 5 (numRefs ≥ 5)"
    ---@field ExtendedVelocity table<number, number> "Extended velocity for level set method, keyed by [dim],  where dim is 2 (2D) or 3 (3D)"

    ---@class elongationParams
    ---@field rate number "Elongation rate: 225 × (10 μm)/(500 s · 40 μM · 6.25 μM) = 0.018 μm/(s·μM²)"
    ---@field VelocityCurvatureMinimun table<number, table<number, number>>   "Minimum curvature threshold for cone area where it will deform , keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field VelocityCurvatureMaximun table<number, table<number, number>>   "Maximum curvature threshold for cone area where it will deform, keyed by [dim][ref_def],  where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6)"
    ---@field inhibitionMinConcentrationAvoidance number "Concentration of Inhibitor to start the avoidance process, value less than this number is elongation"
    ---@field inhibitionMinConcentrationRetraction number "Concentration of Inhibitor to start the retraction process"
    ---@field tubulineConcentrationThreshold number "Concentration of Tubulin to start the growth process; and if the value is less, the interface wont deform"

    ---@class LevelSetParams
    ---@field elongation  elongationParams
    ---@field smooth smoothParams

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
                inhibitionMinConcentrationAvoidance = util.GetParamNumber("-inhibitionAvoidance", 0.0062, "Concentration of Inhibitor to start the avoidance process, value less than this number is elongation"),
                inhibitionMinConcentrationRetraction = util.GetParamNumber("-inhibitionRetraction", 0.0074, "Concentration of Inhibitor to start the retraction process"),
                tubulineConcentrationThreshold = util.GetParamNumber("-tubulineGrowthThreshold", 0.4, "Concentration of Tubulin to start the growth process, value less than this number is retraction"),
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
            },

        -- Smooth parameters 
            smooth = {
                Transporteq = {
                    diffusion = {
                        [2] = {
                            [5] = util.GetParamNumber("-synthICurvaturemax", 0.015,"Diffusion term of the Implicit discretization of the transport equation (2D, numRefs ≤ 5)"),
                            [6] = util.GetParamNumber("-synthICurvaturemax", 0.00001,"Diffusion term of the Implicit discretization of the transport equation (2D, numRefs ≥ 6)"),
                        },
                        [3] = {
                            [5] = util.GetParamNumber("-synthICurvaturemax", 0.01, "Diffusion term of the Implicit discretization of the transport equation (3D, use numRefs ≤ 5)"),
                        },
                    },
                },
                Curvature2D = {
                    isotropic = util.GetParamNumber("-isotropicCurvature2D", 0.007, "Isotropic curvature for 2D level set method"),
                    anisotropic = util.GetParamNumber("-anisotropicCurvature2D", 0.0825, "Anisotropic curvature for 2D level set method"),
                    
                },
                Curvature3D = {
                    isotropicFirst = util.GetParamNumber("-isotropicCurvature3D", 0.00012, "Isotropic curvature for 3D level set method"),
                    anisotropic = util.GetParamNumber("-anisotropicCurvature3D", 0.11, "Anisotropic curvature for 3D level set method"),
                    isotropicSecond = util.GetParamNumber("-isotropicCurvature3DSecond", 0.0001, "Isotropic curvature for 3D level set method (second time)"),
                },
                InterfaceVelocity = {
                    [4] = util.GetParamNumber("-interfaceVelocity", 0.001, "Interface velocity for level set method (refinemente 4 or less)"),
                    [5] = util.GetParamNumber("-interfaceVelocity", 0.00039, "Interface velocity for level set method (refinemente 5 or more)"),  -- SMALl vakue because the multiplication B and T is small when the inhibition is difussed
                },

                ExtendedVelocity = {
                    [2] = util.GetParamNumber("-extendedVelocity2D", 0.005, "Extended velocity for 2D level set method"), -- In 2D 
                    [3] = util.GetParamNumber("-extendedVelocity3D", 0.005, "Extended velocity for 3D level set method"), --- In 3D
                }

            },
        },
    }

    return params
end


---@class GridFunctionGroup
---@field u GridFunction -- Main Gridfunction that have the values of the intracellular simulation (levelset < 0)
---@field u2 GridFunction -- Main Gridfunction that have the values of the extracellular simulation (levelset > 0)
---@field lsf GridFunction -- Gridfunction that define the values of the levelset function. This will define the interface between the intracellular and extracellular simulation. This lsf will be the input of the LSGFsimpleDomainDiscretization() for the intracellular model (negative value is the domain of the intracellular simulation).
---@field lsf_oposite GridFunction -- Gridfunction that define the oposite sign of the valules of the levelset function . This have the same function of define the interface between the intracellular and extracellular simulation. This lsf_oposite will be the input of the LSGFsimpleDomainDiscretization() for the extracellular model (negative value is the domain of the extracellular simulation).
---@field sdf_new GridFunction -- New values of the level set function after the reinitialize it (Used as update in EikonalDisc = HiResFluxBasedLSM() )
---@field sdf_old GridFunction -- Used for solving the sdf_new
---@field sdf_cfl_num_gf GridFunction -- Used to store the signed distance field for CFL number during sdf_new computation (EikonalDisc = HiResFluxBasedLSM())
---@field u_new GridFunction -- After Extending the interface velocity (scalar) (ExtensionDisc = HiResFluxBasedLSM()), this is the result of the velocity that will be used to deform the interface 
---@field u_old GridFunction -- Used for solving the u_new
---@field norm_vel_int GridFunction -- Normal interface velocity obtained using compute_normal_vel(), where obtain the normal velocity of the velocity field (calculated from LSConcentrationDepentVelocity()) 
---@field curvature GridFunction -- Curvature of the level set function, it is computed using  : Kapla_Func = Kappla_LS(); Kapla_Func:compute_kappla(lsf, curvature)
---@field Eikonal_Inside GridFunction -- Values of the Stationary Poisson Equation inside the neuron (were the source term is the growth cones). It would be used to obtain the gradient and them used as: the direction of the tubulin velocity (convection) or neuronal growth cones (velocity fied in linker LSConcentrationDepentVelocity() )
---@field Eikonal_Inside_direc_inter GridFunction -- Values of the Eikanol Equiation inside the neuron were the base (soma) is fixed as a Dirichlet condition. It would be used to obtain the gradient and them used as: the direction of the tubulin velocity (convection) or neuronal growth cones (velocity fied in linker LSConcentrationDepentVelocity() )
---@field velocity GridFunction -- Used to save the information of the velocity field from a Linker. It is used OutGrid = Export(), OutGrid:Save_Linker_to_Gridfuntion(dir, velocity)
---@field flujo_exp GridFunction -- Used to save the information of the flux from a Linker. It is used OutGrid = Export(), OutGrid:Save_Linker_to_Gridfuntion(dir, flujo_exp)
---@field flujo_outside_exp GridFunction -- Used to save the information of the flux from a Linker outside the neuron. It is used OutGrid = Export(), OutGrid:Save_Linker_to_Gridfuntion(dir, flujo_outside_exp)
---@field direc_Tubulin GridFunction -- Used to save the information of the direct tubulin from a Linker. It is used OutGrid = Export(), OutGrid:Save_Linker_to_Gridfuntion(dir, direct_tubulin)   

--- Function that Initializes and returns all required GridFunctions for the simulation
-- Constructs a table of named `GridFunction` objects using the appropriate
-- approximation spaces (e.g., `approxSpace`, `approxSpace2`, `approxSpace_Vel`)
-- required for solving the coupled PDEs in the model.
--
-- The returned table includes:
--   * Primary variables: u, u2
--   * Level set functions: lsf, lsf_oposite
--   * Signed distance fields (SDF): sdf_new, sdf_old, sdf_cfl_num_gf
--   * Solution variables: u_new, u_old
--   * Geometry-related fields: norm_vel_int, curvature
--   * Eikonal-based components: Eikonal_Inside, Eikonal_Inside_direc_inter
--   * Velocity fields: velocity, flujo_exp, flujo_outside_exp, direct_tubulin
--
-- All GridFunctions are initialized but not interpolated or filled.
-- This function assumes that the approximation spaces (`approxSpace`, `approxSpace2`, and `approxSpace_Vel`)
-- are defined globally or in the calling context.
--
-- @usage
-- local gf = GenerateGridfunctions()
-- you can use "gf.lsf" to access the primary GridFunction for the variable lsf
-- use: gf.lsf:set(0.0)
--
-- @return table gf
--   A flat table containing named GridFunctions ready to be used in the simulation
---@return GridFunctionGroup
function GenerateGridfunctions()

    -- Ensure that the approximation spaces are defined
    if not approxSpace then error("approxSpace is not defined. Please define it before calling GenerateGridfunctions.") end
    if not approxSpace2 then error("approxSpace2 is not defined. Please define it before calling GenerateGridfunctions.") end
    if not approxSpace_Vel then error("approxSpace_Vel is not defined. Please define it before calling GenerateGridfunctions.") end
    
    -- Create a table to hold all GridFunctions
    -- Each GridFunction corresponds to a specific variable in the simulation
    -- and is initialized with the appropriate approximation space.
    -- The GridFunctions are not interpolated or filled, they are just initialized.    
    local gf = {
        u = GridFunction(approxSpace), -- Main Gridfunction that have the values of the intracellular simulation (levelset < 0)
        u2 = GridFunction(approxSpace2), -- Main Gridfunction that have the values of the extracellular simulation (levelset > 0)

        lsf = GridFunction(approxSpace2), -- Gridfunction that define the values of the levelset function. This will define the interface between the intracellular and extracellular simulation. This lsf will be the input of the LSGFsimpleDomainDiscretization() for the intracellular model (negative value is the domain of the intracellular simulation).
        lsf_oposite = GridFunction(approxSpace2), -- Gridfunction that define the oposite sign of the valules of the levelset function . This have the same function of define the interface between the intracellular and extracellular simulation. This lsf_oposite will be the input of the LSGFsimpleDomainDiscretization() for the extracellular model (negative value is the domain of the extracellular simulation).

        sdf_new = GridFunction(approxSpace2), -- New values of the level set function after the reinitialize it (Used as update in EikonalDisc = HiResFluxBasedLSM() )
        sdf_old = GridFunction(approxSpace2), -- Used for solving the sdf_new
        sdf_cfl_num_gf = GridFunction(approxSpace2), -- Used to store the signed distance field for CFL number during sdf_new computation (EikonalDisc = HiResFluxBasedLSM())
        u_new = GridFunction(approxSpace2), -- After Extending the interface velocity (scalar) (ExtensionDisc = HiResFluxBasedLSM()), this is the result of the velocity that will be used to deform the interface 
        u_old = GridFunction(approxSpace2), -- Used for solving the u_new

        norm_vel_int = GridFunction(approxSpace2), -- Normal interface velocity obtained using compute_normal_vel(), where obtain the normal velocity of the velocity field (calculated from LSConcentrationDepentVelocity()) 
        curvature = GridFunction(approxSpace2), -- Curvature of the level set function, it is computed using  : Kapla_Func = Kappla_LS(); Kapla_Func:compute_kappla(lsf, curvature)

        Eikonal_Inside = GridFunction(approxSpace2), -- Values of the Stationary Poisson Equation inside the neuron (were the source term is the growth cones). It would be used to obtain the gradient and them used as: the direction of the tubulin velocity (convection) or neuronal growth cones (velocity fied in linker LSConcentrationDepentVelocity() )
        Eikonal_Inside_direc_inter = GridFunction(approxSpace2),-- Values of the Eikanol Equiation inside the neuron were the base (soma) is fixed as a Dirichlet condition. It would be used to obtain the gradient and them used as: the direction of the tubulin velocity (convection) or neuronal growth cones (velocity fied in linker LSConcentrationDepentVelocity() )

        velocity = GridFunction(approxSpace_Vel), -- Used to save the information of the velocity field from a Linker. It is used OutGrid = Export(), OutGrid:Save_Linker_to_Gridfuntion(dir, velocity) -- not best way to export velocity
        flujo_exp = GridFunction(approxSpace_Vel), -- Used to save the information of the velocity field from a Linker. It is used OutGrid = Export(), OutGrid:Save_Linker_to_Gridfuntion(flujo, flujo_exp) -- not best way to export velocity
        flujo_outside_exp = GridFunction(approxSpace_Vel), -- Used to save the information of the velocity field from a Linker. It is used OutGrid = Export(), OutGrid:Save_Linker_to_Gridfuntion(flujo_outside, flujo_outside_exp) -- not best way to export velocity
        direc_Tubulin = GridFunction(approxSpace_Vel), -- Used to save the information of the direct tubulin flux field from a Linker. It is used OutGrid = Export(), OutGrid:Save_Linker_to_Gridfuntion(direction, direct_tubulin) -- not best way to export velocity
    }

    -- Set the initial values for all the GridFunctions to zero
    for _, v in pairs(gf) do
        v:set(0.0)
    end

    return gf
end

