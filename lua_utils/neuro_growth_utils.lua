

------------------------------------------------------------------------------------------------------------
-- Function to generate the initial condition of the level-set function
------------------------------------------------------------------------------------------------------------

--- Generate initial condition for the LSF from the subset skeletum (called inner)
-- Use the function FV1_Convection() and eikanol velocity from pluglin LevelSet; and DirichletBoundary is the radius of the position of the interface.
-- @param ApproxSpace_lsf ApproxSpace of the level-set function 
-- @param gridfunction_ls Gridfunction of the level-set function (initial condition for this simulation is 0):
-- @param params table with the parameters of the simulation, including the radius of the initial geometry, which is used in the DirichletBoundary condition to set the position of the interface. The radius is defined in params.LevelSet.Initialgeometry.radius, and it is keyed by ref_def, where ref_def is 5 (numRefs ≤ 5) or 6 (numRefs ≥ 6). This allows to have different initial radii for different levels of refinement.
-- @return lsf GridFunction after the simulation
function Initialgeometry (ApproxSpace_lsf, gridfunction_ls,params)
    print("------------------------------------------------------------------------------------------")
    print("--------Iniciando los calculos de la interface inicial -----------------------------------")
    print("------------------------------------------------------------------------------------------")
    local lsf = gridfunction_ls
    local approxSpace2 = ApproxSpace_lsf

    if dim == 2 then 
        -- retornarlo a ilusolver <3
        geo_dt = 0.01 -- time step length estaba en 0.01 para ref 4 
        dt_min = 0.005
        geo_NumTimeSteps =100 -- Funcionó bien para hacer la interface en 100. --- 200 para 0.005 ; 100 para 0.01

        if numRefs >= 6 then
            geo_dt = 0.005 -- time step length
            dt_min = 0.0005
            geo_NumTimeSteps = 300

            if numRefs == 7 then
                geo_dt = 0.0025 -- time step length
                dt_min = 0.00005
                geo_NumTimeSteps = 600
            end
        end

    end

    if dim == 3 then 

        -- Funciona Ref 0,1,2,3 con preref 0
        geo_dt = 0.01 -- time step length estaba en 0.01 para ref 4 
        geo_NumTimeSteps =150 -- Funcionó bien para hacer la interface en 100. --- 200 para 0.005 ; 100 para 0.01
        dt_min = 0.005
        -- paralelo debes ponerlo en exactsolver
        if numRefs == 5 then
            geo_dt = 0.005 -- time step length
            dt_min = 0.001
            geo_NumTimeSteps = 300
        end

        if gridName == "Succes.ugx" then
            geo_dt = 0.0001 -- time step length
            geo_NumTimeSteps = 600
            dt_min = 0.00001
        end
    end    

    print("Geometría inicial :  Eikonal Equation Solver Parameters:")
    print("    dt          = " .. geo_dt)
    print("    numTimeSteps= " .. geo_NumTimeSteps)

    approxSpace2:print_statistic()
    --------------------------------------------------------------------------------
    --  Discretization
    --------------------------------------------------------------------------------
    local geo_upwind = FullUpwind()
    local geo_outer_velocity = EikonalVel(lsf)
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- 2. Exterior (phi > 0) :
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    local geo_elemDisc_out = FV1_Convection(geo_upwind, "ca_cyt", "Outer") --- Outer
    geo_elemDisc_out:set_velocity(geo_outer_velocity)
    geo_elemDisc_out:set_source(1) -- squeletum to the inner ..
    geo_elemDisc_out:set_non_sink(false)

    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- 3. Interface (phi = 0) :
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    local geo_dirichletBND = DirichletBoundary()
    if numRefs <= 5 then
        radious_negative = - 1 * params.LevelSet.Initialgeometry.radius[5]
         --- el valor siempre ha sido -0.1
    else
        radious_negative = - 1 * params.LevelSet.Initialgeometry.radius[6]
    end

    geo_dirichletBND:add(radious_negative, "ca_cyt", "Inner") -- subset : Side_Cyt

    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- 4. ------ Global discretization:  set up the global discretization
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    local geo_domainDisc = DomainDiscretization(approxSpace2)
    geo_domainDisc:add(geo_elemDisc_out)
    geo_domainDisc:add(geo_dirichletBND)

    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- SOLVER
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- este es mi smoother 
    ilu = ILU()
    ilu:set_beta(0.80) --- 
    --ilu:set_ordering_algorithm(BoostCuthillMcKeeOrdering())
    ilu:set_inversion_eps(0)
    -- ilu : agrgaer el otro error

    ilutSolver = LinearSolver()
    ilutSolver:set_preconditioner(ilu)
    ilutSolver:set_convergence_check(ConvCheck(10000, 1e-7, 1e-1, false))
  

    Base = LinearSolver()
    Base:set_convergence_check(ConvCheck(2000, 1e-7, 1e-2, false))
    Base:set_preconditioner(ilu)

    exactSolver = LU()

    -- 20260702 FIX: direct LU base solver (exactSolver) is INFEASIBLE on our fine
    -- base grid (~3M cells at base_level=numPreRefs on box_sim_soma_n40) in parallel
    -- -> it hangs ("Using ILUT ... 68921x68921" then no progress). A direct
    -- factorization of a multi-million-cell coarse grid across 128 ranks cannot
    -- complete. Use the BOUNDED iterative ILU base (Base) for 3D: it can never hang
    -- and an approximate coarse solve is fine inside a GMG preconditioner. Keep the
    -- direct solver for 2D (where the coarse grid is genuinely small).
    if dim == 3 then
        baseSolver = Base
    else
        baseSolver = exactSolver
    end
    -- Usar el iluSolver para serie en 2D o 3D; o paralelo en 2D
    --baseSolver=ilutSolver

    --smoother=vanka
    smoother=ilu
    --numSmooth=2*numRefs
    numPreSmooth=2
    numPostSmooth=2
    gmg = GeometricMultiGrid(approxSpace2)
    gmg:set_discretization(geo_domainDisc)
    gmg:set_base_level(numPreRefs)
    gmg:set_base_solver(baseSolver) -- dimitry sugerio el exactsolver
    gmg:set_smoother(smoother)
    gmg:set_cycle_type("V")
    gmg:set_num_presmooth(numPreSmooth)
    gmg:set_num_postsmooth(numPostSmooth)
    gmg:set_rap(true)
    --gmg:set_damp(MinimalResiduumDamping())
    --gmg:set_damp(0.9)
    --gmg:set_damp(MinimalEnergyDamping())
    --gmg:set_debug(dbgWriter)


    BiCGStabSolver = BiCGStab()
    -- 20260702 FIX (research-backed): for 3D use a SINGLE-LEVEL ILU preconditioner,
    -- NOT GMG. GMG's coarse (base) solve on our fine n40 grid hangs, and making the
    -- base iterative turns the preconditioner VARIABLE -> BiCGStab's Krylov recurrence
    -- becomes formally invalid and stalls. A fixed ILU preconditioner is bounded and
    -- cannot hang (approximate interior direction field is fine for growth). 2D keeps GMG.
    if dim == 3 then
        BiCGStabSolver:set_preconditioner(ilu)
    else
        BiCGStabSolver:set_preconditioner(gmg)
    end
    BiCGStabSolver:set_convergence_check(ConvCheck(100, 1e-9, 1e-1, true))
    --BiCGStabSolver:set_compute_fresh_defect_when_finished(true)

    -- create Linear Solver
    convCheck = ConvCheck(1000, 1e-10, 1e-2, true)
    convCheck:set_verbose(false)

    gmgSolver = LinearSolver()
    gmgSolver:set_preconditioner(gmg)
    gmgSolver:set_convergence_check(ConvCheck(100, 1e-8, 1e-02, true))
    -- create Exact solver


    -- choose a solver
    --solver = exactSolver
    --solver = baseSolver
    geo_solver = BiCGStabSolver
    --solver = gmgSolver --- este era el bueno
    --solver = ilutSolver
    --solver = linSolver

    --------------------------------------------------------------------------------
    --	Output facilities
    --------------------------------------------------------------------------------
    --local geo_vtkOut = VTKOutput();
    --geo_vtkOut:clear_selection();
    --geo_vtkOut:select_nodal("ca_cyt", "distFunc"); -- distance function es la variable "distFunc"
    --local geo_vtkOut = VTKOutput();
    --geo_vtkOut:clear_selection();
    --geo_vtkOut:select_nodal("ca_cyt", "distFunc");
    --local VTK_name = "Initial"
     geo_vtkOut = nil
     VTK_name = nil
    --------------------------------------------------------------------------------
    --  Initial condition
    --------------------------------------------------------------------------------
    lsf:set(0);
    if gridName ~= "Succes.ugx" then
        Interpolate(1, lsf, "ca_cyt", "Outer,Outflow,Bottom")
    elseif gridName == "Succes.ugx" then
        Interpolate(1, lsf, "ca_cyt", "Outer,Outflow")
    end

    Interpolate(radious_negative, lsf, "ca_cyt", "Inner")
    print("Geometría inicial : Initial condition created ----------------------------------- \n")
    ---------------------------------c-----------------------------------------------
    --  Timestepping
    --------------------------------------------------------------------------------

    util.SolveLinearTimeProblem(lsf, {
        domainDisc = geo_domainDisc,
        reassemble = true
    }, geo_solver, geo_vtkOut, VTK_name, "ImplEuler", 1, 0, nil, geo_dt, dt_min, nil, false, {
        preProcess = geo_dwOrder
    }, 0, geo_NumTimeSteps);

    print("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    print("-------- Finalizado los calculos de la interface inicial ----------------------------------------------")
    print("------------------------------------------------------------------------------------------------ \n\n\n")
    return lsf
end






-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
--  Function to obtain the direction of the tubuline inside the neuron
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- @param ApproxSpace_lsf ApproxSpace of the level-set function
-- @param gridfunction_gradiente Gridfunction of the gradient of the level-set function, which will be used to obtain the direction of the tubuline inside the neuron. The gradient is obtained from the solution of the Eikonal equation, which is solved in this function.
-- @param gridfunction_ls Gridfunction of the level-set function
-- @param domain Domain of the simulation
-- @param inicial_pregunta boolean that indicates if it is the first time that the function is called, which is used to set the parameters of the Eikonal equation solver. If it is the first time, the parameters are set to generate the gradient in the whole geometry, and if it is not the first time, the parameters are set to generate the gradient only in a part of the geometry, which is enough to have the direction of the tubuline inside the neuron.
-- @param params Table of parameters for the simulation
-- @return gridfunction_gradiente Gridfunction of the gradient of the level-set function, which will be used to obtain the direction of the tubuline inside the neuron. The gradient is obtained from the solution of the Eikonal equation, which is solved in this function.
function GenerateImaginaryMolecule2(ApproxSpace_lsf, gridfunction_gradiente, gridfunction_ls, domain, inicial_pregunta, params)
    print("Start: Eikanol inside")
    
    -- 1. Parameters ----------------------------------------------------------------- 
    local inicial  = inicial_pregunta
    if numRefs == 6 or numRefs == 5 then
        geo_dt = 0.01 -- time step length
        geo_NumTimeSteps = 200
        geo_endTime = 2 -- 0.01 * 200
        if inicial then  -- si es la primera vez, se genera toda la gradiente 
            geo_NumTimeSteps = 300 -- Dependiendo de la geometria debes de modificar este valor para asegurarte que la advencion sea suficiente para que vaya por toda la geometria
            geo_endTime = 3 -- 0.01 * 300
        end        
        if dim == 3 then
            geo_dt = 0.001 -- time step length
            geo_NumTimeSteps = 1000 -- Dependiendo de la geometria debes de modificar este valor para asegurarte que la advencion sea suficiente para que vaya por toda la geometria
            geo_endTime = 1 -- 0.005 * 300
            if inicial then  -- si es la primera vez, se genera toda la gradiente 
                geo_NumTimeSteps = 2000 -- Dependiendo de la geometria debes de modificar este valor para asegurarte que la advencion sea suficiente para que vaya por toda la geometria
                geo_endTime = 2 -- 0.01 * 300
            end  
        end
        dt_min = 0.00001 -- solo prueba para funcion
    end

    if numRefs == 7 then
        geo_dt = 0.005 -- time step length
        geo_NumTimeSteps = 400
        geo_endTime = 2 -- 0.01 * 200
        if inicial then  -- si es la primera vez, se genera toda la gradiente 
            geo_NumTimeSteps = 600 -- Dependiendo de la geometria debes de modificar este valor para asegurarte que la advencion sea suficiente para que vaya por toda la geometria
            geo_endTime = 3 -- 0.01 * 300
        end        
        dt_min = 0.00001 -- solo prueba para funcion
    end

    if numRefs == 4 or numRefs <= 3 then
        geo_dt = 0.1 -- time step length
        geo_NumTimeSteps = 100 -- Dependiendo de la geometria debes de modificar este valor para asegurarte que la advencion sea suficiente para que vaya por toda la geometria
        geo_endTime = 10 -- 0.1 * 100
        
        if inicial then  -- si es la primera vez, se genera toda la gradiente 
            geo_NumTimeSteps = 150 --200 estaba  Dependiendo de la geometria debes de modificar este valor para asegurarte que la advencion sea suficiente para que vaya por toda la geometria
            geo_endTime = 15 -- 0.1 * 150
        end
        
        if dim == 3 then
            geo_dt = 0.01 -- time step length
            geo_NumTimeSteps = 200 -- Dependiendo de la geometria debes de modificar este valor para asegurarte que la advencion sea suficiente para que vaya por toda la geometria
            geo_endTime = 2 -- 0.01 * 200
            if inicial then  -- si es la primera vez, se genera toda la gradiente 
                geo_NumTimeSteps = 300 -- Dependiendo de la geometria debes de modificar este valor para asegurarte que la advencion sea suficiente para que vaya por toda la geometria
                geo_endTime = 3 -- 0.01 * 300
            end  
        end
        dt_min = 0.0001 -- solo prueba para funcion
    end

    -- 20260706: the geo_dt/geo_endTime above are chosen by a numRefs branch that treats
    -- numRefs<=3 as "coarse" (dt 0.01). But the real neuron's box is already fine at ref2, so
    -- dt 0.01 violates CFL for the eikonal advection -> Newton fails -> dt death-spirals and the
    -- solve never finishes. Allow CLI overrides to match dt/endTime to the actual mesh.
    geo_dt      = util.GetParamNumber("-eikDt",      geo_dt,      "eikonal substep dt (override; lower for fine meshes)")
    geo_endTime = util.GetParamNumber("-eikEndTime", geo_endTime, "eikonal end time (override; front must reach the tips)")
    geo_NumTimeSteps = math.ceil(geo_endTime / geo_dt)
    print(string.format("[eikonal] geo_dt=%.5g endTime=%.5g NumTimeSteps=%d", geo_dt, geo_endTime, geo_NumTimeSteps))

    -- 2. Load a Domain -----------------------------------------------------------------
    local dom = domain
    -- 3. Approximation Space: set up approximation space ---------------------------------
    local approxSpace = ApproxSpace_lsf
    ----------------------------------------------------------------------------------------------
    -- Gridfuncctions
    -----------------------------------------------------------------------------------------------
    -- 1. Grid functions for the solution
    local lsf = gridfunction_ls
    -- 2. Grid function for the solution of the solution (mol)
    local u = gridfunction_gradiente
    
    ----------------------------------------------------------------------------------------------
    -- PARTIAL DIFFERENTIAL EQUATIONS 
    ----------------------------------------------------------------------------------------------
    local geo_upwind = FullUpwind()
    local geo_inner_velocity = EikonalVel(u)
        -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- 1. Interior (phi < 0) :
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    local geo_elemDisc_inner = FV1_Convection(geo_upwind, "ca_cyt", "Outer") --- Inner
    geo_elemDisc_inner:set_velocity(geo_inner_velocity)
    geo_elemDisc_inner:set_source(1.0)
    geo_elemDisc_inner:set_diffusion(params.LevelSet.Eikanolinside.diffusion) -- this will smooth the eikonal to make it parallel
    geo_elemDisc_inner:set_non_sink(false)

    if inicial then 
        geo_elemDisc_inner:set_diffusion(params.LevelSet.Eikanolinside.initialdiffusion) -- estaba en 0.06
    end


    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- 3. Interface (phi = 0) :
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    local geo_dirichletBND = DirichletBoundary()
    geo_dirichletBND:add(0, "ca_cyt", "Outflow,Bottom") -- subset : Side_Cyt

    ----------------------------------------------------------------------------
    --- Domain Discretization
    ----------------------------------------------------------------------------
    gradiente_domainDisc = LSGFsimpleDomainDiscretization(approxSpace) -- Free surface: Using the simple extrapolation for the embedded interface
    gradiente_domainDisc:add(geo_elemDisc_inner)
    gradiente_domainDisc:add(geo_dirichletBND)

    gradiente_domainDisc:set_LSF(lsf) -- Set the Level-Set Function
    gradiente_domainDisc:set_Neumann0_on_if_for("ca_cyt") --- indicando not cruzar el dominio 
    gradiente_domainDisc:project_LSF()
    -----------------------------------------------------------------------------------------------
    -- SOLUTION
    -----------------------------------------------------------------------------------------------
    --Interpolate(0, u, "ca_cyt", "Inner")
    --Interpolate(0, u, "ca_cyt", "Outer")
    --Interpolate(0, u, "ca_cyt", "Outflow")
    --Interpolate(0, u, "ca_cyt", "Bottom") -- es como la interface 

    print("Initial condition created ----------------------------------- \n\n")

    -- set up solver (using 'util/solver_util.lua')
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- SOLVER
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    --util.solver.defaults.approxSpace = approxSpace --- yes

    local geo_solverDesc3D = { --- yes
        type = "newton", -- nombre del metodo de discretization 
    
        lineSearch = {
            type = "standard", -- ["standard", "none"]
            maxSteps = 10, -- maximum number of line search steps
            lambdaStart = 1, -- start value for scaling parameter
            lambdaReduce = 0.7, -- reduction factor for scaling parameter
            acceptBest = true, -- check for best solution if true
            checkAll = false -- check all maxSteps steps if true 
        },
    
        convCheck = {
            type = "standard",
            iterations = 50, -- maximum number of iterations
            -- 20260706 FIX: was 5e-5. At ref2 the eikonal's initial Newton defect (~2.7e-5) is
            -- already BELOW 5e-5, so Newton "converges in 0 steps" and the field stays ZERO ->
            -- no axial direction -> real neuron forced onto the outward normal -> radial thickening.
            -- Tightening below the ref2 initial defect forces a real solve (nonzero eikonal at ref2).
            absolute = util.GetParamNumber("-eikNewtonAbs", 1e-7, "eikonal Newton abs tol (5e-5 was too loose -> 0 field at ref2)"),
            reduction = 1e-8, -- reduction factor of defect to be reached; usually 1e-6 - 1e-8
            verbose = true -- print convergence rates if true
        },
    
        linSolver = {
            type = "bicgstab",
            precond = {
                type = "gmg", -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
                approxSpace = approxSpace,
                smoother = {
                    type = "ilu",
                    beta            = 0.8,
                    damping         = 1,
                    sortEps         = 1.e-50,
                    inversionEps    = 0
                    }, -- pre- and postsmoother, only for gmg ["ilu", "ilut", "jac", "gs", "sgs"]cycle = "V", -- gmg-cycle ["V", "F", "W"]
                preSmooth = 3, -- number presmoothing steps
                postSmooth = 3, -- number postsmoothing steps
                rap = true, -- comutes RAP-product instead of assembling if true
                baseLevel = numPreRefs, -- gmg - coarsest level
                -- 20260706 FIX: was "lu". The real-neuron input grid is already fine, so the GMG
                -- base level (numPreRefs=0) has ~69k DOFs. A gathered DIRECT LU on it took 880 s
                -- PER V-cycle -> ~15 min/substep -> the eikonal never finished at ref2/96. An
                -- ITERATIVE base solver (BiCGStab+ILU) solves the coarse system in seconds.
                baseSolver = {
                    type = "bicgstab",
                    precond = "ilu",
                    convCheck = { type = "standard", iterations = 500,
                                  absolute = 1e-9, reduction = 1e-8, verbose = false }
                }
            },
    
            convCheck = {
                type = "standard",
                iterations = 150, -- number of iterations
                absolute = 1e-8, -- absolut value of defact to be reached; usually 1e-8 - 1e-10 (may not be larger than in newton section)
                reduction = 1e-12, -- reduction factor of defect to be reached
                verbose = true -- print convergence rates if true
            }
        }
    } -- definir el error para que ecuaciones quieres ... 

    local geo_solverDesc2D = { --- yes
    type = "newton", -- nombre del metodo de discretization 

    lineSearch = {
        type = "standard", -- ["standard", "none"]
        maxSteps = 8, -- maximum number of line search steps
        lambdaStart = 1, -- start value for scaling parameter
        lambdaReduce = 0.5, -- reduction factor for scaling parameter
        acceptBest = true, -- check for best solution if true
        checkAll = false -- check all maxSteps steps if true 
    },

    convCheck = {
        type = "standard",
        iterations = 50, -- maximum number of iterations
        absolute = 5e-7, -- absolut value of defect to be reached; usually 1e-7 - 1e-9
        reduction = 1e-10, -- reduction factor of defect to be reached; usually 1e-6 - 1e-8
        verbose = true -- print convergence rates if true
    },

    linSolver = {
        type = "bicgstab",
        precond = {
            type = "gmg", -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
            approxSpace = approxSpace,
            smoother = "ilu", -- pre- and postsmoother, only for gmg ["ilu", "ilut", "jac", "gs", "sgs"]
            cycle = "V", -- gmg-cycle ["V", "F", "W"]
            preSmooth = 2, -- number presmoothing steps
            postSmooth = 2, -- number postsmoothing steps
            rap = true, -- comutes RAP-product instead of assembling if true
            baseLevel = numPreRefs, -- gmg - coarsest level
            baseSolver = "lu"
        },

        convCheck = {
            type = "standard",
            iterations = 150, -- number of iterations
            absolute = 1e-8, -- absolut value of defact to be reached; usually 1e-8 - 1e-10 (may not be larger than in newton section)
            reduction = 1e-12, -- reduction factor of defect to be reached
            verbose = true -- print convergence rates if true
        }
    }
    } -- definir el error para que ecuaciones quieres ... 
    if dim == 2 then
        geo_solverDesc = geo_solverDesc2D
    elseif dim == 3 then
        geo_solverDesc = geo_solverDesc3D
    end
    local geo_newtonSolver
    -- 20260706: BUG #4 fix. GMG diverges for the convection-dominated eikonal at usable dt (defect
    -- explodes 2.7e-4->0.96, rate 3.5e3), so only tiny dt is stable -> thousands of substeps ->
    -- infeasible at ref2. Replace GMG with a convection-ROBUST linear solver: GMRES (robust for
    -- nonsymmetric) + ILU (with Cuthill-McKee ordering + beta), NO coarse-grid correction. Toggle
    -- with -eikConvSolver so the default path is untouched and this is fully reversible.
    if util.HasParamOption("-eikConvSolver") then
        local beta    = util.GetParamNumber("-eikIluBeta",     0.5,  "eikonal ILU beta (MILU); higher = stronger")
        local rst     = util.GetParamNumber("-eikGmresRestart",50,   "eikonal GMRES restart depth")
        local liniter = util.GetParamNumber("-eikLinIters",    600,  "eikonal linear max iterations")
        local nabs    = util.GetParamNumber("-eikNewtonAbs",   1e-7, "eikonal Newton abs tol")
        local ilu = ILU(); ilu:set_beta(beta); ilu:set_sort(true); ilu:set_inversion_eps(0)
        if not util.HasParamOption("-eikNoOrder") then
            ilu:set_ordering_algorithm(NativeCuthillMcKeeOrdering())
        end
        local lin = GMRES(rst); lin:set_preconditioner(ilu)
        lin:set_convergence_check(ConvCheck(liniter, 1e-9, 1e-10, false))
        local nwt = NewtonSolver()
        nwt:set_linear_solver(lin)
        nwt:set_convergence_check(ConvCheck(50, nabs, 1e-8, true))
        nwt:set_line_search(StandardLineSearch(10, 1.0, 0.7, true, false))
        geo_newtonSolver = nwt
        print(string.format("[eikonal] CONVECTION solver: GMRES(%d)+ILU(beta=%.2g,CuthillMcKee) newtonAbs=%.1g", rst, beta, nabs))
    else
        geo_newtonSolver = util.solver.CreateSolver(geo_solverDesc) --- yes
    end

    local geo_vtkOut = VTKOutput();
    geo_vtkOut:clear_selection();
    geo_vtkOut:select_nodal("ca_cyt", "eikanol"); -- distance function es la variable "distFunc"
    geo_name = folder .. "eikanol"
    --geo_vtkOut:print(folder .. "eikanol", u, 0, 0) --- antes era time

    if dim == 2 then
        --------------------------------
        --  Setup Time Discretization --
        --------------------------------
        local timeDisc_Gradient = util.CreateTimeDisc(gradiente_domainDisc, "ImplEuler", 1)
        timeDisc_Gradient:set_stage(1)

        --------------------
        --  Algebra setup --
        --------------------
        geo_newtonSolver:init(AssembledOperator(timeDisc_Gradient))

        -------------------
        -- Time stepping --
        -------------------
        local time = 0
        -- create new grid function for old value
        local uold = u:clone() -- create a copy of the solution
        local solTimeSeries = SolutionTimeSeries()
        solTimeSeries:push(uold, time)
        ---------------------------------c-----------------------------------------------
        --  Timestepping
        --------------------------------------------------------------------------------
        for step = 1, geo_NumTimeSteps do --- 
            timeDisc_Gradient:prepare_step(solTimeSeries, geo_dt)
            geo_newtonSolver:prepare(u)
            geo_newtonSolver:apply(u)

            -- update new time
            time = timeDisc_Gradient:future_time()
            --geo_vtkOut:print(folder .. "gradientess", u, step, time) --- antes era time

            -- update time series
            local oldestSol = solTimeSeries:oldest() -- obtain  the reference the oldest solution
            VecAssign(oldestSol, u) -- ssign the neew to the oldes - reinitialize the oldest solution with the new values
            solTimeSeries:push_discard_oldest(oldestSol, time) -- push the oldest solition swith the new values  to the front, oldest sol pointer is poped  fromn end
            --  gradiente_domainDisc:project_LSF()

        end

    elseif dim == 3 then
    -- 20260706: the eikonal solve wrote an "eikanol" VTK frame EVERY substep (~geo_NumTimeSteps
    -- x 96 parts per growth step) -> thousands of useless intermediate files (inode blow-up).
    -- These intermediate dumps are never used (analysis uses the final Eikonal_Inside field).
    -- Pass nil out/name unless -eikanolDebug is set.
    local _eikOut, _eikName = nil, nil
    if util.HasParamOption("-eikanolDebug") then _eikOut, _eikName = geo_vtkOut, geo_name end
    util.SolveNonlinearTimeProblem(u, gradiente_domainDisc, geo_newtonSolver, _eikOut, _eikName , "ImplEuler", 1, 0, geo_endTime, geo_dt, dt_min, 0.5);

    end

    --[[ -----------------
    -- este es mi smoother 
    ilu = ILU()
    ilu:set_beta(0.80) --- 
    --ilu:set_ordering_algorithm(BoostCuthillMcKeeOrdering())
    ilu:set_inversion_eps(0)
    -- ilu : agrgaer el otro error

    ilutSolver = LinearSolver()
    ilutSolver:set_preconditioner(ilu)
    ilutSolver:set_convergence_check(ConvCheck(10000, 1e-7, 1e-1, false))
  

    Base = LinearSolver()
    Base:set_convergence_check(ConvCheck(2000, 1e-7, 1e-2, false))
    Base:set_preconditioner(ilu)
    
    exactSolver = LU()

    --baseSolver=Base
      baseSolver=exactSolver
    --baseSolver=ilutSolver

    --smoother=vanka
    smoother=ilu
    --numSmooth=2*numRefs
    numPreSmooth=2
    numPostSmooth=2
    gmg = GeometricMultiGrid(approxSpace)
    gmg:set_discretization(gradiente_domainDisc)
    gmg:set_base_level(numPreRefs)
    gmg:set_base_solver(baseSolver) -- dimitry sugerio el exactsolver
    gmg:set_smoother(smoother)
    gmg:set_cycle_type("V")
    gmg:set_num_presmooth(numPreSmooth)
    gmg:set_num_postsmooth(numPostSmooth)
    gmg:set_rap(true)
    --gmg:set_damp(MinimalResiduumDamping())
    --gmg:set_damp(0.9)
    --gmg:set_damp(MinimalEnergyDamping())
    --gmg:set_debug(dbgWriter)

    gmgSolver = LinearSolver()
    gmgSolver:set_preconditioner(gmg)
    gmgSolver:set_convergence_check(ConvCheck(100, 1e-8, 1e-02, true))
    -- create Exact solver

    geo_solver = gmgSolver --- este era el bueno

    newtonConvCheck = ConvCheck()
    newtonConvCheck:set_maximum_steps(10000)
    newtonConvCheck:set_minimum_defect(1e-6)
    newtonConvCheck:set_reduction(1e-10)
    newtonConvCheck:set_verbose(true)

    op = AssembledOperator(gradiente_domainDisc)
    op:init()

    newtonSolver = NewtonSolver(op)
    newtonSolver:set_linear_solver(geo_solver)
    newtonSolver:set_convergence_check(newtonConvCheck)  

    newtonSolver:init(op)
    newtonSolver:prepare(u)
    newtonSolver:apply(u)
    --]]
    print("finish: Eikanol inside")
    return u
end

-- =============================================================================================
--  HEAT METHOD for the interior growth direction  (Crane, Weischedel & Wardetzky)
--  20260713 — replaces the eikonal ADVECTION (GenerateImaginaryMolecule2) above.
--
--  WHY: the eikonal solves  dE/dt + div(E*V) = 1,  V = grad(E)/|grad(E)|  by advecting to steady
--  state. That operator is CONVECTION-DOMINATED, so GMG diverges at ref2 (we had to fall back to
--  GMRES+ILU) and it is CFL-limited => hundreds of nonlinear substeps => ~50 min per growth step.
--
--  KEY INSIGHT: every velocity linker NORMALIZES the direction
--    (ls_concentration_dependent_velocity_*linker.h: VecScale(dir, grad, 1/(|grad|+1e-7)))
--  => only the DIRECTION of grad(E) matters; |grad(E)|=1 ("eikonal cleanliness") is irrelevant.
--  So ANY scalar field whose gradient points soma->tip is a valid drop-in.
--
--  THE METHOD (all steps are ELLIPTIC/SPD => GMG is optimal; no CFL, no Newton, no convection):
--    1. screened Poisson:  -t*Lap(u) + u = 0 ,  u = 1 Dirichlet on the SOMA subset ("Bottom"),
--       Neumann-0 at the membrane (LSGF) => u ~ exp(-d_geodesic/sqrt(t)), maximal at the soma.
--    2. X = -grad(u)/|grad(u)|   (points AWAY from the soma = soma->tip)  = EikonalVel(u, scaling=-1)
--    3. Poisson recovery:  Lap(phi) = div(X) ,  phi = 0 on the soma  => phi ~ geodesic distance,
--       grad(phi) points soma->tip. Consumed exactly like the eikonal was.
--
--  *** THE SIGN TRAP ***  UG4's ConvectionDiffusion convention (convection_diffusion_base.h:57-58) is
--        d/dt(m1*c) - div( D*grad(c) - v*c - F ) + r1*c + r2 = f + div(f2)
--  so with set_diffusion(1) it assembles  -Lap(phi) = div(f2).  Crane needs  Lap(phi) = div(X),
--  hence  f2 = -X = +grad(u)/|grad(u)|  =>  set_vector_source( EikonalVel(u):set_scaling(+1) ).
--  Using -1 here silently yields a direction pointing INTO the soma (retraction/thickening that
--  looks like a physics bug). -heatSign flips it; validate_heat_direction reports the sign.
--
--  CHOOSING t:  Crane's default t ~ h^2 is WRONG here. With sqrt(t) ~ h the field decays ~0.38 per
--  cell; at ref2 (L/h ~ 128) the tip value is ~1e-53 — far below the linear solver's residual floor,
--  so grad(u) at the distal tips would be pure Krylov NOISE. Instead set the smoothing length
--        l = sqrt(t) = max(4h, L/10)          =>  t = l^2      (h0=0.125, L=5  =>  l=0.5, t=0.25)
--  which keeps u_tip/u_soma ~ e^-10 = 4.5e-5 (well above a 1e-10 solve) AND resolves the layer
--  (l/h >= 4). Large t is safe: inside a TUBE the level sets of any soma-anchored diffusion are
--  cross-sections, so the gradient stays AXIAL regardless of t — and the linkers normalize anyway.
-- =============================================================================================

--- Solver descriptor for the heat/Poisson solves: the SAME bicgstab+GMG block that
--- GenerateImaginaryMolecule already runs successfully at ref2. Both operators here are SPD
--- (-t*Lap+1 and -Lap+eps), so GMG is optimal — -eikConvSolver (GMRES+ILU) is NOT needed.
-- The soma heat source: a ball centred on the soma. Consumed as LuaUserNumber("HeatSomaSource") by
-- GenerateHeatDirection (which sets the HEAT_SOMA_* globals). Signature is UG4's (x, y, z, t).
HEAT_SOMA_X, HEAT_SOMA_Y, HEAT_SOMA_Z, HEAT_SOMA_R, HEAT_SOMA_S = 0.0, 0.0, 0.0, 0.28, 1.0
function HeatSomaSource(x, y, z, t)
    local dx, dy, dz = x - HEAT_SOMA_X, y - HEAT_SOMA_Y, z - HEAT_SOMA_Z
    if dx * dx + dy * dy + dz * dz <= HEAT_SOMA_R * HEAT_SOMA_R then return HEAT_SOMA_S end
    return 0.0
end

-- Both heat solves are SPD. Step 1 (screened Poisson, strongly diagonally dominant) is happy with
-- BiCGStab+ILU. Step 3 (the Poisson recovery) is only weakly regularized, and T0 showed BiCGStab
-- BREAKING DOWN on it ("alpha = 0") -> phi came back as a garbage partial iterate, which is what wrecked
-- the direction (G=0.12 vs the eikonal's 0.82), NOT the huge ||phi|| (a constant offset cannot change
-- grad(phi)). For SPD systems CG cannot break down -- but it requires a SYMMETRIC preconditioner, so
-- pair it with a symmetric smoother (sgs), never ILU.
function HeatSolverDesc(approxSpace, absTol, krylov, smoother)
    return {
        type = krylov or "bicgstab",
        precond = {
            type = "gmg",
            approxSpace = approxSpace,
            smoother = smoother or "ilu",
            cycle = "V",
            preSmooth = 2,
            postSmooth = 2,
            rap = true,
            baseSolver = (dim == 3) and ((util.GetParamNumber("-robustInteriorSolve", 0, "3D interior GMG base solver: 0=BiCGStab+ILU, 1=LU") == 1) and "lu" or { type = "bicgstab", precond = { type = "ilu" }, convCheck = { type = "standard", iterations = 200, absolute = 1e-9, reduction = 1e-2, verbose = false } }) or "lu",
            baselevel = numPreRefs
        },
        convCheck = {
            type = "standard",
            iterations = 150,
            absolute = absTol or 1e-10,   -- tight: u at the tips is ~4.5e-5, must not be solver noise
            reduction = 1e-12,
            verbose = false
        }
    }
end

--- Drop-in replacement for GenerateImaginaryMolecule2 (same signature).
--- Fills `gridfunction_gradiente` with a field whose gradient points soma -> tip.
function GenerateHeatDirection(ApproxSpace_lsf, gridfunction_gradiente, gridfunction_ls, domain, inicial_pregunta, params)
    print("Start: HEAT METHOD (interior direction)")
    local t0 = os.time()

    local approxSpace = ApproxSpace_lsf
    local lsf         = gridfunction_ls
    local phi         = gridfunction_gradiente   -- result, written in place

    ------------------------------------------------------------------
    -- 0. diffusion time t  (l = sqrt(t) is the smoothing length)
    ------------------------------------------------------------------
    local h0  = util.GetParamNumber("-heatH0",  0.125, "base mesh size h at numRefs=0 (box 5.0 / 40 cells)")
    local h   = h0 / math.pow(2, numRefs)
    local L   = util.GetParamNumber("-heatL",   5.0,   "longest soma->tip geodesic (domain units); only used by the auto-t rule")
    local ell = util.GetParamNumber("-heatLen", -1.0,  "heat length l=sqrt(t); -1 = auto max(4h, L/10)")
    if ell <= 0 then
        local tCLI = util.GetParamNumber("-heatTime", -1.0, "heat diffusion time t; -1 = auto")
        if tCLI > 0 then ell = math.sqrt(tCLI) else ell = math.max(4.0 * h, L / 10.0) end
    end
    local tHeat = ell * ell

    local somaSub  = util.GetParam("-heatSoma", "Bottom", "Dirichlet subset holding the soma value (NOT Outflow!)")
    local heatSign = util.GetParamNumber("-heatSign", 1.0, "ESCAPE HATCH: flip to -1 if the arrows point INTO the soma")
    local shortcut = util.HasParamOption("-heatShortcut", "1-solve variant: return -u directly (grad already soma->tip); skips the Poisson recovery")
    -- Regularisation of the Poisson recovery. 1e-6 left the operator so close to singular that the Krylov
    -- solver broke down (T0). 1e-3 screens over 1/sqrt(reg) ~ 32 >> L=5, so it does NOT distort the field,
    -- but it conditions the matrix. (An additive constant is harmless anyway: grad(phi) ignores it.)
    local reg      = util.GetParamNumber("-heatReg", 1e-3, "reaction-rate regularisation of the Poisson recovery (conditions the near-singular pure-Neumann operator)")

    print(string.format("[heat] t=%.5g  l=sqrt(t)=%.5g  h=%.5g  l/h=%.2f  L/l=%.2f  soma='%s'  sign=%+.0f  %s",
                        tHeat, ell, h, ell/h, L/ell, somaSub, heatSign, shortcut and "SHORTCUT(1 solve)" or "3-step(2 solves)"))
    if ell/h < 2.0 then print("[heat] WARN: l/h < 2 -> the screened-Poisson boundary layer is UNDER-RESOLVED (FV1 may oscillate). Raise -heatLen.") end
    if L/ell > 15.0 then print("[heat] WARN: L/l > 15 -> u at the distal tips ~exp(-15)=3e-7, near the solver noise floor. Raise -heatLen (3-step self-heals).") end
    if L/ell < 1.0  then print("[heat] WARN: L/l < 1 -> u is nearly CONSTANT; the gradient may be round-off dominated. Lower -heatLen.") end

    ------------------------------------------------------------------
    -- STEP 1: screened Poisson   -t*Lap(u) + u = 0 ,  u = somaVal on the soma subset
    --         Neumann-0 at the membrane => GEODESIC (inside-the-neuron) field.  SPD.
    ------------------------------------------------------------------
    local u_heat = GridFunction(approxSpace)
    u_heat:set(0.0)

    ------------------------------------------------------------------
    -- HOW THE SOMA DRIVES THE HEAT: a localized VOLUME SOURCE, *not* a Dirichlet BC.
    --
    -- T0 (job 12811310) proved a nonzero Dirichlet on the soma subset does NOT couple through the
    -- LSGF projection: the assembled RHS came out identically zero, so BiCGStab hit
    --   "Method breakdown tt = 0.000000e+00"  (zero initial residual = nothing to solve) and u_heat == 0.
    -- Note that EVERY working solve in this file (GenerateImaginaryMolecule's set_vector_source, the
    -- eikonal's set_source(1.0)) is driven by a SOURCE, never by a nonzero Dirichlet. So we use the
    -- forcing path that is proven to assemble: a source ball centred on the soma.
    --   -t*Lap(u) + u = S  inside the ball,  0 outside,  Neumann-0 at the membrane.
    -- Physically identical (heat injected at the soma, decaying with GEODESIC distance) and it is also
    -- robust to the other candidate cause -- the soma subset's nodes lying outside the level set --
    -- because the ball is centred on the soma CENTRE, which is interior by construction.
    -- No Dirichlet at all: the screening term (+u) already makes the operator SPD and non-singular.
    ------------------------------------------------------------------
    HEAT_SOMA_X = util.GetParamNumber("-heatSomaX", 0.120585, "soma centre x (default: real pyramidal)")
    HEAT_SOMA_Y = util.GetParamNumber("-heatSomaY", 0.270423, "soma centre y")
    HEAT_SOMA_Z = util.GetParamNumber("-heatSomaZ", -0.339068, "soma centre z")
    -- >= 2h so the ball always captures several nodes, even at ref0 (h = 0.125)
    HEAT_SOMA_R = util.GetParamNumber("-heatSomaR", math.max(0.28, 2.0 * h), "soma source-ball radius")
    -- SHORTCUT: a NEGATIVE source makes u most-negative at the soma and RISING outward, so grad(u)
    -- already points soma->tip and u is a valid drop-in with no recovery solve.
    HEAT_SOMA_S = shortcut and (-1.0 * heatSign) or 1.0

    local heatEq = ConvectionDiffusionFV1("ca_cyt", "Outer")
    heatEq:set_diffusion(tHeat)     -- D = t
    heatEq:set_reaction_rate(1.0)   -- r1 = 1  =>  -t*Lap(u) + u = S   (the +u removes the pure-Laplace degeneracy)
    heatEq:set_source(LuaUserNumber("HeatSomaSource"))   -- the soma ball (globals above)

    local heatDisc = LSGFsimpleDomainDiscretization(approxSpace)
    heatDisc:add(heatEq)
    heatDisc:set_LSF(lsf)
    heatDisc:set_Neumann0_on_if_for("ca_cyt")   -- no-flux at the membrane => geodesic, not Euclidean
    heatDisc:project_LSF()

    print(string.format("[heat] soma ball: c=(%.4g,%.4g,%.4g) R=%.4g S=%+.0f", HEAT_SOMA_X, HEAT_SOMA_Y,
                        HEAT_SOMA_Z, HEAT_SOMA_R, HEAT_SOMA_S))

    local solver1 = util.solver.CreateSolver(HeatSolverDesc(approxSpace, 1e-10))
    local A1 = AssembledLinearOperator(heatDisc)
    local b1 = GridFunction(approxSpace); b1:set(0.0)
    heatDisc:adjust_solution(u_heat)
    heatDisc:assemble_linear(A1, b1)
    local nb = VecNorm(b1)
    solver1:init(A1, u_heat)
    solver1:apply(u_heat, b1)

    local nu = VecNorm(u_heat)
    print(string.format("[heat] step1 (screened Poisson) done   ||b|| = %.6g   ||u_heat|| = %.6g", nb, nu))
    if nb < 1e-14 then
        print("[heat] FATAL: the assembled RHS is ZERO -> the soma source ball caught no interior node. " ..
              "Check -heatSomaX/Y/Z/R against the geometry (this is the T0 failure mode).")
    end

    ------------------------------------------------------------------
    -- SHORTCUT (1 solve): return u directly; grad(u) already points soma->tip (somaVal = -1)
    ------------------------------------------------------------------
    if shortcut then
        phi:assign(u_heat)
        gradiente_domainDisc = heatDisc   -- GLOBAL: consumed by GridFuncLSGradientData in the model
        print(string.format("finish: HEAT METHOD (shortcut, 1 solve) in %d s", os.time() - t0))
        return phi
    end

    ------------------------------------------------------------------
    -- STEP 2+3 (DEFAULT): Poisson  Lap(phi) = div(X),  X = -grad(u)/|grad(u)|
    --   UG4 assembles  -Lap(phi) = div(f2)   =>   f2 = -X = +grad(u)/|grad(u)|
    --   => vector source = EikonalVel(u_heat) with scaling = +1   (NOT -1 — see the SIGN TRAP above)
    --
    --   Graceful degradation: where u_heat has decayed into solver noise, EikonalVel returns the ZERO
    --   vector (normvel_util.h: |grad|>1e-15 else 0), so the Poisson gets no forcing there => phi is
    --   HARMONIC there => in a tube with Neumann side-walls a harmonic function is LINEAR along the
    --   tube => the direction is STILL AXIAL. (The shortcut has no such safety net at the distal tips.)
    ------------------------------------------------------------------
    local Xsrc = EikonalVel(u_heat)
    Xsrc:set_scaling(heatSign)          -- +1

    local poisEq = ConvectionDiffusionFV1("ca_cyt", "Outer")
    poisEq:set_diffusion(1.0)
    poisEq:set_vector_source(Xsrc)
    if reg > 0 then poisEq:set_reaction_rate(reg) end   -- 1e-6: decay length 1e3 >> L=5 (harmless),
                                                        -- but makes disconnected fragments non-singular

    local poisDir = DirichletBoundary()
    poisDir:add(0.0, "ca_cyt", somaSub)   -- phi = 0 at the soma (the distance anchor)

    local poisDisc = LSGFsimpleDomainDiscretization(approxSpace)
    poisDisc:add(poisEq)
    poisDisc:add(poisDir)
    poisDisc:set_LSF(lsf)
    poisDisc:set_Neumann0_on_if_for("ca_cyt")
    poisDisc:project_LSF()

    phi:set(0.0)
    -- CG + a SYMMETRIC smoother: this solve is SPD and BiCGStab broke down on it in T0 (see HeatSolverDesc).
    local solver2 = util.solver.CreateSolver(HeatSolverDesc(approxSpace, 1e-8,
                        util.GetParam("-heatKrylov", "cg", "Krylov for the Poisson recovery: cg|bicgstab"),
                        util.GetParam("-heatSmoother", "sgs", "GMG smoother for the recovery (sgs = symmetric, required by cg)")))
    local A2 = AssembledLinearOperator(poisDisc)
    local b2 = GridFunction(approxSpace); b2:set(0.0)
    poisDisc:adjust_solution(phi)
    poisDisc:assemble_linear(A2, b2)
    solver2:init(A2, phi)
    solver2:apply(phi, b2)

    gradiente_domainDisc = poisDisc   -- GLOBAL: the disc matching the RETURNED field
    print(string.format("[heat] step3 (Poisson recovery) done   ||phi|| = %.6g", VecNorm(phi)))
    print(string.format("finish: HEAT METHOD (3-step, 2 solves) in %d s", os.time() - t0))
    return phi
end

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
--  Function to obtain the direction of the tubuline inside the neuron
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
function GenerateImaginaryMolecule(ApproxSpace_lsf, gridfunction_ls, domain)
    print("Start: Stationary Difusion inside")

    -- 2. Load a Domain -----------------------------------------------------------------
    local dom = domain

    -- 3. Approximation Space: set up approximation space ---------------------------------
    local approxSpace = ApproximationSpace(dom)
    approxSpace:add_fct("ca_cyt", "Lagrange", 1) -- element 1 
    approxSpace:init_levels()
    approxSpace:init_surfaces()
    approxSpace:init_top_surface()
    ----------------------------------------------------------------------------------------------
    -- Gridfuncctions
    -----------------------------------------------------------------------------------------------
    -- 1. Grid functions for the solution
    local lsf = gridfunction_ls
    -- 2. Grid function for the solution of the solution (mol)
    local u = GridFunction(approxSpace) --- no necesito crear una matri<  --- yes
    u:set(0.0) -- set initial value to 0
    ------------------
    -- Curvature
    -------------------
    local curvature = GridFunction(ApproxSpace_lsf) -- level-set function
    Kapla_Func = Kappla_LS()
    Kapla_Func:compute_kappla(lsf, curvature)

    ----------------------------------------------------------------------------------------------
    -- PARTIAL DIFFERENTIAL EQUATIONS 
    ----------------------------------------------------------------------------------------------
    -------------------
    -- Linker of flux of the molecule in the top of the branchs to genrate the gradient.
    -------------------
    -- ============================================================================
    -- GROWTH-CONE DETECTION (curvature gate) -- AFFECTS ALL MODELS (defined here,
    -- shared by every Model_*.lua via GenerateImaginaryMolecule). LSInfluxTopBranch
    -- selects growth cones where |curvature| in [interval_min, interval_max] and
    -- seeds the interior-diffusion source there; its gradient sets growth + tubulin
    -- direction. So a WRONG gate -> wrong cones -> wrong interior field -> wrong growth.
    --
    -- CALIBRATION 2026-07-02 vs ground-truth SWC tips (sim=9, real=20). FINDING: the
    -- curvature intervals below are RESOLUTION-DEPENDENT and UNRELIABLE on coarse
    -- meshes -- tip and shaft curvature distributions OVERLAP:
    --   sim  numRef0: tip |H| median 6.8 vs shaft p95 21.8  -> best curvature F1 ~0.55
    --   real numRef0: tip |H| median 175 vs shaft p95 593   -> best curvature F1 ~0.35
    -- A fixed [min,max] therefore over/under-selects. Consequences:
    --   1) CALIBRATE the interval PER-MESH against ground-truth tips; do NOT trust
    --      the hard-coded values below across resolutions/geometries.
    --   2) Detection needs numRef>=2: at numRef0 the real neuron's thin dendrites
    --      FRAGMENT (interface breaks) so even a topological detector fails.
    -- Tooling + ground truth (Mac side):
    --   3D_evaluator/scripts/calibrate_growth_cones_20260702.py  (sweeps |H| threshold
    --      AND skeleton-endpoint spur-prune vs the SWC tips; reports F1 + value to use)
    --   3D_evaluator/reference_trees/{sim_neuron,real_pyramidal}_tips_*.csv  (truth)
    --   3D_evaluator/reports/RESEARCH_3D_growth_methods_20260702.md  (why + fixes)
    -- Mesh-stabler cross-check: skeleton-endpoint detector (sim F1 0.75 > curvature 0.55).
    -- ============================================================================
    flujo = LSInfluxTopBranch() --- this is the flux that will go inside only in the boundary (membrana) - Linker
    flujo:set_Curvature(GridFunctionNumberData(curvature, "ca_cyt")) -- curvature
    flujo:set_LevelSet_gradient(GridFunctionGradientData(lsf, "ca_cyt")) -- Gradiente del LS para obtener la direccion del flujo
    if dim == 2 then
        flujo:set_interval_min(3.7) -- 10 para que sea mas selectivo ? -- original 5.5
        flujo:set_interval_max(25) -- 18.5
        if numRefs <= 5 then
            flujo:set_interval_min(3.7) -- 5.5 is the previous value
            flujo:set_interval_max(25) -- 20
        else -- refinamiento 6
            flujo:set_interval_min(10) -- 5.5 is the previous value
            flujo:set_interval_max(20) -- 20
        end
    elseif dim == 3 then
        flujo:set_interval_min(20)
        flujo:set_interval_max(27)
    end
    flujo:set_magnitud_influx(-1.0) -- 1.0

    --- DIfusion: 
    DiffusionEq = ConvectionDiffusionFV1("ca_cyt", "Outer") -- DIFUSION PROCESS
    DiffusionEq:set_diffusion(1) 
    DiffusionEq:set_vector_source(flujo)

    dirichletBND = DirichletBoundary()
    dirichletBND:add(0.0, "ca_cyt", "Outflow,Bottom")

    --NeumannBND = NeumannBoundary("ca_cyt","fv1")
    --NeumannBND:add(0.0, "Inner","Outer")

    ----------------------------------------------------------------------------
    --- Domain Discretization
    ----------------------------------------------------------------------------
    gradiente_domainDisc = LSGFsimpleDomainDiscretization(approxSpace) -- Free surface: Using the simple extrapolation for the embedded interface
    flujo:set_domain_discretizacion(gradiente_domainDisc) -- domain discretization
    gradiente_domainDisc:add(dirichletBND)
    gradiente_domainDisc:add(DiffusionEq)

    gradiente_domainDisc:set_LSF(lsf) -- Set the Level-Set Function
    gradiente_domainDisc:set_Neumann0_on_if_for("ca_cyt") --- indicando not cruzar el dominio 
    gradiente_domainDisc:project_LSF()
    -----------------------------------------------------------------------------------------------
    -- SOLUTION
    -----------------------------------------------------------------------------------------------
    Interpolate(0.0, u, "ca_cyt", "Inner")
    Interpolate(0.0, u, "ca_cyt", "Outer")
    Interpolate(0.0, u, "ca_cyt", "Outflow,Bottom")
    print("Initial condition created ----------------------------------- \n\n")

    -- set up solver (using 'util/solver_util.lua')
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- SOLVER
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- set up solver (using 'util/solver_util.lua')
    solverDesc = {
        type = "bicgstab",
        precond = {
            type = "gmg", -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
            approxSpace = approxSpace,
            smoother = "ilu", -- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
            cycle = "V", -- gmg-cycle ["V", "F", "W"]
            preSmooth = 2, -- number presmoothing steps
            postSmooth = 2, -- number postsmoothing steps,
            rap = true, -- use the Galerkin formular in the hierarchy
            baseSolver = (dim == 3) and ((util.GetParamNumber("-robustInteriorSolve", 0, "20260705: 3D base solver for interior/smoothing GMG: 0=bounded BiCGStab+ILU (fast, can stall/abort on deformed geometry), 1=direct LU (robust, slower)") == 1) and "lu" or { type = "bicgstab", precond = { type = "ilu" }, convCheck = { type = "standard", iterations = 200, absolute = 1e-9, reduction = 1e-2, verbose = false } }) or "lu", -- 20260705: -robustInteriorSolve toggles direct LU base (was: 20260703 dim==3 bounded iterative base to avoid full ILUT(0) factorization; 2D keeps direct "lu")
            baselevel = numPreRefs -- the base level of the hierarchy
        },
        convCheck = {
            type = "standard",
            iterations = 150,
            absolute = 1e-5,
            reduction = 1e-12,
            verbose = true
        }
    }

    local solver = util.solver.CreateSolver(solverDesc)

    --local solver = LU() --- tbhis is the exact solver -- no esta paralelizado, debes usar agglomartingSolver

    print("\nsolving...")
    local A = AssembledLinearOperator(gradiente_domainDisc)
    local b = GridFunction(approxSpace)
    b:set(0.0)

    gradiente_domainDisc:adjust_solution(u)
    gradiente_domainDisc:assemble_linear(A, b)

    solver:init(A, u)
    solver:apply(u, b)
    print("Finish: Stationary Difusion inside")

    return u
end

--  Function to obtain the direction of the tubuline inside the neuron
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
function GenerateImaginaryMolecule3(ApproxSpace_lsf, gridfunction_ls, domain)
    print("Start: Stationary Difusion of the additional concentration to obtain the internal direction")

    -- 2. Load a Domain -----------------------------------------------------------------
    local dom = domain

    -- 3. Approximation Space: set up approximation space ---------------------------------
    local approxSpace = ApproximationSpace(dom)
    approxSpace:add_fct("ca_cyt", "Lagrange", 1) -- element 1 
    approxSpace:init_levels()
    approxSpace:init_surfaces()
    approxSpace:init_top_surface()
    ----------------------------------------------------------------------------------------------
    -- Gridfuncctions
    -----------------------------------------------------------------------------------------------
    -- 1. Grid functions for the solution
    local lsf = gridfunction_ls
    -- 2. Grid function for the solution of the solution (mol)
    local u = GridFunction(approxSpace) --- no necesito crear una matri<  --- yes
    u:set(0.0) -- set initial value to 0

    ----------------------------------------------------------------------------------------------
    -- PARTIAL DIFFERENTIAL EQUATIONS 
    ----------------------------------------------------------------------------------------------

    --- DIfusion: 
    DiffusionEq = ConvectionDiffusionFV1("ca_cyt", "Outer") -- DIFUSION PROCESS
    DiffusionEq:set_diffusion(1) 

    dirichletBND = DirichletBoundary()
    dirichletBND:add(0.0, "ca_cyt", "Outflow")
    dirichletBND:add(100.0, "ca_cyt", "Bottom")

    NeumannBND = NeumannBoundary("ca_cyt","fv1")
    NeumannBND:add(100.0, "Bottom","Outer")

    ----------------------------------------------------------------------------
    --- Domain Discretization
    ----------------------------------------------------------------------------
    domainDisc = LSGFsimpleDomainDiscretization(approxSpace) -- Free surface: Using the simple extrapolation for the embedded interface
    domainDisc:add(dirichletBND)
    domainDisc:add(DiffusionEq)
    domainDisc:add(NeumannBND)

    domainDisc:set_LSF(lsf) -- Set the Level-Set Function
    domainDisc:set_Neumann0_on_if_for("ca_cyt") --- indicando not cruzar el dominio 
    domainDisc:project_LSF()
    -----------------------------------------------------------------------------------------------
    -- SOLUTION
    -----------------------------------------------------------------------------------------------
    Interpolate(0.0, u, "ca_cyt", "Inner")
    Interpolate(0.0, u, "ca_cyt", "Outer")
    Interpolate(0.0, u, "ca_cyt", "Outflow")
    Interpolate(100.0, u, "ca_cyt", "Bottom") -- es como la interface
    print("Initial condition created ----------------------------------- \n\n")

    -- set up solver (using 'util/solver_util.lua')
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- SOLVER
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- set up solver (using 'util/solver_util.lua')
    solverDesc = {
        type = "bicgstab",
        precond = {
            type = "gmg", -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
            approxSpace = approxSpace,
            smoother = "ilu", -- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
            cycle = "V", -- gmg-cycle ["V", "F", "W"]
            preSmooth = 2, -- number presmoothing steps
            postSmooth = 2, -- number postsmoothing steps,
            rap = true, -- use the Galerkin formular in the hierarchy
            baseSolver = (dim == 3) and ((util.GetParamNumber("-robustInteriorSolve", 0, "20260705: 3D base solver for interior/smoothing GMG: 0=bounded BiCGStab+ILU (fast, can stall/abort on deformed geometry), 1=direct LU (robust, slower)") == 1) and "lu" or { type = "bicgstab", precond = { type = "ilu" }, convCheck = { type = "standard", iterations = 200, absolute = 1e-9, reduction = 1e-2, verbose = false } }) or "lu", -- 20260705: -robustInteriorSolve toggles direct LU base (was: 20260703 dim==3 bounded iterative base to avoid full ILUT(0) factorization; 2D keeps direct "lu")
            baselevel = numPreRefs -- the base level of the hierarchy
        },
        convCheck = {
            type = "standard",
            iterations = 150,
            absolute = 1e-5,
            reduction = 1e-12,
            verbose = true
        }
    }

    --local solver = util.solver.CreateSolver(solverDesc)

    local solver = LU() --- tbhis is the exact solver

    print("\nsolving...")
    local A = AssembledLinearOperator(domainDisc)
    local b = GridFunction(approxSpace)

    domainDisc:adjust_solution(u)
    domainDisc:assemble_linear(A, b)

    solver:init(A, u)
    solver:apply(u, b)
    print("Finish: Stationary Difusion of the additional concentration to obtain the internal direction")

    return u
end
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
--  Function to Smmooth of Curvature
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
function SmoothCutrvature(ApproxSpace_lsf, gridfunction_curvature, domain, Difusion_smooth)

    -- 2. My inputs-----------------------------------------------------------------
    local dom = domain
    local approxSpace = ApproxSpace_lsf
    local curvature = gridfunction_curvature
    local difusion = Difusion_smooth -- the value of the difussion coeficiente 

    ----------------------------------------------------------------------------
    --- Domain Discretization
    ----------------------------------------------------------------------------
    LSTransportUpwind = FullUpwind()
    SmoothCurvature_Imp = FV1_Convection(LSTransportUpwind, "ca_cyt", "Outer")
    SmoothCurvature_Imp:set_diffusion(difusion)
    SmoothCurvature_Imp:set_non_sink(false)

    --- 2. Domain Discretization del transporte del LS implicit method
    CurvatureDomainDisc = DomainDiscretization(approxSpace)
    CurvatureDomainDisc:add(SmoothCurvature_Imp)

    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    --- Solver for the implicit solution of the Transport disc: 
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    solverDesc2 = {
        type = "bicgstab", -- "linear", "bicgstab"
        precond = {
            type = "gmg", -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
            approxSpace = approxSpace,
            smoother = "ilu", -- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
            cycle = "V", -- gmg-cycle ["V", "F", "W"]
            preSmooth = 2, -- number presmoothing steps
            postSmooth = 2, -- number postsmoothing steps,
            rap = true, -- use the Galerkin formular in the hierarchy
            baseSolver = (dim == 3) and ((util.GetParamNumber("-robustInteriorSolve", 0, "20260705: 3D base solver for interior/smoothing GMG: 0=bounded BiCGStab+ILU (fast, can stall/abort on deformed geometry), 1=direct LU (robust, slower)") == 1) and "lu" or { type = "bicgstab", precond = { type = "ilu" }, convCheck = { type = "standard", iterations = 200, absolute = 1e-9, reduction = 1e-2, verbose = false } }) or "lu", -- 20260705: -robustInteriorSolve toggles direct LU base (was: 20260703 dim==3 bounded iterative base to avoid full ILUT(0) factorization; 2D keeps direct "lu")
            baselevel = numPreRefs -- the base level of the hierarchy
        },
        convCheck = {
            type = "standard",
            iterations = 100,
            absolute = 1e-8,
            reduction = 1e-12,
            verbose = true
        }
    }
    solver2 = util.solver.CreateSolver(solverDesc2)
    --local solver = LU() 

    print("\nsolving the smooth for the curvature...")
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- 4. Move the LSF : Fenomeno fisico para que se mueva 
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    local geo_dt = 0.01
    local geo_NumTimeSteps = 5

    util.SolveLinearTimeProblem(curvature, {
        domainDisc = CurvatureDomainDisc,
        reassemble = true
    }, solver2, nil, nil, "ImplEuler", 1, 0, nil, geo_dt, 0.001 , nil, false, nil, 0, geo_NumTimeSteps)

    return curvature
end

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
--  Function to Smmooth of Curvature
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
function SmoothCutrvatureAnisotropic(ApproxSpace_lsf, gridfunction_ls, gridfunction_curvature, domain , Difusion_smooth)

    -- 2. My inputs-----------------------------------------------------------------
    local dom = domain
    local approxSpace = ApproxSpace_lsf
    local lsf = gridfunction_ls
    local curvature = gridfunction_curvature
    local difusion = Difusion_smooth -- the value of the difussion coeficiente

    ----------------------------------------------------------------------------
    --- Domain Discretization
    ----------------------------------------------------------------------------
    TensorLink = LSTensorLinker()
    TensorLink:set_gradient_stationary_difussion(GridFunctionGradientData(lsf, "ca_cyt"))
    TensorLink:set_diffusion_coeff(difusion)

    local SmoothCurvature_Imp = ConvectionDiffusionFV1("ca_cyt", "Outer") -- DIFUSION PROCESS
    SmoothCurvature_Imp:set_diffusion(TensorLink) 

    --- 2. Domain Discretization del transporte del LS implicit method
    CurvatureDomainDisc = DomainDiscretization(approxSpace)
    --TensorLink:set_domain_discretizacion(CurvatureDomainDisc)
    CurvatureDomainDisc:add(SmoothCurvature_Imp)

    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    --- Solver for the implicit solution of the Transport disc: 
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    solverDesc2 = {
        type = "bicgstab", -- "linear", "bicgstab"
        precond = {
            type = "gmg", -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
            approxSpace = approxSpace,
            smoother = "ilu", -- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
            cycle = "V", -- gmg-cycle ["V", "F", "W"]
            preSmooth = 2, -- number presmoothing steps
            postSmooth = 2, -- number postsmoothing steps,
            rap = true, -- use the Galerkin formular in the hierarchy
            baseSolver = (dim == 3) and ((util.GetParamNumber("-robustInteriorSolve", 0, "20260705: 3D base solver for interior/smoothing GMG: 0=bounded BiCGStab+ILU (fast, can stall/abort on deformed geometry), 1=direct LU (robust, slower)") == 1) and "lu" or { type = "bicgstab", precond = { type = "ilu" }, convCheck = { type = "standard", iterations = 200, absolute = 1e-9, reduction = 1e-2, verbose = false } }) or "lu", -- 20260705: -robustInteriorSolve toggles direct LU base (was: 20260703 dim==3 bounded iterative base to avoid full ILUT(0) factorization; 2D keeps direct "lu")
            baselevel = numPreRefs -- the base level of the hierarchy
        },
        convCheck = {
            type = "standard",
            iterations = 100,
            absolute = 1e-8,
            reduction = 1e-12,
            verbose = true
        }
    }
    solver2 = util.solver.CreateSolver(solverDesc2)
    --local solver = LU() 

    print("\nsolving the smooth for the curvature...")
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- 4. Move the LSF : Fenomeno fisico para que se mueva 
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    local geo_dt = 0.01
    local geo_NumTimeSteps = 5
    local geo_dtmin = 0.001

    if numRefs == 7 then
        geo_dt = 0.001
        geo_NumTimeSteps = 50
        geo_dtmin = 0.0001
    end

    util.SolveLinearTimeProblem(curvature, {
        domainDisc = CurvatureDomainDisc,
        reassemble = true
    }, solver2, nil, nil, "ImplEuler", 1, 0, nil, geo_dt, geo_dtmin , nil, false, nil, 0, geo_NumTimeSteps)

    return curvature
end




-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
--  Function to Smmooth of the velocity of the interface:
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
function Smooth_Normal_Interface_Velocity(ApproxSpace_lsf, gridfunction_normal_velocity, domain, Difusion_smooth)

    -- 2. My inputs-----------------------------------------------------------------
    local dom = domain
    local approxSpace = ApproxSpace_lsf
    local InterfaceVelocity = gridfunction_normal_velocity
    local difusion = Difusion_smooth -- the value of the difussion coeficiente 

    ----------------------------------------------------------------------------
    --- Domain Discretization
    ----------------------------------------------------------------------------
    LSTransportUpwind = FullUpwind()
    SmoothVelNorm_Imp = FV1_Convection(LSTransportUpwind, "ca_cyt", "Outer")
    SmoothVelNorm_Imp:set_diffusion(difusion)
    SmoothVelNorm_Imp:set_non_sink(false)

    --- 2. Domain Discretization del transporte del LS implicit method
    InterfaceVelocityDomainDisc = DomainDiscretization(approxSpace)
    InterfaceVelocityDomainDisc:add(SmoothVelNorm_Imp)

    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    --- Solver for the implicit solution of the Transport disc: 
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    solverDesc2 = {
        type = "bicgstab", -- "linear", "bicgstab"
        precond = {
            type = "gmg", -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
            approxSpace = approxSpace,
            smoother = "ilu", -- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
            cycle = "V", -- gmg-cycle ["V", "F", "W"]
            preSmooth = 2, -- number presmoothing steps
            postSmooth = 2, -- number postsmoothing steps,
            rap = true, -- use the Galerkin formular in the hierarchy
            baseSolver = (dim == 3) and ((util.GetParamNumber("-robustInteriorSolve", 0, "20260705: 3D base solver for interior/smoothing GMG: 0=bounded BiCGStab+ILU (fast, can stall/abort on deformed geometry), 1=direct LU (robust, slower)") == 1) and "lu" or { type = "bicgstab", precond = { type = "ilu" }, convCheck = { type = "standard", iterations = 200, absolute = 1e-9, reduction = 1e-2, verbose = false } }) or "lu", -- 20260705: -robustInteriorSolve toggles direct LU base (was: 20260703 dim==3 bounded iterative base to avoid full ILUT(0) factorization; 2D keeps direct "lu")
            baselevel = numPreRefs -- the base level of the hierarchy
        },
        convCheck = {
            type = "standard",
            iterations = 100,
            absolute = 1e-8,
            reduction = 1e-12,
            verbose = true
        }
    }
    local solver2 = util.solver.CreateSolver(solverDesc2)

    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- 4. Move the LSF : Fenomeno fisico para que se mueva 
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    local geo_dt = 0.01
    local geo_NumTimeSteps = 5

    util.SolveLinearTimeProblem(InterfaceVelocity, {
        domainDisc = InterfaceVelocityDomainDisc,
        reassemble = true
    }, solver2, nil, nil, "ImplEuler", 1, 0, nil, geo_dt, 0.001, nil, false, nil, 0, geo_NumTimeSteps)
    print("\nsolving the smooth ..\n\n")
    return InterfaceVelocity
end









--- Compute the initial inhibition values outside the neuron. It will avoid start from 0
-- @param num_Steps_Initial_inh number of steps to compute the initial inhibition values
-- @param gridfunction_inhibition GridFunction to store the inhibition values
-- @return gridfunction_inhibition
function GenerateInitialInhibition (num_Steps_Initial_inh, gridfunction_inhibition)
    local dt = 0.015
    local gridfunction_inhibition = gridfunction_inhibition
    print("Initial Values: Compute Solution of the model outside the neuron ....")

    --[[
    local solverDesc = {
        type = "bicgstab",
        precond = {
            type = "gmg", -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
            approxSpace = approxSpace3,
            smoother = "ilu", -- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
            cycle = "V", -- gmg-cycle ["V", "F", "W"]
            preSmooth = 2, -- number presmoothing steps
            postSmooth = 2, -- number postsmoothing steps,
            rap = true, -- use the Galerkin formular in the hierarchy
            baseSolver = (dim == 3) and ((util.GetParamNumber("-robustInteriorSolve", 0, "20260705: 3D base solver for interior/smoothing GMG: 0=bounded BiCGStab+ILU (fast, can stall/abort on deformed geometry), 1=direct LU (robust, slower)") == 1) and "lu" or { type = "bicgstab", precond = { type = "ilu" }, convCheck = { type = "standard", iterations = 200, absolute = 1e-9, reduction = 1e-2, verbose = false } }) or "lu", -- 20260705: -robustInteriorSolve toggles direct LU base (was: 20260703 dim==3 bounded iterative base to avoid full ILUT(0) factorization; 2D keeps direct "lu")
            baselevel = numPreRefs -- the base level of the hierarchy
        },
        convCheck = {
            type = "standard",
            iterations = 150,
            absolute = 1e-5,
            reduction = 1e-12,
            verbose = true
        }
    }

    local solver = util.solver.CreateSolver(solverDesc)

    --local solver = LU() --- tbhis is the exact solver -- no esta paralelizado, debes usar agglomartingSolver

    print("\nsolving...")
    local A = AssembledLinearOperator(domainDisc_outside)
    local b = GridFunction(approxSpace3)
    b:set(0.0)

    domainDisc_outside:adjust_solution(gridfunction_inhibition)
    domainDisc_outside:assemble_linear(A, b)

    solver:init(A, gridfunction_inhibition)
    solver:apply(gridfunction_inhibition, b)
    --]]
    for step = 1, num_Steps_Initial_inh do 
        timeDisc_Concentration2:prepare_step(solTimeSeries2, dt)
        newtonSolver2:prepare(gridfunction_inhibition)
        newtonSolver2:apply(gridfunction_inhibition)

        -- update new time
        time2 = timeDisc_Concentration2:future_time()
        
        -- update time series
        local oldestSol2 = solTimeSeries2:oldest() -- obtain  the reference the oldest solution
        VecAssign(oldestSol2, gridfunction_inhibition) -- ssign the neew to the oldes - reinitialize the oldest solution with the new values
        solTimeSeries2:push_discard_oldest(oldestSol2, time2) -- push the oldest solition swith the new values  to the front, oldest sol pointer is poped  fromn end
        
    end

    print("Calculated the initial value: outside the neuron ....")
    return gridfunction_inhibition
end











---- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
--  Function to obtain the intensity of the velocity to make it divergence free (gradA * V = 0)
--  grad (V) is not 0 , I will add A.
--  Then I will use this to have : grad (A * V) = 0 
--   Debo cumplir : A * grad V + V * grad A = 0 -> A * grad V = - V * grad A

-- dA/dt + V * grad A = 0 -> con esto aseguro que la divergencia sea 0 : grad(A*V)=0
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
function Divergence_Free_intensity_Calculation(ApproxSpace_lsf, gridfunction_intensity, GradientDataVelocity, gridfunction_ls, domain,label , name)
    print("Start: Intensity calculation")
    
    -- 1. Parameters ----------------------------------------------------------------- 
    if numRefs == 6 or numRefs == 5 then
        geo_dt = 0.01 -- time step length
        geo_NumTimeSteps = 700  
    end

    if numRefs == 4 then
        geo_dt = 0.1 -- time step length
        geo_NumTimeSteps = 500 -- Dependiendo de la geometria debes de modificar este valor para asegurarte que la advencion sea suficiente para que vaya por toda la geometria       
    end

    -- 2. Load a Domain -----------------------------------------------------------------
    local dom = domain
    -- 3. Approximation Space: set up approximation space ---------------------------------
    local approxSpace = ApproxSpace_lsf
    ----------------------------------------------------------------------------------------------
    -- Gridfuncctions
    -----------------------------------------------------------------------------------------------
    -- 1. Grid functions of the interface
    local lsf = gridfunction_ls
    -- 2. Grid function for the solution of the solution (mol)
    local u = gridfunction_intensity
    
    ----------------------------------------------------------------------------------------------
    -- PARTIAL DIFFERENTIAL EQUATIONS 
    ----------------------------------------------------------------------------------------------
    local divergenceEq = FV1_Convection(FullUpwind(), "ca_cyt", "Outer")
    divergenceEq:set_diffusion(0.0)
    divergenceEq:set_source(0.0)
    divergenceEq:set_non_sink(label)


    local direction = LSTubulinVelocity()
    direction:set_gradient_stationary_difussion(GradientDataVelocity) -- Gradiente del LS para obtener la direccion del flujo
    direction:set_magnitud(1.0)
    
    function ourVelocityField2d(x, y, t)
        return	0, 1.0
    end
    --velocityField = LuaUserVector("ourVelocityField2d")
    divergenceEq:set_velocity(direction)
            
    ----------------------------------------------------------------------------
    local geo_dirichletBND = DirichletBoundary()
    geo_dirichletBND:add(1, "ca_cyt", "Bottom") -- subset : Side_Cyt

    --local NeumannBND_Soma =  NeumannBoundary("t","fv1")
    --NeumannBND_Soma:add(-1, "Bottom","Outer") -- positivo quita y negativo agrega

    ----------------------------------------------------------------------------
    --- Domain Discretization
    ----------------------------------------------------------------------------
    intensity_domainDisc = LSGFsimpleDomainDiscretization(approxSpace) -- Free surface: Using the simple extrapolation for the embedded interface
    direction:set_domain_discretizacion(intensity_domainDisc)
    intensity_domainDisc:add(divergenceEq)
    --intensity_domainDisc:add(convec)
    intensity_domainDisc:add(geo_dirichletBND)

    intensity_domainDisc:set_LSF(lsf) -- Set the Level-Set Function
    intensity_domainDisc:set_Neumann0_on_if_for("ca_cyt") --- indicando not cruzar el dominio 
    --intensity_domainDisc:set_Dirichlet_on_if_for("ca_cyt",0.0)
    intensity_domainDisc:project_LSF()

    -----------------------------------------------------------------------------------------------
    -- SOLUTION
    -----------------------------------------------------------------------------------------------
    Interpolate(0.0, u, "ca_cyt", "Outer,Inner,Outflow") 
    
    print("Initial condition created ----------------------------------- \n\n")

    -- set up solver (using 'util/solver_util.lua')
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    -- SOLVER
    -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
    --util.solver.defaults.approxSpace = approxSpace --- yes

    local geo_solverDesc = { --- yes
        type = "newton", -- nombre del metodo de discretization 
    
        lineSearch = {
            type = "standard", -- ["standard", "none"]
            maxSteps = 8, -- maximum number of line search steps
            lambdaStart = 1, -- start value for scaling parameter
            lambdaReduce = 0.5, -- reduction factor for scaling parameter
            acceptBest = true, -- check for best solution if true
            checkAll = false -- check all maxSteps steps if true 
        },
    
        convCheck = {
            type = "standard",
            iterations = 50, -- maximum number of iterations
            absolute = 5e-7, -- absolut value of defect to be reached; usually 1e-7 - 1e-9
            reduction = 1e-10, -- reduction factor of defect to be reached; usually 1e-6 - 1e-8
            verbose = true -- print convergence rates if true
        },
    
        linSolver = {
            type = "bicgstab",
            precond = {
                type = "gmg", -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
                approxSpace = approxSpace,
                smoother = "ilu", -- pre- and postsmoother, only for gmg ["ilu", "ilut", "jac", "gs", "sgs"]
                cycle = "V", -- gmg-cycle ["V", "F", "W"]
                preSmooth = 2, -- number presmoothing steps
                postSmooth = 2, -- number postsmoothing steps
                rap = true, -- comutes RAP-product instead of assembling if true
                baseLevel = numPreRefs, -- gmg - coarsest level
                baseSolver = "lu"
            },
    
            convCheck = {
                type = "standard",
                iterations = 150, -- number of iterations
                absolute = 1e-8, -- absolut value of defact to be reached; usually 1e-8 - 1e-10 (may not be larger than in newton section)
                reduction = 1e-12, -- reduction factor of defect to be reached
                verbose = true -- print convergence rates if true
            }
        }
    } -- definir el error para que ecuaciones quieres ... 


    local geo_newtonSolver = util.solver.CreateSolver(geo_solverDesc) --- yes

    --------------------------------
    --  Setup Time Discretization --
    --------------------------------
    local timeDisc_Intensity = util.CreateTimeDisc(intensity_domainDisc, "ImplEuler", 1)
    timeDisc_Intensity:set_stage(1)

    --------------------
    --  Algebra setup --
    --------------------
    geo_newtonSolver:init(AssembledOperator(timeDisc_Intensity))

    -------------------
    -- Time stepping --
    -------------------
    local time = 0
    -- create new grid function for old value
    local uold = u:clone() -- create a copy of the solution
    local solTimeSeries = SolutionTimeSeries()
    solTimeSeries:push(uold, time)

    --------------------------------------------------------------------------------
    --	Output facilities
    --------------------------------------------------------------------------------
    local geo_vtkOut = VTKOutput();
    geo_vtkOut:clear_selection();
    geo_vtkOut:select_nodal("ca_cyt", "intensity"); -- distance function es la variable "distFunc"
    

    CreateDirectory(name)
    diver_folder = name .."/" -- You will use out:print(folder .. "file", u, time, step)


    geo_vtkOut:print(diver_folder .. "evaluation", u, 0, 0) --- antes era time
    out2:print(diver_folder .. "LSF", lsf, 0, 0) -- plot the initial condition of inferface

        -- print velocity: 
    local Vgeo_vtkOut = VTKOutput();
    Vgeo_vtkOut:clear_selection();
    Vgeo_vtkOut:select(direction, "direc_divergence"); -- distance function es la variable "distFunc"
    Vgeo_vtkOut:print(diver_folder .. "direc", u, 0, 0) --- antes era time
    ---------------------------------c-----------------------------------------------
    --  Timestepping
    --------------------------------------------------------------------------------
    for step = 1, geo_NumTimeSteps do --- 

        print("Time step " .. step .. " at time " .. time)

        timeDisc_Intensity:prepare_step(solTimeSeries, geo_dt)
        geo_newtonSolver:prepare(u)
        geo_newtonSolver:apply(u)

        -- update new time
        time = timeDisc_Intensity:future_time()
        geo_vtkOut:print(diver_folder .. "evaluation", u, step, time) --- antes era time
        Vgeo_vtkOut:print(folder .. "direc", u, step, time) --- antes era time
        --out2:print(diver_folder .. "LSF", lsf, step, time) -- plot the initial condition of inferface   


        -- update time series
        local oldestSol = solTimeSeries:oldest() -- obtain  the reference the oldest solution
        VecAssign(oldestSol, u) -- ssign the neew to the oldes - reinitialize the oldest solution with the new values
        solTimeSeries:push_discard_oldest(oldestSol, time) -- push the oldest solition swith the new values  to the front, oldest sol pointer is poped  fromn end
        --  intensity_domainDisc:project_LSF()
    end

    print("finish: Intensity of the velocity to make it divergence free")
    return u
end




--- Submit a SLURM post-processing job to pack (merge + optionally delete)
--  a range of UG4 VTK outputs.
--
--  This is meant to be called from the simulation (preferably only on rank 0)
--  every N timesteps, to avoid accumulating millions of partitioned files.
--
--  The SLURM script receives:
--    ROOT : case root directory (outside "merged/")
--    T0   : first timestep index to process (inclusive)
--    T1   : last timestep index to process  (inclusive)
--    SNAPSHOT_INTERVAL : output stride used by UG4 VTK printing
--
--  Example:
--    Submit_pack_job(0, 49)   -- process Output_t000000 ... Output_t000049
--
-- @param t0 integer First timestep index to process (inclusive)
-- @param t1 integer Last timestep index to process (inclusive)
-- @param post_slurm_script string Path to the SLURM script that will do the packing
-- @param root string Root directory of the case
-- @return integer exit_code Return code from os.execute (0 typically means success)
function Submit_pack_job(t0, t1, snapshot_interval, post_slurm_script, root)
    if ProcRank and ProcRank() ~= 0 then -- Only submit from rank 0
        return 0
    end

    local cmd = string.format(
        'sbatch --export=ALL,ROOT="%s",T0=%d,T1=%d,SNAPSHOT_INTERVAL=%d "%s"',
        root, t0, t1, snapshot_interval, post_slurm_script
    )
    print("[sbatch] "..cmd)
    return os.execute(cmd)
end
