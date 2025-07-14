

------------------------------------------------------------------------------------------------------------
-- Function to generate the initial condition of the level-set function
------------------------------------------------------------------------------------------------------------

--- Generate initial condition for the LSF from the subset skeletum (called inner)
-- Use the function FV1_Convection() and eikanol velocity from pluglin LevelSet; and DirichletBoundary is the radius of the position of the interface.
-- @param ApproxSpace_lsf ApproxSpace of the level-set function 
-- @param gridfunction_ls Gridfunction of the level-set function (initial condition for this simulation is 0):
-- @return lsf GridFunction after the simulation
function Initialgeometry (ApproxSpace_lsf, gridfunction_ls)
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
        geo_dirichletBND:add(-0.1, "ca_cyt", "Inner") -- subset : Side_Cyt
         --- el valor siempre ha sido -0.1
    else
        geo_dirichletBND:add(-0.06, "ca_cyt", "Inner") -- subset : Side_Cyt
    end
    
   
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

    --baseSolver=Base
    -- usar el exactSolver para paralelo en 3D
    baseSolver=exactSolver
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
    BiCGStabSolver:set_preconditioner(gmg)
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

    Interpolate(-0.1, lsf, "ca_cyt", "Inner")
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
function GenerateImaginaryMolecule2(ApproxSpace_lsf, gridfunction_gradiente, gridfunction_ls, domain, inicial_pregunta)
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

    if numRefs == 4 or numRefs == 3 then
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
    geo_elemDisc_inner:set_diffusion(0.05) -- this will smooth the eikonal to make it parallel
    geo_elemDisc_inner:set_non_sink(false)

    if inicial then 
        geo_elemDisc_inner:set_diffusion(0.1) -- estaba en 0.06
    end

    if dim == 3 and inicial then
        geo_elemDisc_inner:set_diffusion(0.1) -- this will smooth the eikonal to make it parallel
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
            absolute = 5e-5, -- absolut value of defect to be reached; usually 1e-7 - 1e-9
            reduction = 1e-8, -- reduction factor of defect to be reached; usually 1e-6 - 1e-8
            verbose = false -- print convergence rates if true
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
                baseSolver = "lu"
            },
    
            convCheck = {
                type = "standard",
                iterations = 150, -- number of iterations
                absolute = 1e-8, -- absolut value of defact to be reached; usually 1e-8 - 1e-10 (may not be larger than in newton section)
                reduction = 1e-12, -- reduction factor of defect to be reached
                verbose = false -- print convergence rates if true
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
        verbose = false -- print convergence rates if true
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
            verbose = false -- print convergence rates if true
        }
    }
    } -- definir el error para que ecuaciones quieres ... 
    if dim == 2 then
        geo_solverDesc = geo_solverDesc2D
    elseif dim == 3 then
        geo_solverDesc = geo_solverDesc3D
    end
    local geo_newtonSolver = util.solver.CreateSolver(geo_solverDesc) --- yes

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
    util.SolveNonlinearTimeProblem(u, gradiente_domainDisc, geo_newtonSolver, geo_vtkOut, geo_name , "ImplEuler", 1, 0, geo_endTime, geo_dt, dt_min, 0.5); 

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
            baseSolver = "lu", -- direct coarse grid solver
            baselevel = numPreRefs -- the base level of the hierarchy
        },
        convCheck = {
            type = "standard",
            iterations = 150,
            absolute = 1e-5,
            reduction = 1e-12,
            verbose = false
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
            baseSolver = "lu", -- direct coarse grid solver
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
            baseSolver = "lu", -- direct coarse grid solver
            baselevel = numPreRefs -- the base level of the hierarchy
        },
        convCheck = {
            type = "standard",
            iterations = 100,
            absolute = 1e-8,
            reduction = 1e-12,
            verbose = false
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
function SmoothCutrvatureAnisotropic(ApproxSpace_lsf, gridfunction_curvature, domain , Difusion_smooth)

    -- 2. My inputs-----------------------------------------------------------------
    local dom = domain
    local approxSpace = ApproxSpace_lsf
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
            baseSolver = "lu", -- direct coarse grid solver
            baselevel = numPreRefs -- the base level of the hierarchy
        },
        convCheck = {
            type = "standard",
            iterations = 100,
            absolute = 1e-8,
            reduction = 1e-12,
            verbose = false
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
            baseSolver = "lu", -- direct coarse grid solver
            baselevel = numPreRefs -- the base level of the hierarchy
        },
        convCheck = {
            type = "standard",
            iterations = 100,
            absolute = 1e-8,
            reduction = 1e-12,
            verbose = false
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
            baseSolver = "lu", -- direct coarse grid solver
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
            verbose = false -- print convergence rates if true
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
                verbose = false -- print convergence rates if true
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