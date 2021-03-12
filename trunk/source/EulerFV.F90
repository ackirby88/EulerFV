Program EulerFV
    use inputs_module, only: timeScheme,timeSteps
    use globals_module, only: numFields,numElem
    use inputs_initialization_module
    use read_mesh_module
    use geometry_module
    use solution_initialization_module
    use evolve_solution_module
    implicit none

    call inputs_initialization()
    call read_mesh()
    call geometry_initialization()
    call solution_initialization(numFields,numElem)
    call evolve_solution(timeScheme,timeSteps)
end program