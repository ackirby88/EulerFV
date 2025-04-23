module grid_data_type_module
    use my_kinddefs
    implicit none

    public :: node_type
    public :: elem_type
    public :: edge_type

    type node_type
        real(dp) :: x ! x-coordinate
        real(dp) :: y ! y-coordinate
    end type node_type

    type elem_type
        real(dp)             :: area,perimeter,oneOarea
        real(dp)             :: centroids(2) ! (x,y) coordinate
        integer(i4), pointer :: vtxList(:)   ! list of vertices
    end type elem_type

    type edge_type
        integer(i4) :: n1, n2         ! start/end node id
        integer(i4) :: e1, e2         ! left/right element id
        integer(i4) :: edgeType       ! interior/boundary edge
        real(dp)    :: C2Mvector(2,2) ! (x/y,Left/Right)
    end type edge_type
 end module grid_data_type_module

 module globals_module
    use my_kinddefs
    use grid_data_type_module
    implicit none

    !  Solution Data
    integer(i4),parameter :: numFields = 4, dim = 2

    integer(i4) :: alloc_stat
    integer(i4) :: printCount = 0
    integer(i4) :: printCountRes = 0
    integer(i4) :: printCountCoeffs = 0

    ! Node data
    integer(i4)              :: numNode     ! number of nodes
    type(node_type), pointer :: nodeList(:) ! array of nodes

    ! Element data
    integer(i4)              :: numElem     ! number of elements
    real(dp),        pointer :: elem_h(:)   ! array of element size parameter
    type(elem_type), pointer :: elemList(:) ! array of elements

    ! Edge data
    integer(i4)              :: numEdge             ! number of edges
    integer(i4)              :: numInterior         ! number of interior edges
    integer(i4)              :: numBoundary         ! number of boundary edges
    type(edge_type), pointer :: edgeList(:)         ! array of edges
    integer(i4),     pointer :: interiorEdgeList(:) ! array of interior edges
    integer(i4),     pointer :: boundaryEdgeList(:) ! array of boundary edges
    real(dp),        pointer :: edgeNormals(:,:)    ! (x/y,numEdge): edge normal vectors
end module globals_module

module inputs_module
    use my_kinddefs
    implicit none

    character(80) :: input_name      ! Name of the input file

    character(80) :: stats_name      ! Name of the mesh stats file
    character(80) :: nodes_name      ! Name of the mesh nodes file
    character(80) :: cells_name      ! Name of the mesh cells file
    character(80) :: edges_name      ! Name of the mesh edges file

    character(80) :: checkpoint_name ! Name of checkpoint file to initialize solution
    character(20) :: output_file     ! Name of the output file
    integer(i4)   :: output_freq     ! Output file frequency
    integer(i4)   :: checkpoint_freq ! Checkpoint file frequency
    integer(i4)   :: logging_freq=10 ! Logging frequency (L2 norm)
    logical       :: restart         ! Restart solution flag from checkpoint

    character(20) :: limiterScheme   ! Name of limiter

    integer(i4)   :: spatial_order   ! Spatial order of accuracy
    integer(i4)   :: outer_bc        ! In/Out boundary condition
    integer(i4)   :: wall_bc         ! Wall boundary condition
    integer(i4)   :: fluxScheme      ! Flux scheme
    integer(i4)   :: timeScheme      ! Temporal discretization
    integer(i4)   :: timeSteps       ! Number of time steps
    logical       :: tacc            ! Time accurate flag

    real(dp)      :: gamma,gm1       ! Ratio of specific heats for air(1.4), gamma-1
    real(dp)      :: oneOgm1,gOgm1   ! 1/gm1, gamma/gm1
    real(dp)      :: fmach           ! Free-stream Mach number
    real(dp)      :: alpha           ! Angle of attack in degrees
    real(dp)      :: CFL             ! CFL number

    real(dp)      :: pressure0       ! Free-stream pressure
    real(dp)      :: density0        ! Free-stream density
    real(dp)      :: rho0            ! Free-stream density
    real(dp)      :: rhou0           ! Free-stream x-momentum
    real(dp)      :: rhov0           ! Free-stream y-momentum
    real(dp)      :: rhoe0           ! Free-stream energy density
    real(dp)      :: u0,v0           ! Free-stream velocities
    real(dp)      :: C0              ! Free-stream speed of sound
end module inputs_module

module solution_module
    use my_kinddefs
    implicit none

    real(dp),pointer :: Q(:,:)    ! (numFields,numElem): solution array
    real(dp),pointer :: R(:,:)    ! (numFields,numElem): residual array
    real(dp),pointer :: dQ(:,:,:) ! (numFields,x/y,numElem): solution gradient array
    real(dp),pointer :: dt(:)     ! (numElem): time step size
    contains

    subroutine solution_allocation(nf,nelem)
        use my_kinddefs
        implicit none

        integer(i4),intent(in) :: nf
        integer(i4),intent(in) :: nelem

        allocate(Q(nf,nelem))
        allocate(R(nf,nelem))
        allocate(dQ(nf,2,nelem))
        allocate(dt(nelem))
    end subroutine

    subroutine solution_deallocation()
        use my_kinddefs
        implicit none

        deallocate(Q)
        deallocate(R)
        deallocate(dQ)
        deallocate(dt)
    end subroutine
end module