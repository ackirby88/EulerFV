module read_mesh_module
    implicit none
    contains

    subroutine read_mesh()
        use my_kinddefs
        implicit none

        call read_mesh_stats()  !Get numNode, numElem, numEdge
        call read_nodes()       !Get x and y coordinates of the nodes
        call read_cells()       !Get the nodes of the elements
        call read_edges()       !Get n1,n2,e1,e2: the associated nodes and elements with each edge
    end subroutine

    subroutine read_mesh_stats()
        use my_kinddefs
        use inputs_module,  only: stats_name
        use globals_module, only: numNode,numElem,numEdge
        implicit none

        integer(i4) :: funit

        funit = 22
        open(UNIT = funit,FILE = stats_name,STATUS = 'old')
        read(funit,*) numNode, numElem, numEdge
        close(funit)

        call mesh_allocate()
    end subroutine read_mesh_stats

    subroutine mesh_allocate()
        use globals_module
        implicit none

        allocate(nodeList(numNode))
        allocate(elemList(numElem))
        allocate(edgeList(numEdge))
    end subroutine mesh_allocate

    subroutine read_nodes()
        use my_kinddefs
        use inputs_module,  only: nodes_name
        use globals_module, only: numNode,nodeList
        implicit none

        integer(i4) :: funit,i

        funit = 22
        open(UNIT = funit,FILE = nodes_name,STATUS = 'old')

        do i = 1,numNode
            read(funit,*) nodeList(i)%x, &  ! x-coordinate
                          nodeList(i)%y     ! y-coordinate
        end do
        close(funit)
    end subroutine read_nodes

    subroutine read_cells()
        use my_kinddefs
        use inputs_module,  only: cells_name
        use globals_module, only: numElem,elemList
        implicit none

        integer(i4) :: funit,i

        funit = 22
        open(UNIT = funit,FILE = cells_name,STATUS = 'old')

        ! read in the element node list (triangle: 3 per volume)
        do i = 1,numElem
           allocate(elemList(i)%vtxList(3))

           read(funit,*) elemList(i)%vtxList(1), &  ! node 1
                         elemList(i)%vtxList(2), &  ! node 2
                         elemList(i)%vtxList(3)     ! node 3
        end do
        close(funit)
    end subroutine read_cells

    subroutine read_edges()
        use my_kinddefs
        use inputs_module,  only: edges_name
        use globals_module, only: numEdge,numInterior,numBoundary
        use globals_module, only: edgeList,interiorEdgeList,boundaryEdgeList
        implicit none

        integer(i4) :: funit,i,inCount,extCount

        funit = 22
        open(UNIT = funit,FILE = edges_name,STATUS = 'old')

        numInterior = 0
        extCount = 0

       !Get n1,n2,e1,e2
        do i = 1,numEdge
            read(funit,*) edgeList(i)%n1,edgeList(i)%n2, &  ! start/end nodes
                          edgeList(i)%e1,edgeList(i)%e2, &  ! left/right elements
                          edgeList(i)%edgeType              ! edge type: interior(0), boundary(1,9)

            ! update number of interior edges
            if(edgeList(i)%edgeType == 0) numInterior = numInterior + 1
        end do
        close(funit)

        ! number of boundary edges
        numBoundary = numEdge - numInterior

        allocate(interiorEdgeList(numInterior))
        allocate(boundaryEdgeList(numBoundary))

        inCount = 1
        extCount = 1
        do i = 1,numEdge
            select case(edgeList(i)%edgeType)
                case(0) !Interior edge
                    interiorEdgeList(inCount) = i
                    inCount = inCount + 1

                case(1,9) !Boundary edge
                    boundaryEdgeList(extCount) = i
                    extCount = extCount + 1

                case default
                    print*,'Edge Type does not match any boundary conditions(0,1,9): type=',edgeList(i)%edgeType,' Stopping.'
                    stop
            end select
        end do
    end subroutine read_edges

end module read_mesh_module