#define _X_ 1
#define _Y_ 2

Module geometry_module
    implicit none
    contains

    subroutine geometry_initialization()
        use my_kinddefs
        use globals_module
        implicit none

        ! allocate element h size
        allocate(elem_h(numElem))

        call geometry_cell_area(numNode,numElem,nodeList,elemList,elem_h)
        call geometry_cell_centroid(numNode,numElem,nodeList,elemList)
        call geometry_cell_centroid_midpt_vector(numNode,numElem,numEdge,numInterior,numBoundary,&
                                                 nodeList,elemList,edgeList,interiorEdgeList,boundaryEdgeList)
        call geometry_edge_normals(numNode,numEdge,nodeList,edgeList,edgeNormals)
    end subroutine

    subroutine geometry_cell_area(numNode,numElem,nodeList,elemList,elem_h)
        use my_kinddefs
        use grid_data_type_module, only: node_type,elem_type
        implicit none

        integer(i4),       intent(in) :: numNode
        integer(i4),       intent(in) :: numElem
        type(node_type),   intent(in) :: nodeList(numNode)
        type(elem_type),intent(inout) :: elemList(numElem)
        real(dp),         intent(out) :: elem_h(numElem)

        integer(i4) :: triNumber
        integer(i4) :: v1,v2,v3
        real(dp)    :: a1,a2,b1,b2
        real(dp)    :: L1,L2,L3
        real(dp)    :: crossProduct

        ! Triangle Orientation:
        !           L2
        !   V1 * - - - - - * V3
        !       \         /
        !        \       /
        !      L3 \     / L1
        !          \   /
        !           \ /
        !            *
        !           V2


        do triNumber = 1,numElem
            ! get the 3 nodes of the element
            v1 = elemList(triNumber)%vtxList(1)
            v2 = elemList(triNumber)%vtxList(2)
            v3 = elemList(triNumber)%vtxList(3)

            ! calculate the edges lengths
            L1 = sqrt((nodeList(v3)%x - nodeList(v2)%x)**2 + (nodeList(v3)%y - nodeList(v2)%y)**2)
            L2 = sqrt((nodeList(v1)%x - nodeList(v3)%x)**2 + (nodeList(v1)%y - nodeList(v3)%y)**2)
            L3 = sqrt((nodeList(v2)%x - nodeList(v1)%x)**2 + (nodeList(v2)%y - nodeList(v1)%y)**2)

            !compute cross product to ensure orientation
            a1 = nodeList(v3)%x - nodeList(v2)%x ! vector a: x-coordinate
            a2 = nodeList(v3)%y - nodeList(v2)%y ! vector a: y-coordinate
            b1 = nodeList(v1)%x - nodeList(v2)%x ! vector b: x-coordinate
            b2 = nodeList(v1)%y - nodeList(v2)%y ! vector b: y-coordinate

            crossProduct = a1*b2 - a2*b1

            ! switch since it is oriented the wrong way
            if (crossProduct < 0) then
                print*,'cross product was negative. Quitting.'
                stop
            end if

            ! triangle area = 0.5 * crossProduct
            elemList(triNumber)%area  = 0.5 * crossProduct
            elemList(triNumber)%oneOarea = 1.0/elemList(triNumber)%area
            elemList(triNumber)%perimeter = l1 + l2 + l3

            elem_h(triNumber) = 1.0 / elemList(triNumber)%perimeter
        end do
    end subroutine

    subroutine geometry_cell_centroid(numNode,numElem,nodeList,elemList)
        use my_kinddefs
        use grid_data_type_module, only: node_type,elem_type
        implicit none

        integer(i4),       intent(in) :: numNode
        integer(i4),       intent(in) :: numElem
        type(node_type),   intent(in) :: nodeList(numNode)
        type(elem_type),intent(inout) :: elemList(numElem)


        integer(i4) :: v1,v2,v3
        integer(i4) :: i

        do i = 1,numElem
            v1 = elemList(i)%vtxList(1)
            v2 = elemList(i)%vtxList(2)
            v3 = elemList(i)%vtxList(3)

            elemList(i)%centroids(_X_) = onethird * (nodeList(v1)%x + nodeList(v2)%x + nodeList(v3)%x)
            elemList(i)%centroids(_Y_) = onethird * (nodeList(v1)%y + nodeList(v2)%y + nodeList(v3)%y)
        end do
    end subroutine

    subroutine geometry_cell_centroid_midpt_vector(numNode,numElem,numEdge,numInterior,numBoundary,&
                                                   nodeList,elemList,edgeList,interiorEdgeList,boundaryEdgeList)
        use my_kinddefs
        use grid_data_type_module, only: node_type,elem_type,edge_type
        implicit none

        integer(i4),       intent(in) :: numNode
        integer(i4),       intent(in) :: numElem
        integer(i4),       intent(in) :: numEdge
        integer(i4),       intent(in) :: numInterior
        integer(i4),       intent(in) :: numBoundary
        type(node_type),   intent(in) :: nodeList(numNode)
        type(elem_type),   intent(in) :: elemList(numElem)
        type(edge_type),intent(inout) :: edgeList(numEdge)
        integer(i4),       intent(in) :: interiorEdgeList(numInterior)
        integer(i4),       intent(in) :: boundaryEdgeList(numBoundary)

        real(dp)    :: midPoint(2)
        integer(i4) :: leftTri,rightTri
        integer(i4) :: edge_id
        integer(i4) :: i

        do i = 1,numInterior
            edge_id = interiorEdgeList(i)
            edgeList(edge_id)%C2Mvector(:,2) = 0.0

            leftTri  = edgeList(edge_id)%e1
            rightTri = edgeList(edge_id)%e2

            ! midpoint of edge
            midPoint(_X_) = half * (nodeList(edgeList(edge_id)%n2)%x + nodeList(edgeList(edge_id)%n1)%x) ! x-coordinate
            midPoint(_Y_) = half * (nodeList(edgeList(edge_id)%n2)%y + nodeList(edgeList(edge_id)%n1)%y) ! y-coordinate

            ! left vector: centroid to midpoint vector
            edgeList(edge_id)%C2Mvector(_X_,1) = midPoint(_X_) - elemList(leftTri)%centroids(_X_) ! x-coordinate
            edgeList(edge_id)%C2Mvector(_Y_,1) = midPoint(_Y_) - elemList(leftTri)%centroids(_Y_) ! y-coordinate

            ! right vector: centroid to midpoint vector
            edgeList(edge_id)%C2Mvector(_X_,2) = midPoint(_X_) - elemList(rightTri)%centroids(_X_) ! x-coordinate
            edgeList(edge_id)%C2Mvector(_Y_,2) = midPoint(_Y_) - elemList(rightTri)%centroids(_Y_) ! y-coordinate
        end do

        do i = 1,numBoundary
            edge_id = boundaryEdgeList(i)
            leftTri = edgeList(edge_id)%e1

            ! midpoint of edge
            midPoint(_X_) = half * (nodeList(edgeList(edge_id)%n2)%x + nodeList(edgeList(edge_id)%n1)%x) ! x-coordinate
            midPoint(_Y_) = half * (nodeList(edgeList(edge_id)%n2)%y + nodeList(edgeList(edge_id)%n1)%y) ! y-coordinate

            !Left R vector: centroid to midpoint vector
            edgeList(edge_id)%C2Mvector(_X_,1) = midPoint(_X_) - elemList(leftTri)%centroids(_X_)
            edgeList(edge_id)%C2Mvector(_Y_,1) = midPoint(_Y_) - elemList(leftTri)%centroids(_Y_)
        end do
    end subroutine

    subroutine geometry_edge_normals(numNode,numEdge,nodeList,edgeList,edgeNormals)
        use my_kinddefs
        use grid_data_type_module, only: node_type,edge_type
        implicit none

        integer(i4),       intent(in) :: numNode
        integer(i4),       intent(in) :: numEdge
        type(node_type),   intent(in) :: nodeList(numNode)
        type(edge_type),   intent(in) :: edgeList(numEdge)
        real(dp),pointer,  intent(out):: edgeNormals(:,:)

        integer(i4) :: n1,n2
        integer(i4) :: i

        allocate(edgeNormals(2,numEdge))

        do i = 1,numEdge
            n1 = edgeList(i)%n1 !starting node of edge
            n2 = edgeList(i)%n2 !ending   node of edge

            edgeNormals(1,i) = nodeList(n2)%y - nodeList(n1)%y ! nx =   y2 - y1
            edgeNormals(2,i) = nodeList(n1)%x - nodeList(n2)%x ! ny = -(x2 - x1)
        end do
    end subroutine

end module