module gradient_module
    implicit none
    contains

    subroutine gradient(Q,grad)
        use my_kinddefs
        use globals_module, only: numFields
        use globals_module, only: numElem,elemList
        use globals_module, only: numInterior,numBoundary
        use globals_module, only: edgeList,interiorEdgeList,boundaryEdgeList
        use globals_module, only: edgeNormals
        use inputs_module,  only: outer_bc,wall_bc
        use bc_module
        implicit none

        real(dp),intent(in) :: Q(numFields,numElem)
        real(dp),intent(out):: grad(numFields,2,numElem)

        real(dp) :: QBarDotN(numFields),QB(numFields)
        integer(i4) :: i,edgeNum,leftTri,rightTri,bc_type

        ! initialize
        grad = 0.0

        do i = 1,numInterior
            edgeNum = interiorEdgeList(i)
            leftTri  = edgeList(edgeNum)%e1
            rightTri = edgeList(edgeNum)%e2

            QBarDotN = half * (Q(:,leftTri) + Q(:,rightTri))*edgeNormals(1,edgeNum)
            grad(:,1,leftTri)  = grad(:,1,leftTri)  + QBarDotN
            grad(:,1,rightTri) = grad(:,1,rightTri) - QBarDotN

            QBarDotN = half * (Q(:,leftTri) + Q(:,rightTri))*edgeNormals(2,edgeNum)
            grad(:,2,leftTri)  = grad(:,2,leftTri)  + QBarDotN
            grad(:,2,rightTri) = grad(:,2,rightTri) - QBarDotN
        end do

        do i = 1,numBoundary
            edgeNum = boundaryEdgeList(i)
            leftTri = edgeList(edgeNum)%e1
            bc_type = edgeList(edgeNum)%edgeType

            !form QB
            select case(bc_type)
                case(9) ! outer boundary
                    call bc(outer_bc, edgeNormals(:,edgeNum), Q(:,leftTri), QB)
                case(1) ! wall boundary
                    call bc(wall_bc, edgeNormals(:,edgeNum), Q(:,leftTri), QB)
            end select

            QBarDotN = half * (Q(:,leftTri) + QB)*edgeNormals(1,edgeNum)
            grad(:,1,leftTri)  = grad(:,1,leftTri)  + QBarDotN

            QBarDotN = half * (Q(:,leftTri) + QB)*edgeNormals(2,edgeNum)
            grad(:,2,leftTri)  = grad(:,2,leftTri)  + QBarDotN
        end do

        do i = 1, numElem
            grad(:,:,i) = grad(:,:,i) / elemList(i)%area
        end do
    end subroutine
end module