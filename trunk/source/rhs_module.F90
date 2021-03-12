module rhs_module
    implicit none
    contains

    subroutine rhs(order,Q,R)
        use my_kinddefs
        use bc_module
        use flux_module
        use gradient_module
        use solution_module, only: dQ
        use globals_module, only: numFields
        use globals_module, only: numElem,numEdge,numInterior,numBoundary
        use globals_module, only: edgeList,interiorEdgeList,boundaryEdgeList
        use globals_module, only: edgeNormals
        use inputs_module,  only: fluxScheme,outer_bc,wall_bc
        use limiter_module, only: limiter
        implicit none

        integer(i4),intent(in) :: order                 ! spatial order of accuracy
        real(dp),   intent(in) :: Q(numFields,numElem)  ! solution vector
        real(dp),   intent(out):: R(numFields,numElem)  ! residual vector

        real(dp)    :: qB(numFields)
        real(dp)    :: flux(numFields)
        real(dp)    :: fluxLow(numFields)
        real(dp)    :: fluxHigh(numFields)
        real(dp)    :: qLMid(numFields),qRMid(numFields)
        real(dp)    :: phi
        integer(i4) :: edge_id,leftTri,rightTri,edgeNum,bc_type

        ! initialize
        R = 0.0

        ! ==================== !
        ! Interior Edge Fluxes !
        ! ==================== !
        if (order == 2) then
            ! ===================== !
            ! Second Order Residual !
            ! ===================== !
            ! flow gradients
            call gradient(Q,dQ)

            do edge_id = 1,numInterior
                edgeNum = interiorEdgeList(edge_id)
                leftTri  = edgeList(edgeNum)%e1
                rightTri = edgeList(edgeNum)%e2

                ! limit gradients
                call limiter(dQ(:,:,leftTri),dQ(:,:,rightTri),phi)

                ! construct the midpoint solutions using the gradient extensions
                qLMid(:) = Q(:,leftTri)  + 0.8*matmul(dQ(:,:,leftTri) ,edgeList(edgeNum)%C2Mvector(:,1))
                qRMid(:) = Q(:,rightTri) + 0.8*matmul(dQ(:,:,rightTri),edgeList(edgeNum)%C2Mvector(:,2))

                call flux_riemann(fluxScheme,numFields,edgeNormals(:,edgeNum),qLMid,qRMid,fluxHigh)
                call flux_riemann(fluxScheme,numFields,edgeNormals(:,edgeNum),Q(:,leftTri),Q(:,rightTri),fluxLow)
                flux = fluxLow - phi*(fluxLow - fluxHigh)

                R(:,leftTri)  = R(:,leftTri)  + flux(:)
                R(:,rightTri) = R(:,rightTri) - flux(:)
            end do
        else
            ! ==================== !
            ! First Order Residual !
            ! ==================== !
            do edge_id = 1,numInterior
                edgeNum = interiorEdgeList(edge_id)
                leftTri  = edgeList(edgeNum)%e1
                rightTri = edgeList(edgeNum)%e2

                call flux_riemann(fluxScheme,numFields,edgeNormals(:,edgeNum),Q(:,leftTri),Q(:,rightTri),flux)

                R(:,leftTri)  = R(:,leftTri)  + flux(:)
                R(:,rightTri) = R(:,rightTri) - flux(:)
            end do
        end if

        ! ==================== !
        ! Boundary Edge Fluxes !
        ! ==================== !
        do edge_id = 1,numBoundary
            edgeNum = boundaryEdgeList(edge_id)
            bc_type = edgeList(edgeNum)%edgeType
            leftTri = edgeList(edgeNum)%e1

            select case(bc_type)
                case(9) !outer boundary
                         ! [4] characteristic
                    call bc(outer_bc, edgeNormals(:,edgeNum),Q(:,leftTri),qB)
                    call flux_bc(numFields,qB,edgeNormals(:,edgeNum),flux)

                    R(:,leftTri) = R(:,leftTri) + flux(:)
                case(1) !wall boundary
                         ! [2] slip wall
                         ! [3] nonslip wall
                    call bc(wall_bc,edgeNormals(:,edgeNum),Q(:,leftTri),qB)
                    call flux_bc(numFields,qB,edgeNormals(:,edgeNum),flux)

                    R(:,leftTri) = R(:,leftTri) + flux(:)
            end select
        end do
    end subroutine
end module