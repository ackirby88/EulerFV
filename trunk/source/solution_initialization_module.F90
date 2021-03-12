module solution_initialization_module
     implicit none
     contains

    subroutine solution_initialization(nf,nelem)
        use my_kinddefs
        use solution_module
        use inputs_module, only: rho0,rhou0,rhov0,rhoe0
        use inputs_module, only: restart
        implicit none

        integer(i4),intent(in) :: nf
        integer(i4),intent(in) :: nelem

        ! allocate solution and residual arrays
        call solution_allocation(nf,nelem)

        if (restart) then
            call read_solution(nf,nelem,Q)
        else
            call ic_freestream(nf,nelem,rho0,rhou0,rhov0,rhoe0,Q)
        end if
    end subroutine

    subroutine read_solution(nf,numElem,Q)
        use my_kinddefs
        use inputs_module, only: checkpoint_name
        implicit none

        integer(i4),intent(in) :: nf
        integer(i4),intent(in) :: numElem
        real(dp),   intent(out):: Q(nf,numElem)

        character(80) :: comments
        integer(i4)   :: funit,elem_id,numCells

        funit = 22
        open(unit = funit, file = checkpoint_name, ACTION='READ', status = 'old')
        read(funit,*) comments
        read(funit,*) numCells

        if (numCells .ne. numElem) then
            print*,'The input solution file does not contain same number of elements: stopping'
            stop
        end if

         ! read solution
         do elem_id = 1,numElem
            read(funit,*) Q(1,elem_id), & ! rho
                          Q(2,elem_id), & ! rhou
                          Q(3,elem_id), & ! rhov
                          Q(4,elem_id)    ! rhoE
        end do
        close(funit)
        print*,'Restarting solution complete.'
    end subroutine

    subroutine ic_freestream(nf,npts,rho0,rhou0,rhov0,rhoe0,Q)
        use my_kinddefs
        implicit none

        integer(i4),intent(in)  :: nf,npts
        real(dp),   intent(in)  :: rho0
        real(dp),   intent(in)  :: rhou0
        real(dp),   intent(in)  :: rhov0
        real(dp),   intent(in)  :: rhoe0
        real(dp),   intent(out) :: Q(nf,npts)

        integer(i4) :: pt

        Q = 0.0
        do pt = 1,npts
            Q(1,pt) = rho0
            Q(2,pt) = rhou0
            Q(3,pt) = rhov0
            Q(4,pt) = rhoe0
        end do
    end subroutine

end module