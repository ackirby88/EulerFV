module output_module
    implicit none
    contains

    subroutine output_solution(project,solCoeffs)
        use my_kinddefs
        use globals_module
        use gradient_module
        use solution_module, only: dQ
        implicit none

        character(20),intent(in) :: project
        real(dp),     intent(in) :: solCoeffs(numFields,numElem)

        !---> Local Variables
        character(80) :: filename
        character(80) :: file_num
        integer(i4)   :: iunit
        integer(i4)   :: io
        integer(i4)   :: n,n1,n2,n3

        real(dp) :: ddensity,u,v,vort,oneOrho,rho,pressure,gamma,kappa,e,mach
        real(dp) :: dudx,dvdx,dudy,dvdy
        real(dp) :: vv1,vv2,vv3

        gamma = 1.4_dp
        kappa = gamma - 1.0_dp

        print*,'# Writing to File...',printCount
        ! ================== !
        ! solution gradients !
        ! ================== !
        call gradient(solCoeffs,dQ)

        ! ================ !
        ! file information !
        ! ================ !
        write(file_num,'(I4.4)') printCount
        printCount = printCount + 1
        filename = 'WRK/' // adjustr(trim(project)) // '_' // &
                   adjustr(trim(file_num)) // '.vtu'

        iunit = 20
        open(unit = iunit, file = filename, status = 'replace', iostat = io)
        if (io /= 0) then
           print*,'ERROR: Opening soln vtu file'
        else

            write(iunit,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
            write(iunit,*)' <UnstructuredGrid>'
            write(iunit,*)'  <Piece NumberOfPoints="',3*numElem,'" NumberOfCells="',numElem,'">'
            write(iunit,*)'   <CellData>'
            write(iunit,*)'    <DataArray type="Float32" Name="Mach" Format="ascii">'
                do n = 1,numElem
                    rho      = solCoeffs(1,n);
                    oneOrho  = 1.0_dp / rho
                    u        = solCoeffs(2,n)*oneOrho
                    v        = solCoeffs(3,n)*oneOrho;
                    e        = solCoeffs(4,n)*oneOrho;
                    pressure = kappa*(solCoeffs(4,n) - 0.5_dp*rho*(u*u + v*v))
                    mach     = sqrt(u*u + v*v)/sqrt(gamma*pressure*oneOrho)
                    write(iunit,*)'    ',mach
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Float32" Name="Pressure" Format="ascii">'
                do n = 1,numElem
                    rho      = solCoeffs(1,n);
                    oneOrho  = 1.0_dp / rho
                    u        = solCoeffs(2,n)*oneOrho
                    v        = solCoeffs(3,n)*oneOrho;
                    e        = solCoeffs(4,n)*oneOrho;
                    pressure = kappa*(solCoeffs(4,n) - 0.5_dp*rho*(u*u + v*v))
                    write(iunit,*)'    ',pressure
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Float32" Name="Density" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',solCoeffs(1,n)
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Float32" Name="Vel-U" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',solCoeffs(2,n)/solCoeffs(1,n)
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Float32" Name="Vel-V" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',solCoeffs(3,n)/solCoeffs(1,n)
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Float32" Name="Energy Density" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',solCoeffs(4,n)
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Float32" Name="Density Gradient" Format="ascii">'
                do n = 1,numElem
                    ddensity = sqrt(dQ(1,1,n)*dQ(1,1,n) + dQ(1,2,n)*dQ(1,2,n))
                    write(iunit,*)'    ',ddensity
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Float32" Name="Vorticity" Format="ascii">'
                do n = 1, numElem
                    oneOrho = 1.0 / solCoeffs(1,n)
                    u = solCoeffs(2,n)*oneOrho
                    v = solCoeffs(3,n)*oneOrho

                    dvdx = (dQ(3,1,n) - v*dQ(1,1,n))*oneOrho
                    dudy = (dQ(2,2,n) - u*dQ(1,2,n))*oneOrho

                    vv1 = dvdx
                    vv2 = dudy
                    vv3 = dvdx-dudy
                    vort = sqrt(vv1*vv1 + vv2*vv2 + vv3*vv3)
                    write(iunit,*)'    ',vort
                end do
            write(iunit,*)'    </DataArray>'
    !       write(iunit,*)'    <DataArray type="Int32" Name="BC" Format="ascii">'
    !           do elem_id = 1, numElem
    !               isBoundElem = .false.
    !               do edgeCount = 1, extCount-1
    !                    bdedge =  boundaryEdgeList(edgeCount)
    !                    boundaryTri = edgeList(bdedge)%e1
    !                    if (elem_id == boundaryTri) then
    !                        isBoundElem = .true.
    !                        bc_type = bcFlag(bdedge)
    !                    end if
    !               end do
    !
    !               select case(isBoundElem)
    !                   case(.true.)
    !                        write(iunit,*)'    ',bc_type
    !                   case(.false.)
    !                        write(iunit,*)'    ',0
    !               end select
    !           end do
    !       write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Int32" Name="Element Number" Format="ascii">'
                do n = 1, numElem
                    write(iunit,*)'    ',n
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'   </CellData>'
            write(iunit,*)'   <Points>'
            write(iunit,*)'    <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
                do n = 1,numElem
                    n1 = elemList(n)%vtxList(1)
                    n2 = elemList(n)%vtxList(2)
                    n3 = elemList(n)%vtxList(3)

                    write(iunit,*)'    ', nodeList(n1)%x, nodeList(n1)%y, 0.0 &
                                        , nodeList(n2)%x, nodeList(n2)%y, 0.0 &
                                        , nodeList(n3)%x, nodeList(n3)%y, 0.0
                end do
            write(iunit,*) '    </DataArray>'
            write(iunit,*) '   </Points>'
            write(iunit,*) '   <Cells>'
            write(iunit,*) '    <DataArray type="Int32" Name="connectivity" Format="ascii">'
                do n = 1,numElem
                   write(iunit,*)'    ',3*n-3,3*n-2,3*n-1
                end do
            write(iunit,*) '    </DataArray>'
            write(iunit,*) '    <DataArray type="Int32" Name="offsets" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',3*n
                end do
            write(iunit,*) '    </DataArray>'
            write(iunit,*) '    <DataArray type="Int32" Name="types" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',5
                end do
            write(iunit,*) '    </DataArray>'
            write(iunit,*) '   </Cells>'
            write(iunit,*) '  </Piece>'
            write(iunit,*) ' </UnstructuredGrid>'
            write(iunit,*) '</VTKFile>'
        end if
        close(iunit)
    end subroutine

    subroutine output_residual(project,res)
        use my_kinddefs
        use globals_module
        implicit none

        character(*),intent(in) :: project
        real(dp),    intent(in) :: res(numFields,numElem)

        !---> Local Variables
        character(80) :: filename
        character(80) :: file_num
        integer(i4)   :: iunit
        integer(i4)   :: n,io

        real(dp)    :: eps = 1.0e-13
        integer(i4) :: elem_id,boundaryTri,edgeCount,bdedge,bc_type,n1,n2,n3
        logical     :: isBoundElem

        write(file_num,'(I4.4)') printCountRes

        printCountRes = printCountRes + 1
        filename = 'WRK/' // adjustr(trim(project)) // '_Res_' // &
                             adjustr(trim(file_num)) // '.vtu'

        iunit = 22
        open(unit = iunit, file = filename, status = 'replace', iostat = io)
        if (io /= 0) then
           print*,'ERROR: Opening soln vtu file'
        else

            write(iunit,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
            write(iunit,*)' <UnstructuredGrid>'
            write(iunit,*)'  <Piece NumberOfPoints="',3*numElem,'" NumberOfCells="',numElem,'">'
            write(iunit,*)'   <CellData>'
            write(iunit,*)'    <DataArray type="Float32" Name="Res_Rho" Format="ascii">'
                do n = 1,numElem
                write(iunit,*)'    ',abs(res(1,n))+eps
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Float32" Name="Res_U" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',abs(res(2,n))+eps
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Float32" Name="Res_V" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',abs(res(3,n))+eps
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Float32" Name="Res_rhoE" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',abs(res(4,n))+eps
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'    <DataArray type="Int32" Name="Element Number" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',n
                end do
            write(iunit,*)'    </DataArray>'
            write(iunit,*)'   </CellData>'
            write(iunit,*)'   <Points>'
            write(iunit,*)'    <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
                do n = 1,numElem
                    n1 = elemList(n)%vtxList(1)
                    n2 = elemList(n)%vtxList(2)
                    n3 = elemList(n)%vtxList(3)
                    write(iunit,*)'    ', nodeList(n1)%x, nodeList(n1)%y, 0.0 &
                                        , nodeList(n2)%x, nodeList(n2)%y, 0.0 &
                                        , nodeList(n3)%x, nodeList(n3)%y, 0.0
                end do
            write(iunit,*) '    </DataArray>'
            write(iunit,*) '   </Points>'
            write(iunit,*) '   <Cells>'
            write(iunit,*) '    <DataArray type="Int32" Name="connectivity" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',3*n-3,3*n-2,3*n-1
                end do
            write(iunit,*) '    </DataArray>'
            write(iunit,*) '    <DataArray type="Int32" Name="offsets" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',3*n
                end do
            write(iunit,*) '    </DataArray>'
            write(iunit,*) '    <DataArray type="Int32" Name="types" Format="ascii">'
                do n = 1,numElem
                    write(iunit,*)'    ',5
                end do
            write(iunit,*) '    </DataArray>'
            write(iunit,*) '   </Cells>'
            write(iunit,*) '  </Piece>'
            write(iunit,*) ' </UnstructuredGrid>'
            write(iunit,*) '</VTKFile>'
        end if
        close(iunit)
    end subroutine

    subroutine writeResNorm(project,norm)
        use my_kinddefs
        implicit none

        character(*),intent(in) :: project
        real(dp),    intent(in) :: norm

        !---> Local Variables
        character(80) :: filename
        character(80) :: file_num
        integer(i4)   :: io,tm,elem_id
        integer(i4)   :: iunit,i

        iunit = 22
        filename = 'WRK/' // adjustr(trim(project)) // '_l2hist.txt'
        open(unit = iunit, file = filename, status = 'UNKNOWN', ACCESS='APPEND', iostat = io)
            write(iunit,*) norm
        close(iunit)
    end subroutine

    subroutine checkpoint(project,solCoeffs)
        use my_kinddefs
        use globals_module, only: printCountCoeffs,numElem
        implicit none

        character(*),intent(in) :: project
        real(dp),    intent(in) :: solCoeffs(:,:) !res(numEulerVars,numElem)

        !---> Local Variables
        character(80) :: filename
        character(80) :: file_num
        integer(i4)   :: io,tm,elem_id
        integer(i4)   :: iunit

        write(file_num,'(I4.4)') printCountCoeffs

        iunit = 22
        printCountCoeffs = printCountCoeffs + 1
        filename = 'WRK/' // adjustr(trim(project)) // '_' // adjustr(trim(file_num)) // '.checkpoint'
        open(unit = iunit, file = filename, status = 'replace', iostat = io)

        if (io /= 0) then
           print*,'ERROR: Opening soln vtu file'
        else
            write(iunit,*)' Solution Coefficients for file: ',filename
            write(iunit,*) numElem
            do elem_id = 1, numElem
                    write(iunit,*) solCoeffs(1,elem_id),solCoeffs(2,elem_id),&
                                   solCoeffs(3,elem_id),solCoeffs(4,elem_id)
            end do
        end if
        close(iunit)
    end subroutine

    subroutine nacaSurfacePressure(project,solCoeffs)
        use my_kinddefs
        use inputs_module,  only: gamma,pressure0,density0,u0,v0,gm1
        use globals_module, only: edgeList,numBoundary,numEdge,boundaryEdgeList,nodeList
        implicit none

        character(*),intent(in) :: project
        real(dp),    intent(in) :: solCoeffs(:,:) !res(numEulerVars,numElem)

        integer(i4) :: i,edgeNum,bc_type,leftTri
        real(dp)  :: p,x_ave
        real(dp)  :: rho,u,v

        !---> Local Variables
        character(80) :: filename
        character(80) :: file_num
        integer(i4)   :: io,tm,elem_id
        integer(i4)   :: iunit

        iunit = 22
        filename = 'WRK/' // adjustr(trim(project)) // '_Cp.dat'
        open(unit = iunit, file = filename, status = 'replace', iostat = io)

        if (io /= 0) then
           print*,'ERROR: Opening Cp File'
        else
            write(iunit,*) 'TITLE = "Surface Pressure Coefficient: (x,C_p) ordered pairs" '
            write(iunit,*) 'Variables = x, Cp'

            do i = 1,numBoundary
                edgeNum = boundaryEdgeList(i)
                bc_type = edgeList(edgeNum)%edgeType
                leftTri = edgeList(edgeNum)%e1
                x_ave = half*(nodeList(edgeList(edgeNum)%n1)%x + nodeList(edgeList(edgeNum)%n2)%x)

                 if (bc_type == 1) then !wall boundary
                     !----> Calculate pressure: print x_ave,pressure
                        rho = solCoeffs(1,leftTri)
                        u   = solCoeffs(2,leftTri)/rho
                        v   = solCoeffs(3,leftTri)/rho
                        p   = gm1*(solCoeffs(4,leftTri) - 0.5_dp*rho*(u*u + v*v))

                        write(iunit,*) x_ave, -((p-pressure0)/(half*density0*(u0*u0 + v0*v0)))
                 end if
           end do
        end if
        close(iunit)
  end subroutine
end module
