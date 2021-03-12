module time_step_module
    implicit none
    contains

    subroutine time_step(tacc,nf,nelem,Q,elem_h,cfl,dt)
        use my_kinddefs
        use inputs_module,  only: gamma,gm1
        implicit none

        logical,    intent(in) :: tacc
        integer(i4),intent(in) :: nf
        integer(i4),intent(in) :: nelem
        real(dp),   intent(in) :: Q(nf,nelem)
        real(dp),   intent(in) :: elem_h(nelem)
        real(dp),   intent(in) :: cfl
        real(dp),   intent(out):: dt(nelem)

        real(dp)    :: pOverRho,c
        real(dp)    :: u,v,oneOrho
        real(dp)    :: min_dt
        integer(i4) :: triNumber

        ! default max
        min_dt = 1.0E+10

        do triNumber = 1,nelem
            oneOrho = 1.0/Q(1,triNumber)
            u   = Q(2,triNumber)*oneOrho
            v   = Q(3,triNumber)*oneOrho

            pOverRho = gm1 * (Q(4,triNumber)*oneOrho - 0.5*(u*u + v*v))
            c = sqrt(gamma* pOverRho)

            dt(triNumber) = cfl * elem_h(triNumber)/(sqrt(u*u + v*v) + c)
            min_dt = min(min_dt,dt(triNumber))
        end do

        ! time accurate: overwrite first dt
        if(tacc) dt(1) = min_dt
    end subroutine
end module