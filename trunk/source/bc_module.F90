module bc_module
    implicit none
    contains

    subroutine bc(bc_type,nxy,q,qb)
        use my_kinddefs
        use inputs_module
        use globals_module
        implicit none

        integer(i4),intent(in) :: bc_type
        real(dp),intent(in)    :: nxy(2)
        real(dp),intent(in)    :: q(4)
        real(dp),intent(out)   :: qb(4)

        !---> Local Variables
        real(dp) :: oneOmag,nx,ny
        real(dp) :: rho,oneOrho,u,v,p,c
        real(dp) :: oneOrho0

        real(dp) :: rhob,ub,vb
        real(dp) :: un,uref
        real(dp) :: riem_p,riem_m,unb,cb,sb,pb

        qb = 0.0
        select case(bc_type)
            case(1) !free-stream
                ! ============== !
                ! interior state !
                ! ============== !
                rho = q(1)
                oneOrho = 1.0/rho
                u = q(2)*oneOrho
                v = q(3)*oneOrho

                ! ============== !
                ! boundary state !
                ! ============== !
                qb(1) = rho
                qb(2) = rho*u0
                qb(3) = rho*v0
                qb(4) = q(4) - 0.5*rho*(u*u + v*v)
                return
            case(2) !slip wall
                ! =========== !
                ! face normal !
                ! =========== !
                oneOmag = 1.0 / sqrt((nxy(1)*nxy(1) + nxy(2)*nxy(2)))
                nx = nxy(1)*oneOmag
                ny = nxy(2)*oneOmag

                ! ============== !
                ! interior state !
                ! ============== !
                rho = q(1)
                oneOrho = 1.0/rho
                u = q(2)*oneOrho
                v = q(3)*oneOrho

                ! ================= !
                ! boundary velocity !
                ! ================= !
                ub = u - (u*nx + v*ny)*nx
                vb = v - (u*nx + v*ny)*ny

                ! ============== !
                ! boundary state !
                ! ============== !
                qb(1) = rho
                qb(2) = rho*ub
                qb(3) = rho*vb
                qb(4) = q(4)
                return
            case(3) !no-slip wall
                ! ============== !
                ! interior state !
                ! ============== !
                rho = q(1)
                oneOrho = 1.0/rho
                u = q(2)*oneOrho
                v = q(3)*oneOrho

                ! ============== !
                ! boundary state !
                ! ============== !
                qb(1) = rho
                qb(2) = 0.0
                qb(3) = 0.0
                qb(4) = q(4) - 0.5*rho*(u*u + v*v)
                return
            case(4) !characteristic
                ! =========== !
                ! face normal !
                ! =========== !
                oneOmag = 1.0 / sqrt((nxy(1)*nxy(1) + nxy(2)*nxy(2)))
                nx = nxy(1)*oneOmag
                ny = nxy(2)*oneOmag

                ! ============== !
                ! interior state !
                ! ============== !
                rho = q(1)
                oneOrho = 1.0/rho
                u = q(2)*oneOrho
                v = q(3)*oneOrho
                p = gm1*(q(4) - 0.5*rho*(u*u + v*v))
                c = sqrt(gamma*p*oneOrho)

                ! ============= !
                ! initial state !
                ! ============= !
                oneOrho0 = 1.0/rho0

                ! =============== !
                ! normal velocity !
                ! =============== !
                un = u*nx + v*ny

                ! ================== !
                ! reference velocity !
                ! ================== !
                uref = u0*nx + v0*ny

                ! ================== !
                ! Riemann invariants !
                ! ================== !
                riem_p = un   + 2.0*oneOgm1*c
                riem_m = uref - 2.0*oneOgm1*C0
                unb    = 0.5     *(riem_p + riem_m)
                cb     = 0.25*gm1*(riem_p - riem_m)

                ! ==================== !
                ! check flow condition !
                ! ==================== !
                if (abs(unb) < abs(cb)) then
                    !subsonic
                    if (unb < 0.0_dp) then
                        !subsonic inflow
                        sb = pressure0*oneOrho0**gamma
                        ub = u0 + (unb-uref)*nx
                        vb = v0 + (unb-uref)*ny
                    else
                        !subsonic outflow
                        sb = p*oneOrho**gamma
                        ub = u + (unb-un)*nx
                        vb = v + (unb-un)*ny
                    end if

                    rhob = (cb*cb/sb/gamma)**(oneOgm1)
                    pb   = sb*rhob**gamma

                    qb(1) = rhob
                    qb(2) = rhob*ub
                    qb(3) = rhob*vb
                    qb(4) = pb*oneOgm1 + 0.5*rhob*(ub*ub + vb*vb)
                else
                    !supersonic
                    if (unb < 0.0) then
                        ! supersonic inflow: imposed from outside
                        qb(1) = rho0
                        qb(2) = rhou0
                        qb(3) = rhov0
                        qb(4) = rhoe0
                    else
                        ! supersonic outflow: imposed from inside
                        qb(1:4) = q(1:4)
                    end if
                end if
                return
        end select
    end subroutine
end module