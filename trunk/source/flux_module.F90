module flux_module
    implicit none
    contains

    subroutine flux_riemann(scheme,nfields,normals,qL,qR,flux)
        use my_kinddefs
        implicit none

        integer(i4),intent(in) :: scheme
        integer(i4),intent(in) :: nfields
        real(dp),   intent(in) :: normals(2)
        real(dp),   intent(in) :: qL(nfields)
        real(dp),   intent(in) :: qR(nfields)
        real(dp),   intent(out):: flux(nfields)

        select case(scheme)
            case(1)
                call flux_roe(nfields,normals,qL,qR,flux)
            case(2)
                call flux_lax_friedrichs(nfields,normals,qL,qR,flux)
            case default
                print*,'Choose either [1] Roe, [2] LaxF for flux'
                stop
        end select
    end subroutine

    subroutine flux_bc(nfields,q,normals,flux)
        use my_kinddefs
        use inputs_module, only: gm1
        implicit none

        integer(i4),intent(in) :: nfields
        real(dp),   intent(in) :: q(nfields)
        real(dp),   intent(in) :: normals(2)
        real(dp),   intent(out):: flux(nfields)

        real(dp) :: fluxX(4),fluxY(4)
        real(dp) :: oneOrho,u,v,p
        real(dp) :: nx,ny

        ! normal vector (keep dimensional)
        nx = normals(1)
        ny = normals(2)

        oneOrho = 1.0/q(1)
        u = q(2)*oneOrho
        v = q(3)*oneOrho
        p = gm1*(q(4) - 0.5*q(1)*(u*u + v*v))

        !convective flux x component
        fluxX(1) = q(1)*u
        fluxX(2) = q(2)*u + p
        fluxX(3) = q(3)*u
        fluxX(4) = q(4)*u + p*u

        !convective flux y component
        fluxY(1) = q(1)*v
        fluxY(2) = q(2)*v
        fluxY(3) = q(3)*v + p
        fluxY(4) = q(4)*v + p*v

        flux = fluxX*nx + fluxY*ny
    end subroutine

    subroutine flux_roe(nf,normals,ql,qr,flux)
        use my_kinddefs
        use inputs_module, only: gm1
        implicit none

        integer(i4),intent(in) :: nf
        real(dp),   intent(in) :: normals(2)
        real(dp),   intent(in) :: ql(nf)
        real(dp),   intent(in) :: qr(nf)
        real(dp),   intent(out):: flux(nf)

        real(dp) :: oneOmag,mag,nx,ny
        real(dp) :: rhol,oneOrhol,ul,vl,pl,hl,unorml,fl(4)
        real(dp) :: rhor,oneOrhor,ur,vr,pr,hr,unormr,fr(4)
        real(dp) :: ubar,vbar,hbar,qqbar,cbar,unormbar
        real(dp) :: fact,A,B,G1,G2,S1,S2,C1,C2
        real(dp) :: eigu,eigumc,eigupc
        real(dp) :: qjump(4),diss(4)
        real(dp) :: eps0 = 1E-2
        real(dp) :: eps

        ! =========== !
        ! face normal !
        ! =========== !
        mag = sqrt(normals(1)*normals(1) + normals(2)*normals(2))
        oneOmag = 1.0/mag
        nx = normals(1)*oneOmag
        ny = normals(2)*oneOmag

        ! ========== !
        ! left state !
        ! ========== !
        rhol = ql(1)
        oneOrhol = 1.0/rhol
        ul  = ql(2)*oneOrhol
        vl  = ql(3)*oneOrhol
        pl  = gm1*(ql(4) - 0.5*rhol*(ul*ul + vl*vl))
        hl  = (ql(4) + pl)*oneOrhol

        ! =========== !
        ! right state !
        ! =========== !
        rhor = qr(1)
        oneOrhor = 1.0/rhor
        ur  = qr(2)*oneOrhor
        vr  = qr(3)*oneOrhor
        pr  = gm1*(qr(4) - 0.5*rhor*(ur*ur + vr*vr))
        hr  = (qr(4) + pr)*oneOrhor

        ! ====================== !
        ! face normal velocities !
        ! ====================== !
        unorml = ul*nx + vl*ny
        unormr = ur*nx + vr*ny

        ! =============== !
        ! jump conditions !
        ! =============== !
        qjump = qr - ql

        ! ================= !
        ! Roe average state !
        ! ================= !
        fact = sqrt(rhor/rhol)
        A    = 1.0/(1.0 + fact)
        B    = fact/(1.0 + fact)

        ubar     = ul*A + ur*B
        vbar     = vl*A + vr*B
        hbar     = hl*A + hr*B
        qqbar    = 0.5*(ubar*ubar + vbar*vbar)
        cbar     = sqrt(gm1*(hbar - qqbar))
        unormbar = ubar*nx + vbar*ny

        ! =========== !
        ! eigenvalues !
        ! =========== !
        eigu   = abs(unormbar)
        eigupc = abs(unormbar + cbar)
        eigumc = abs(unormbar - cbar)

        ! ================== !
        ! Harten entropy fix !
        ! ================== !
        eps = cbar*eps0
        if(eigu   < eps) eigu   = 0.5*(eps + eigu*eigu/eps)
        if(eigupc < eps) eigupc = 0.5*(eps + eigupc*eigupc/eps)
        if(eigumc < eps) eigumc = 0.5*(eps + eigumc*eigumc/eps)

        ! ======== !
        ! Roe flux !
        ! ======== !
        G1 = gm1/(cbar*cbar) * (qjump(1)*qqbar    - qjump(2)*ubar - qjump(3)*vbar + qjump(4))
        G2 =                  (-qjump(1)*unormbar + qjump(2)*nx   + qjump(3)*ny)
        S1 = 0.5*(eigupc + eigumc)
        S2 = 0.5*(eigupc - eigumc)
        C1 = G2*S2/cbar + G1*(S1 - eigu)
        C2 = G1*S2*cbar + G2*(S1 - eigu)

        diss(1) = (                   C1 + qjump(1)*eigu)
        diss(2) = (      nx*C2 + ubar*C1 + qjump(2)*eigu)
        diss(3) = (      ny*C2 + vbar*C1 + qjump(3)*eigu)
        diss(4) = (unormbar*C2 + hbar*C1 + qjump(4)*eigu)

        ! ============= !
        ! native fluxes !
        ! ============= !
        fl(1) = rhol*unorml
        fl(2) = rhol*unorml*ul + nx*pl
        fl(3) = rhol*unorml*vl + ny*pl
        fl(4) = (ql(4) + pl) * unorml

        fr(1) = rhor*unormr
        fr(2) = rhor*unormr*ur + nx*pr
        fr(3) = rhor*unormr*vr + ny*pr
        fr(4) = (qr(4) + pr) * unormr

        ! ============ !
        ! Riemann Flux !
        ! ============ !
        flux = 0.5*(fl + fr) - 0.5*diss

        ! rescale flux to incorporate edge length for surface integration
        flux = mag*flux
    end subroutine

    subroutine flux_lax_friedrichs(nf,normals,ql,qr,flux)
        use my_kinddefs
        use inputs_module, only: gamma,gm1
        implicit none

        integer(i4),intent(in) :: nf
        real(dp),   intent(in) :: normals(2)
        real(dp),   intent(in) :: ql(nf)
        real(dp),   intent(in) :: qr(nf)
        real(dp),   intent(out):: flux(nf)

        real(dp) :: oneOmag,mag,nx,ny
        real(dp) :: sm
        real(dp) :: rhol,oneOrhol,ul,vl,pl,cl,unorml,fl(4)
        real(dp) :: rhor,oneOrhor,ur,vr,pr,cr,unormr,fr(4)

        ! =========== !
        ! face normal !
        ! =========== !
        mag = sqrt(normals(1)*normals(1) + normals(2)*normals(2))
        oneOmag = 1.0/mag
        nx = normals(1)*oneOmag
        ny = normals(2)*oneOmag

        ! ========== !
        ! left state !
        ! ========== !
        rhol = ql(1)
        oneOrhol = 1.0/rhol
        ul  = ql(2)*oneOrhol
        vl  = ql(3)*oneOrhol
        pl  = gm1*(ql(4) - 0.5*rhol*(ul*ul + vl*vl))
        cl  = sqrt(gamma*pl*oneOrhol)

        ! =========== !
        ! right state !
        ! =========== !
        rhor = qr(1)
        oneOrhor = 1.0/rhor
        ur  = qr(2)*oneOrhor
        vr  = qr(3)*oneOrhor
        pr  = gm1*(qr(4) - 0.5*rhor*(ur*ur + vr*vr))
        cr  = sqrt(gamma*pr*oneOrhor)

        ! ====================== !
        ! face normal velocities !
        ! ====================== !
        unorml = ul*nx + vl*ny
        unormr = ur*nx + vr*ny

        ! ============= !
        ! native fluxes !
        ! ============= !
        fl(1) = rhol*unorml
        fl(2) = rhol*unorml*ul + nx*pl
        fl(3) = rhol*unorml*vl + ny*pl
        fl(4) = (ql(4) + pl) * unorml

        fr(1) = rhor*unormr
        fr(2) = rhor*unormr*ur + nx*pr
        fr(3) = rhor*unormr*vr + ny*pr
        fr(4) = (qr(4) + pr) * unormr

        ! ================ !
        ! dissipative flux !
        ! ================ !
        ! max eigenvalue
        sm = max(abs(unorml) + cl,abs(unormr) + cr)

        ! ============ !
        ! Riemann Flux !
        ! ============ !
        flux = 0.5*(fl + fr) - 0.5*sm*(qr - ql)

        ! rescale flux to incorporate edge length for surface integration
        flux = flux*mag
    end subroutine
end module