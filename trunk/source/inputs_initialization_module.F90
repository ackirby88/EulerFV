module inputs_initialization_module
    implicit none
    contains

    subroutine inputs_initialization()
        use my_kinddefs
        use inputs_module
        implicit none

        character(80) :: file_command

        ! initialize flow variables
        call initialize_inputs()
        call initialize_flow_info()

        ! make WRK directory
        file_command = 'mkdir WRK'
        call system(file_command)
    end subroutine

    subroutine initialize_flow_info()
        use inputs_module
        implicit none

        gm1 = gamma - 1.0
        gOgm1 = gamma / gm1
        oneOgm1 = 1.0 / gm1

        C0 = sqrt(gamma * pressure0 / density0)
        U0 = C0*fmach*cos(alpha*pi/180.0)
        V0 = C0*fmach*sin(alpha*pi/180.0)

        rho0  = density0
        rhou0 = density0*u0
        rhov0 = density0*v0
        rhoE0 = (pressure0 * oneOgm1) + 0.5*density0*(u0*u0 + v0*v0)
    end subroutine

    subroutine initialize_inputs()
        use my_kinddefs
        use inputs_module,  only: input_name
        implicit none

        integer(i4)   :: argc
        character(30) :: fmt

        fmt = "(A,A)"
        argc = command_argument_count()
        if (argc == 1) then
            call get_command_argument(1,input_name)
            write(*,fmt) "Input file read: ",trim(input_name)
        else
            print*," Please provide input file: STOPPING"
            stop
        end if

        call read_inputs(input_name)
    end subroutine

    subroutine read_inputs(file_name)
        use my_kinddefs
        use inputs_module
        implicit none

        character(80),intent(in) :: file_name
        integer(i4) :: funit,i,ind

        namelist/INPUTS/    stats_name,nodes_name,cells_name,edges_name,         &
                            gamma,fmach,alpha,density0,pressure0,                &
                            restart,checkpoint_name,                             &
                            limiterScheme,fluxScheme,spatial_order,              &
                            outer_bc,wall_bc,                                    &
                            output_file,output_freq,checkpoint_freq,logging_freq,&
                            CFL,timeSteps,timeScheme,tacc

        funit = 22
        open(unit=funit,file=file_name,status='old',form='formatted')
        read(funit,INPUTS)
        close(funit)

        write(*,*)
        write(*,*) "=================================================== "
        write(*,*) "! ----------------- Mesh Fields ----------------- !"
        write(*,*) "stats_name: ",stats_name
        write(*,*) "nodes_name: ",nodes_name
        write(*,*) "cells_name: ",cells_name
        write(*,*) "edges_name: ",edges_name
        write(*,*) "restart: ",restart
        write(*,*) "checkpoint_name: ",checkpoint_name
        write(*,*) "! --------------- Flow Variables ---------------- !"
        write(*,*) "gamma: ",gamma
        write(*,*) "fmach: ",fmach
        write(*,*) "alpha: ",alpha
        write(*,*) "density0: ",density0
        write(*,*) "pressure0: ",pressure0
        write(*,*) "limiterScheme: ",limiterScheme
        write(*,*) "fluxScheme: ",fluxScheme
        write(*,*) "spatial_order: ",spatial_order
        write(*,*) "outer_bc: ",outer_bc
        write(*,*) "wall_bc: ",wall_bc
        write(*,*)
        write(*,*) "! ---------- Plot/Checkpoint Options ------------ !"
        write(*,*) "logging_freq: ",logging_freq
        write(*,*) "output_file: ",output_file
        write(*,*) "output_freq: ",output_freq
        write(*,*) "checkpoint_freq: ",checkpoint_freq
        write(*,*)
        write(*,*) "! ----------- Time Stepping Options ------------- !"
        write(*,*) "CFL: ",CFL
        write(*,*) "timeSteps: ",timeSteps
        write(*,*) "timeScheme: ",timeScheme
        write(*,*) "tacc: ",tacc
        write(*,*) "=================================================== "
        write(*,*)
    end subroutine

end module