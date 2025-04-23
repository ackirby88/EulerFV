module evolve_solution_module
    implicit none
    contains

    subroutine evolve_solution(time_scheme,num_steps)
        use my_kinddefs
        use solution_module, only: Q,R
        implicit none

        integer(i4),intent(in) :: time_scheme
        integer(i4),intent(in) :: num_steps

        select case(time_scheme)
            case(1)
                call RK1(num_steps,Q,R)
        end select

    end subroutine

    subroutine RK1(timeSteps,Q,R)
        use my_kinddefs
        use rhs_module
        use norm_module
        use output_module
        use time_step_module
        use inputs_module, only: spatial_order,cfl,tacc
        use inputs_module, only: output_freq,checkpoint_freq,output_file
        use globals_module, only: elem_h,numFields,numElem
        use solution_module, only: dt
        use c_printer_module
        implicit none

        integer(i4),intent(in) :: timeSteps
        real(dp),intent(inout) :: Q(numFields,numElem)
        real(dp),intent(inout) :: R(numFields,numElem)

        integer(i4) :: timeStep
        integer(i4) :: elem_id
        real(dp) :: t1,t2
        real(dp) :: norm,normQ

        integer(i4) :: vecSize
        integer(i4) :: Rmode,Qmode

        Rmode = 1
        Qmode = 2     
        vecSize = numFields*numElem

        ! ======================================================== !
        call output_solution(output_file,Q)
        call cpu_time(t1)
            do timeStep = 1,timeSteps

                !------ spatial residual ------!
                call rhs(spatial_order,Q,R)
                call L2_norm(numFields,numElem,R,norm)
                call L2_norm(numFields,numElem,Q,normQ)

                !call c_printer(Rmode,R,vecSize,timeStep)
                !call c_printer(Qmode,Q,vecSize,timeStep)

                !------ update solution ------!
                call time_step(tacc,numFields,numElem,Q,elem_h,cfl,dt)
                !dt(1) = 90.0_dp

                if (tacc) then
                    do elem_id = 1,numElem
                        Q(:,elem_id) = Q(:,elem_id) - dt(1)*R(:,elem_id)
                    end do
                else
                    do elem_id = 1,numElem
                        Q(:,elem_id) = Q(:,elem_id) - dt(elem_id)*R(:,elem_id)
                    end do
                end if

                !--------------- Output Solution and Residual ---------------!
                print*,"Time Step: ",timeStep, "L2 Norm: ",norm, "Q norm:", normQ, "dt:",dt(1)
                if (MOD(timeStep,output_freq) == 0) then
                    call cpu_time(t2)
                    print*,'CPU-Time per DOF:', (t2-t1)/timeStep/(numElem)

                    call output_solution(output_file,Q)
                    call output_residual(output_file,R)
                    call writeResNorm(output_file,norm)
                    call nacaSurfacePressure(output_file,Q)
                end if
                !------------------- Checkpoint Solution --------------------!
                if(MOD(timeStep,checkpoint_freq) == 0) call checkpoint(output_file,Q)
            end do
        call cpu_time(t2)
        print*,'CPU-Time per DOF:', (t2-t1)/timesteps/numElem
        ! ======================================================== !

        call output_solution(output_file,Q)
        call output_residual(output_file,R)
        call writeResNorm(output_file,norm)
        call nacaSurfacePressure(output_file,Q)
    end subroutine
end module
