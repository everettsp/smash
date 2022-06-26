!%    This module `mw_optimize` encapsulates all SMASH optimize (type, subroutines, functions)
module mw_optimize
    
    use m_common, only: sp, dp, lchar, np, ns
    use mw_setup, only: SetupDT
    use mw_mesh, only: MeshDT
    use mw_input_data, only: Input_DataDT
    use mw_parameters, only: ParametersDT, &
    & parameters_derived_type_to_matrix, matrix_to_parameters_derived_type
    use mw_states, only: StatesDT
    use mw_output, only: OutputDT
    
    use mw_run, only: direct_model
    
    implicit none
    
    public :: optimize_sbs
    
    private :: transformation
    
    contains
    
        subroutine optimize_sbs(setup, mesh, input_data, parameters, states, output, cost)
        
            implicit none
            
            type(SetupDT), intent(in) :: setup
            type(MeshDT), intent(in) :: mesh
            type(Input_DataDT), intent(in) :: input_data
            type(ParametersDT), intent(inout) :: parameters
            type(StatesDT), intent(inout) :: states
            type(OutputDT), intent(inout) :: output
            real(sp), intent(out) :: cost
            
            type(StatesDT) :: init_states
            integer, dimension(mesh%nrow, mesh%ncol) :: mask_ac
            integer, dimension(2) :: loc_ac
            integer :: iter, nfg, ia, iaa, iam, jf, jfa, jfaa, j, p, pp
            real(sp), dimension(mesh%nrow, mesh%ncol, np) :: parameters_matrix
            real(sp), dimension(np) :: x, x_tf, y, y_tf, z_tf, sdx, lb_tf, ub_tf
            real(sp) :: gx, ga, clg, ddx, dxn, nop, f
            
            where (mesh%global_active_cell .eq. 1)
            
                mask_ac = 1
            
            else where
            
                mask_ac = 0
            
            end where
        
            loc_ac = maxloc(mask_ac)
            
            init_states = states
            
            call direct_model(setup, mesh, input_data, parameters, states, output, cost)
            
            call parameters_derived_type_to_matrix(parameters, parameters_matrix)
            
            x = parameters_matrix(loc_ac(1), loc_ac(2), :)
            
            x_tf = transformation(setup, x)
            lb_tf = transformation(setup, setup%lb_parameters)
            ub_tf = transformation(setup, setup%ub_parameters)

            gx = cost
            ga = gx
            clg = 0.7_sp ** (1._sp / real(np))
            z_tf = x_tf
            sdx = 0._sp
            ddx = 0.64_sp
            dxn = ddx
            ia = 0
            iaa = 0
            iam = 0
            jfaa = 0
            nfg = 0
            
            nop = count(setup%optim_parameters .eq. 1)
            
            do iter=1, int(setup%maxiter * nop + 1)
                
                if (dxn .gt. ddx) dxn = ddx
                if (ddx .gt. 2._sp) ddx = dxn
                
                do p=1, np
                
                    if (setup%optim_parameters(p) .eq. 1) then
                    
                        y_tf = x_tf
                        
                        do 7 j=1, 2
                        
                            jf = 2 * j - 3
                            if (p .eq. iaa .and. jf .eq. -jfaa) goto 7
                            if (x_tf(p) .le. lb_tf(p) .and. jf .lt. 0) goto 7
                            if (x_tf(p) .ge. ub_tf(p) .and. jf .gt. 0) goto 7
                            
                            y_tf(p) = x_tf(p) + jf * ddx
                            if (y_tf(p) .lt. lb_tf(p)) y_tf(p) = lb_tf(p)
                            if (y_tf(p) .gt. ub_tf(p)) y_tf(p) = ub_tf(p)
                            
                            y = inv_transformation(setup, y_tf)
                            
                            do pp=1, np
                            
                                where (mask_ac .eq. 1)
                                
                                    parameters_matrix(:,:,pp) = y(pp)
                                
                                end where
                                
                            end do
                            
                            call matrix_to_parameters_derived_type(parameters_matrix, parameters)
                            
                            states = init_states
                            
                            call direct_model(setup, mesh, input_data, parameters, states, output, cost)
                            
                            f = cost
                            nfg = nfg + 1
                            
                            if (f .lt. gx) then
                            
                                z_tf = y_tf
                                gx = f
                                ia = p
                                jfa = jf
                            
                            end if
                            
                        7 continue

                    end if
                
                end do
                
                iaa = ia
                jfaa = jfa
                
                if (ia .ne. 0) then
                
                    x_tf = z_tf
                    x = inv_transformation(setup, x_tf)
                    
                    do p=1, np
                            
                        where (mask_ac .eq. 1)
                        
                            parameters_matrix(:,:,p) = x(p)
                        
                        end where
                        
                    end do
                    
                    call matrix_to_parameters_derived_type(parameters_matrix, parameters)
                    
                    sdx = clg * sdx
                    sdx(ia) = (1._sp - clg) * real(jfa) * ddx + clg * sdx(ia)
                    
                    iam = iam + 1
                    
                    if (iam .gt. 2 * nop) then
                        
                        ddx = ddx * 2._sp
                        iam = 0
                    
                    end if
                    
                    if (gx .lt. ga - 2) then
                    
                        ga = gx
                    
                    end if
                    
                else
                
                    ddx = ddx / 2._sp
                    iam = 0
                
                end if
                
                if (iter .gt. 4 * nop) then
                
                    do p=1, np
                    
                        if (setup%optim_parameters(p) .eq. 1) then
                        
                            y_tf(p) = x_tf(p) + sdx(p)
                            if (y_tf(p) .lt. lb_tf(p)) y_tf(p) = lb_tf(p)
                            if (y_tf(p) .gt. ub_tf(p)) y_tf(p) = ub_tf(p)
                        
                        end if
                    
                    end do
                    
                    y = inv_transformation(setup, y_tf)
                    
                    do p=1, np
                    
                        where (mask_ac .eq. 1)
                        
                            parameters_matrix(:,:,p) = y(p)
                        
                        end where
                    
                    end do
                    
                    call matrix_to_parameters_derived_type(parameters_matrix, parameters)
                    
                    states = init_states
                    
                    call direct_model(setup, mesh, input_data, parameters, states, output, cost)
                    
                    f = cost
                    
                    nfg = nfg + 1
                    
                    if (f .lt. gx) then
                    
                        gx = f
                        jfaa = 0
                        x_tf = y_tf
                        
                        x = inv_transformation(setup, x_tf)
                        
                        do p=1, np
                        
                            where (mask_ac .eq. 1) 
                            
                                parameters_matrix(:,:,p) = x(p)
                            
                            end where
                        
                        end do
                        
                        call matrix_to_parameters_derived_type(parameters_matrix, parameters)
                        
                        if (gx .lt. ga - 2) then
                        
                            ga = gx
                            
                        end if
                    
                    end if
                
                end if
                
                ia = 0
                
            end do

        end subroutine optimize_sbs
        
        
        function transformation(setup, x) result(x_tf)
        
            implicit none
            
            type(SetupDT), intent(in) :: setup
            real(sp), dimension(np), intent(in) :: x
            real(sp), dimension(np) :: x_tf
            
            integer :: p
            
            do p=1, np
            
                if (setup%lb_parameters(p) .lt. 0._sp) then
                    
                    x_tf(p) = sinh(x(p))
                    
                else if (setup%lb_parameters(p) .ge. 0._sp .and. setup%ub_parameters(p) .le. 1._sp) then
                    
                    x_tf(p) = log(x(p) / (1._sp - x(p)))
                    
                else
                
                    x_tf(p) = log(x(p))
                
                end if
 
            end do

        end function transformation
        
        
        function inv_transformation(setup, x_tf) result(x)
        
            implicit none
            
            type(SetupDT), intent(in) :: setup
            real(sp), dimension(np), intent(in) :: x_tf
            real(sp), dimension(np) :: x
            
            integer :: p
            
            do p=1, np
            
                if (setup%lb_parameters(p) .lt. 0._sp) then
                    
                    x(p) = asinh(x_tf(p))
                    
                else if (setup%lb_parameters(p) .ge. 0._sp .and. setup%ub_parameters(p) .le. 1._sp) then
                    
                    x(p) = exp(x_tf(p)) / (1._sp + exp(x_tf(p)))
                    
                else
                
                    x(p) = exp(x_tf(p))
                
                end if
 
            end do
            
        end function inv_transformation

end module mw_optimize
