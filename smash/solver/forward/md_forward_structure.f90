!%      This module `md_forward_structure` encapsulates all SMASH forward_structure.
!%      This module is differentiated.
!%
!%      contains
!%
!%      [1] gr_a_forward
!%      [2] gr_b_forward
!%      [2] gr_c_forward
!%      [2] gr_d_forward
!%      [3] vic_a_forward

module md_forward_structure

    use md_constant !% only: sp
    use mwd_setup !% only: SetupDT
    use mwd_mesh !% only: MeshDT
    use mwd_input_data !% only: Input_DataDT
    use mwd_parameters !% only: ParametersDT
    use mwd_states !% only: StatesDT
    use mwd_output !% only: OutputDT
    use md_gr_operator !% only: gr_interception, gr_production, gr_exchange, &
    !% & gr_transfer
    use md_vic_operator !% only: vic_infiltration, vic_vertical_transfer, vic_interflow, vic_baseflow
    use md_routing_operator !% only: upstream_discharge, linear_routing

    implicit none

contains

    subroutine gr_a_forward(setup, mesh, input_data, parameters, states, output)

        implicit none

        !% =================================================================================================================== %!
        !%   Derived Type Variables (shared)
        !% =================================================================================================================== %!

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(in) :: parameters
        type(StatesDT), intent(inout) :: states
        type(OutputDT), intent(inout) :: output

        !% =================================================================================================================== %!
        !%   Local Variables (private)
        !% =================================================================================================================== %!
        real(sp), dimension(mesh%nrow, mesh%ncol) :: q
        real(sp) :: prcp, pet, ei, pn, en, pr, perc, l, prr, prd, &
        & qr, qd, qt, qup, qrout
        integer :: t, i, row, col, k, g

        !% =================================================================================================================== %!
        !%   Begin subroutine
        !% =================================================================================================================== %!

        do t = 1, setup%ntime_step !% [ DO TIME ]

            do i = 1, mesh%nrow*mesh%ncol !% [ DO SPACE ]

                !% =============================================================================================================== %!
                !%   Local Variables Initialisation for time step (t) and cell (i)
                !% =============================================================================================================== %!

                ei = 0._sp
                pn = 0._sp
                en = 0._sp
                pr = 0._sp
                perc = 0._sp
                l = 0._sp
                prr = 0._sp
                prd = 0._sp
                qr = 0._sp
                qd = 0._sp
                qup = 0._sp
                qrout = 0._sp

                !% =========================================================================================================== %!
                !%   Cell indice (i) to Cell indices (row, col) following an increasing order of flow accumulation
                !% =========================================================================================================== %!

                if (mesh%path(1, i) .gt. 0 .and. mesh%path(2, i) .gt. 0) then !% [ IF PATH ]

                    row = mesh%path(1, i)
                    col = mesh%path(2, i)
                    if (setup%sparse_storage) k = mesh%rowcol_to_ind_sparse(row, col)

                    !% ======================================================================================================= %!
                    !%   Global/Local active cell
                    !% ======================================================================================================= %!

                    if (mesh%active_cell(row, col) .eq. 1 .and. mesh%local_active_cell(row, col) .eq. 1) then !% [ IF ACTIVE CELL ]

                        if (setup%sparse_storage) then

                            prcp = input_data%sparse_prcp(k, t)
                            pet = input_data%sparse_pet(k, t)

                        else

                            prcp = input_data%prcp(row, col, t)
                            pet = input_data%pet(row, col, t)

                        end if

                        if (prcp .ge. 0 .and. pet .ge. 0) then !% [ IF PRCP GAP ]

                            !% =============================================================================================== %!
                            !%   Interception module
                            !% =============================================================================================== %!

                            ei = min(pet, prcp)

                            pn = max(0._sp, prcp - ei)

                            en = pet - ei

                            !% =============================================================================================== %!
                            !%   Production module
                            !% =============================================================================================== %!

                            call gr_production(pn, en, parameters%cp(row, col), 1000._sp, &
                            & states%hp(row, col), pr, perc)

                            !% =============================================================================================== %!
                            !%   Exchange module
                            !% =============================================================================================== %!

                            call gr_exchange(parameters%exc(row, col), states%hft(row, col), l)

                        end if !% [ END IF PRCP GAP ]

                        !% =================================================================================================== %!
                        !%   Transfer module
                        !% =================================================================================================== %!

                        prr = 0.9_sp*(pr + perc) + l
                        prd = 0.1_sp*(pr + perc)

                        call gr_transfer(5._sp, prcp, prr, parameters%cft(row, col), states%hft(row, col), qr)

                        qd = max(0._sp, prd + l)

                        qt = (qr + qd)

                        !% =================================================================================================== %!
                        !%   Routing module
                        !% =================================================================================================== %!

                        call upstream_discharge(setup%dt, mesh%dx, mesh%nrow,&
                        &  mesh%ncol, mesh%flwdir, mesh%flwacc, row, col, q, qup)

                        call linear_routing(setup%dt, qup, parameters%lr(row, col), states%hlr(row, col), qrout)

                        q(row, col) = (qt + qrout*real(mesh%flwacc(row, col) - 1))&
                                     & *mesh%dx*mesh%dx*0.001_sp/setup%dt

                        !% =================================================================================================== %!
                        !%   Store simulated net rainfall on domain (optional)
                        !%   The net rainfall over a surface is a fictitious quantity that corresponds to
                        !%   the part of the rainfall water depth that actually causes runoff.
                        !% =================================================================================================== %!

                        if (setup%save_net_prcp_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_net_prcp_domain(k, t) = qt

                            else

                                output%net_prcp_domain(row, col, t) = qt

                            end if

                        end if

                        !% =================================================================================================== %!
                        !%   Store simulated discharge on domain (optional)
                        !% =================================================================================================== %!

                        if (setup%save_qsim_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_qsim_domain(k, t) = q(row, col)

                            else

                                output%qsim_domain(row, col, t) = q(row, col)

                            end if

                        end if

                    end if !% [ END IF ACTIVE CELL ]

                end if !% [ END IF PATH ]

            end do !% [ END DO SPACE ]

            !% =============================================================================================================== %!
            !%   Store simulated discharge at gauge
            !% =============================================================================================================== %!

            do g = 1, mesh%ng

                output%qsim(g, t) = q(mesh%gauge_pos(g, 1), mesh%gauge_pos(g, 2))

            end do

        end do !% [ END DO TIME ]

    end subroutine gr_a_forward

    subroutine gr_b_forward(setup, mesh, input_data, parameters, states, output)

        implicit none

        !% =================================================================================================================== %!
        !%   Derived Type Variables (shared)
        !% =================================================================================================================== %!

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(in) :: parameters
        type(StatesDT), intent(inout) :: states
        type(OutputDT), intent(inout) :: output

        !% =================================================================================================================== %!
        !%   Local Variables (private)
        !% =================================================================================================================== %!
        real(sp), dimension(mesh%nrow, mesh%ncol) :: q
        real(sp) :: prcp, pet, ei, pn, en, pr, perc, l, prr, prd, &
        & qr, qd, qt, qup, qrout
        integer :: t, i, row, col, k, g

        !% =================================================================================================================== %!
        !%   Begin subroutine
        !% =================================================================================================================== %!

        do t = 1, setup%ntime_step !% [ DO TIME ]

            do i = 1, mesh%nrow*mesh%ncol !% [ DO SPACE ]

                !% =============================================================================================================== %!
                !%   Local Variables Initialisation for time step (t) and cell (i)
                !% =============================================================================================================== %!

                ei = 0._sp
                pn = 0._sp
                en = 0._sp
                pr = 0._sp
                perc = 0._sp
                l = 0._sp
                prr = 0._sp
                prd = 0._sp
                qr = 0._sp
                qd = 0._sp
                qup = 0._sp
                qrout = 0._sp

                !% =========================================================================================================== %!
                !%   Cell indice (i) to Cell indices (row, col) following an increasing order of flow accumulation
                !% =========================================================================================================== %!

                if (mesh%path(1, i) .gt. 0 .and. mesh%path(2, i) .gt. 0) then !% [ IF PATH ]

                    row = mesh%path(1, i)
                    col = mesh%path(2, i)
                    if (setup%sparse_storage) k = mesh%rowcol_to_ind_sparse(row, col)

                    !% ======================================================================================================= %!
                    !%   Global/Local active cell
                    !% ======================================================================================================= %!

                    if (mesh%active_cell(row, col) .eq. 1 .and. mesh%local_active_cell(row, col) .eq. 1) then !% [ IF ACTIVE CELL ]

                        if (setup%sparse_storage) then

                            prcp = input_data%sparse_prcp(k, t)
                            pet = input_data%sparse_pet(k, t)

                        else

                            prcp = input_data%prcp(row, col, t)
                            pet = input_data%pet(row, col, t)

                        end if

                        if (prcp .ge. 0 .and. pet .ge. 0) then !% [ IF PRCP GAP ]

                            !% =============================================================================================== %!
                            !%   Interception module
                            !% =============================================================================================== %!

                            call gr_interception(prcp, pet, parameters%ci(row, col), states%hi(row, col), pn, ei)

                            en = pet - ei

                            !% =============================================================================================== %!
                            !%   Production module
                            !% =============================================================================================== %!

                            call gr_production(pn, en, parameters%cp(row, col), 1000._sp, &
                            & states%hp(row, col), pr, perc)

                            !% =============================================================================================== %!
                            !%   Exchange module
                            !% =============================================================================================== %!

                            call gr_exchange(parameters%exc(row, col), states%hft(row, col), l)

                        end if !% [ END IF PRCP GAP ]

                        !% =================================================================================================== %!
                        !%   Transfer module
                        !% =================================================================================================== %!

                        prr = 0.9_sp*(pr + perc) + l
                        prd = 0.1_sp*(pr + perc)

                        call gr_transfer(5._sp, prcp, prr, parameters%cft(row, col), states%hft(row, col), qr)

                        qd = max(0._sp, prd + l)

                        qt = (qr + qd)

                        !% =================================================================================================== %!
                        !%   Routing module
                        !% =================================================================================================== %!

                        call upstream_discharge(setup%dt, mesh%dx, mesh%nrow,&
                        &  mesh%ncol, mesh%flwdir, mesh%flwacc, row, col, q, qup)

                        call linear_routing(setup%dt, qup, parameters%lr(row, col), states%hlr(row, col), qrout)

                        q(row, col) = (qt + qrout*real(mesh%flwacc(row, col) - 1))&
                                     & *mesh%dx*mesh%dx*0.001_sp/setup%dt

                        !% =================================================================================================== %!
                        !%   Store simulated net rainfall on domain (optional)
                        !%   The net rainfall over a surface is a fictitious quantity that corresponds to
                        !%   the part of the rainfall water depth that actually causes runoff.
                        !% =================================================================================================== %!

                        if (setup%save_net_prcp_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_net_prcp_domain(k, t) = qt

                            else

                                output%net_prcp_domain(row, col, t) = qt

                            end if

                        end if

                        !% =================================================================================================== %!
                        !%   Store simulated discharge on domain (optional)
                        !% =================================================================================================== %!

                        if (setup%save_qsim_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_qsim_domain(k, t) = q(row, col)

                            else

                                output%qsim_domain(row, col, t) = q(row, col)

                            end if

                        end if

                    end if !% [ END IF ACTIVE CELL ]

                end if !% [ END IF PATH ]

            end do !% [ END DO SPACE ]

            !% =============================================================================================================== %!
            !%   Store simulated discharge at gauge
            !% =============================================================================================================== %!

            do g = 1, mesh%ng

                output%qsim(g, t) = q(mesh%gauge_pos(g, 1), mesh%gauge_pos(g, 2))

            end do

        end do !% [ END DO TIME ]

    end subroutine gr_b_forward

    subroutine gr_c_forward(setup, mesh, input_data, parameters, states, output)

        implicit none

        !% =================================================================================================================== %!
        !%   Derived Type Variables (shared)
        !% =================================================================================================================== %!

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(in) :: parameters
        type(StatesDT), intent(inout) :: states
        type(OutputDT), intent(inout) :: output

        !% =================================================================================================================== %!
        !%   Local Variables (private)
        !% =================================================================================================================== %!
        real(sp), dimension(mesh%nrow, mesh%ncol) :: q
        real(sp) :: prcp, pet, ei, pn, en, pr, perc, l, prr, prl, prd, &
        & qr, ql, qd, qt, qup, qrout
        integer :: t, i, row, col, k, g

        !% =================================================================================================================== %!
        !%   Begin subroutine
        !% =================================================================================================================== %!

        do t = 1, setup%ntime_step !% [ DO TIME ]

            do i = 1, mesh%nrow*mesh%ncol !% [ DO SPACE ]

                !% =============================================================================================================== %!
                !%   Local Variables Initialisation for time step (t) and cell (i)
                !% =============================================================================================================== %!

                ei = 0._sp
                pn = 0._sp
                en = 0._sp
                pr = 0._sp
                perc = 0._sp
                l = 0._sp
                prr = 0._sp
                prl = 0._sp
                prd = 0._sp
                qr = 0._sp
                ql = 0._sp
                qd = 0._sp
                qup = 0._sp
                qrout = 0._sp

                !% =========================================================================================================== %!
                !%   Cell indice (i) to Cell indices (row, col) following an increasing order of flow accumulation
                !% =========================================================================================================== %!

                if (mesh%path(1, i) .gt. 0 .and. mesh%path(2, i) .gt. 0) then !% [ IF PATH ]

                    row = mesh%path(1, i)
                    col = mesh%path(2, i)
                    if (setup%sparse_storage) k = mesh%rowcol_to_ind_sparse(row, col)

                    !% ======================================================================================================= %!
                    !%   Global/Local active cell
                    !% ======================================================================================================= %!

                    if (mesh%active_cell(row, col) .eq. 1 .and. mesh%local_active_cell(row, col) .eq. 1) then !% [ IF ACTIVE CELL ]

                        if (setup%sparse_storage) then

                            prcp = input_data%sparse_prcp(k, t)
                            pet = input_data%sparse_pet(k, t)

                        else

                            prcp = input_data%prcp(row, col, t)
                            pet = input_data%pet(row, col, t)

                        end if

                        if (prcp .ge. 0 .and. pet .ge. 0) then !% [ IF PRCP GAP ]

                            !% =============================================================================================== %!
                            !%   Interception module
                            !% =============================================================================================== %!

                            call gr_interception(prcp, pet, parameters%ci(row, col), states%hi(row, col), pn, ei)

                            en = pet - ei

                            !% =============================================================================================== %!
                            !%   Production module
                            !% =============================================================================================== %!

                            call gr_production(pn, en, parameters%cp(row, col), 1000._sp, &
                            & states%hp(row, col), pr, perc)

                            !% =============================================================================================== %!
                            !%   Exchange module
                            !% =============================================================================================== %!

                            call gr_exchange(parameters%exc(row, col), states%hft(row, col), l)

                        end if !% [ END IF PRCP GAP ]

                        !% =================================================================================================== %!
                        !%   Transfer module
                        !% =================================================================================================== %!

                        prr = 0.9_sp*0.6_sp*(pr + perc) + l
                        prl = 0.9_sp*0.4_sp*(pr + perc)
                        prd = 0.1_sp*(pr + perc)

                        call gr_transfer(5._sp, prcp, prr, parameters%cft(row, col), states%hft(row, col), qr)

                        call gr_transfer(5._sp, prcp, prl, parameters%cst(row, col), states%hst(row, col), ql)

                        qd = max(0._sp, prd + l)

                        qt = (qr + ql + qd)

                        !% =================================================================================================== %!
                        !%   Routing module
                        !% =================================================================================================== %!

                        call upstream_discharge(setup%dt, mesh%dx, mesh%nrow,&
                        &  mesh%ncol, mesh%flwdir, mesh%flwacc, row, col, q, qup)

                        call linear_routing(setup%dt, qup, parameters%lr(row, col), states%hlr(row, col), qrout)

                        q(row, col) = (qt + qrout*real(mesh%flwacc(row, col) - 1))&
                                     & *mesh%dx*mesh%dx*0.001_sp/setup%dt

                        !% =================================================================================================== %!
                        !%   Store simulated net rainfall on domain (optional)
                        !%   The net rainfall over a surface is a fictitious quantity that corresponds to
                        !%   the part of the rainfall water depth that actually causes runoff.
                        !% =================================================================================================== %!

                        if (setup%save_net_prcp_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_net_prcp_domain(k, t) = qt

                            else

                                output%net_prcp_domain(row, col, t) = qt

                            end if

                        end if

                        !% =================================================================================================== %!
                        !%   Store simulated discharge on domain (optional)
                        !% =================================================================================================== %!

                        if (setup%save_qsim_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_qsim_domain(k, t) = q(row, col)

                            else

                                output%qsim_domain(row, col, t) = q(row, col)

                            end if

                        end if

                    end if !% [ END IF ACTIVE CELL ]

                end if !% [ END IF PATH ]

            end do !% [ END DO SPACE ]

            !% =============================================================================================================== %!
            !%   Store simulated discharge at gauge
            !% =============================================================================================================== %!

            do g = 1, mesh%ng

                output%qsim(g, t) = q(mesh%gauge_pos(g, 1), mesh%gauge_pos(g, 2))

            end do

        end do !% [ END DO TIME ]

    end subroutine gr_c_forward

    subroutine gr_d_forward(setup, mesh, input_data, parameters, states, output)

        implicit none

        !% =================================================================================================================== %!
        !%   Derived Type Variables (shared)
        !% =================================================================================================================== %!

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(in) :: parameters
        type(StatesDT), intent(inout) :: states
        type(OutputDT), intent(inout) :: output

        !% =================================================================================================================== %!
        !%   Local Variables (private)
        !% =================================================================================================================== %!
        real(sp), dimension(mesh%nrow, mesh%ncol) :: q
        real(sp) :: prcp, pet, ei, pn, en, pr, perc, prr, qr, qt, qup, qrout
        integer :: t, i, row, col, k, g

        !% =================================================================================================================== %!
        !%   Begin subroutine
        !% =================================================================================================================== %!

        do t = 1, setup%ntime_step !% [ DO TIME ]

            do i = 1, mesh%nrow*mesh%ncol !% [ DO SPACE ]

                !% =============================================================================================================== %!
                !%   Local Variables Initialisation for time step (t) and cell (i)
                !% =============================================================================================================== %!

                ei = 0._sp
                pn = 0._sp
                en = 0._sp
                pr = 0._sp
                perc = 0._sp
                prr = 0._sp
                qr = 0._sp
                qup = 0._sp
                qrout = 0._sp

                !% =========================================================================================================== %!
                !%   Cell indice (i) to Cell indices (row, col) following an increasing order of flow accumulation
                !% =========================================================================================================== %!

                if (mesh%path(1, i) .gt. 0 .and. mesh%path(2, i) .gt. 0) then !% [ IF PATH ]

                    row = mesh%path(1, i)
                    col = mesh%path(2, i)
                    if (setup%sparse_storage) k = mesh%rowcol_to_ind_sparse(row, col)

                    !% ======================================================================================================= %!
                    !%   Global/Local active cell
                    !% ======================================================================================================= %!

                    if (mesh%active_cell(row, col) .eq. 1 .and. mesh%local_active_cell(row, col) .eq. 1) then !% [ IF ACTIVE CELL ]

                        if (setup%sparse_storage) then

                            prcp = input_data%sparse_prcp(k, t)
                            pet = input_data%sparse_pet(k, t)

                        else

                            prcp = input_data%prcp(row, col, t)
                            pet = input_data%pet(row, col, t)

                        end if

                        if (prcp .ge. 0 .and. pet .ge. 0) then !% [ IF PRCP GAP ]

                            !% =============================================================================================== %!
                            !%   Interception module
                            !% =============================================================================================== %!

                            ei = min(pet, prcp)

                            pn = max(0._sp, prcp - ei)

                            en = pet - ei

                            !% =============================================================================================== %!
                            !%   Production module
                            !% =============================================================================================== %!

                            call gr_production(pn, en, parameters%cp(row, col), 1000._sp, &
                            & states%hp(row, col), pr, perc)

                        end if !% [ END IF PRCP GAP ]

                        !% =================================================================================================== %!
                        !%   Transfer module
                        !% =================================================================================================== %!

                        prr = pr + perc

                        call gr_transfer(5._sp, prcp, prr, parameters%cft(row, col), states%hft(row, col), qr)

                        qt = qr

                        !% =================================================================================================== %!
                        !%   Routing module
                        !% =================================================================================================== %!

                        call upstream_discharge(setup%dt, mesh%dx, mesh%nrow,&
                        &  mesh%ncol, mesh%flwdir, mesh%flwacc, row, col, q, qup)

                        call linear_routing(setup%dt, qup, parameters%lr(row, col), states%hlr(row, col), qrout)

                        q(row, col) = (qt + qrout*real(mesh%flwacc(row, col) - 1))&
                                     & *mesh%dx*mesh%dx*0.001_sp/setup%dt

                        !% =================================================================================================== %!
                        !%   Store simulated net rainfall on domain (optional)
                        !%   The net rainfall over a surface is a fictitious quantity that corresponds to
                        !%   the part of the rainfall water depth that actually causes runoff.
                        !% =================================================================================================== %!

                        if (setup%save_net_prcp_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_net_prcp_domain(k, t) = qt

                            else

                                output%net_prcp_domain(row, col, t) = qt

                            end if

                        end if

                        !% =================================================================================================== %!
                        !%   Store simulated discharge on domain (optional)
                        !% =================================================================================================== %!

                        if (setup%save_qsim_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_qsim_domain(k, t) = q(row, col)

                            else

                                output%qsim_domain(row, col, t) = q(row, col)

                            end if

                        end if

                    end if !% [ END IF ACTIVE CELL ]

                end if !% [ END IF PATH ]

            end do !% [ END DO SPACE ]

            !% =============================================================================================================== %!
            !%   Store simulated discharge at gauge
            !% =============================================================================================================== %!

            do g = 1, mesh%ng

                output%qsim(g, t) = q(mesh%gauge_pos(g, 1), mesh%gauge_pos(g, 2))

            end do

        end do !% [ END DO TIME ]

    end subroutine gr_d_forward



    subroutine gr_g_forward(setup, mesh, input_data, parameters, states, output)

        implicit none

        !% =================================================================================================================== %!
        !%   Derived Type Variables (shared)
        !% =================================================================================================================== %!

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(in) :: parameters
        type(StatesDT), intent(inout) :: states
        type(OutputDT), intent(inout) :: output

        !% =================================================================================================================== %!
        !%   Local Variables (private)
        !% =================================================================================================================== %!
        real(sp), dimension(mesh%nrow, mesh%ncol) :: q
        real(sp) :: prcp, pet, ei, pn, en, pr, perc, l, prr, prd, &
        & qr, qd, qt, qup, qrout
        integer :: t, i, row, col, k, g

        !% =================================================================================================================== %!
        !%   Begin subroutine
        !% =================================================================================================================== %!

        do t = 1, setup%ntime_step !% [ DO TIME ]

            do i = 1, mesh%nrow*mesh%ncol !% [ DO SPACE ]

                !% =============================================================================================================== %!
                !%   Local Variables Initialisation for time step (t) and cell (i)
                !% =============================================================================================================== %!

                ei = 0._sp
                pn = 0._sp
                en = 0._sp
                pr = 0._sp
                perc = 0._sp
                l = 0._sp
                prr = 0._sp
                prd = 0._sp
                qr = 0._sp
                qd = 0._sp
                qup = 0._sp
                qrout = 0._sp

                !% =========================================================================================================== %!
                !%   Cell indice (i) to Cell indices (row, col) following an increasing order of flow accumulation
                !% =========================================================================================================== %!

                if (mesh%path(1, i) .gt. 0 .and. mesh%path(2, i) .gt. 0) then !% [ IF PATH ]

                    row = mesh%path(1, i)
                    col = mesh%path(2, i)
                    if (setup%sparse_storage) k = mesh%rowcol_to_ind_sparse(row, col)

                    !% ======================================================================================================= %!
                    !%   Global/Local active cell
                    !% ======================================================================================================= %!

                    if (mesh%active_cell(row, col) .eq. 1 .and. mesh%local_active_cell(row, col) .eq. 1) then !% [ IF ACTIVE CELL ]

                        if (setup%sparse_storage) then

                            prcp = input_data%sparse_prcp(k, t)
                            pet = input_data%sparse_pet(k, t)

                        else

                            prcp = input_data%prcp(row, col, t)
                            pet = input_data%pet(row, col, t)

                        end if

                        if (prcp .ge. 0 .and. pet .ge. 0) then !% [ IF PRCP GAP ]

                            !% =============================================================================================== %!
                            !%   Interception module
                            !% =============================================================================================== %!

                            call gr_interception(prcp, pet, parameters%ci(row, col), states%hi(row, col), pn, ei)

                            en = pet - ei

                            !% =============================================================================================== %!
                            !%   Production module
                            !% =============================================================================================== %!

                            call gr_production(pn, en, parameters%cp(row, col), 1000._sp, &
                            & states%hp(row, col), pr, perc)

                            !% =============================================================================================== %!
                            !%   Exchange module
                            !% =============================================================================================== %!

                            call gr_exchange(parameters%exc(row, col), states%hft(row, col), l)

                        end if !% [ END IF PRCP GAP ]

                        !% =================================================================================================== %!
                        !%   Transfer module
                        !% =================================================================================================== %!

                        prr = 0.9_sp*(pr + perc) + l
                        prd = 0.1_sp*(pr + perc)

                        call gr_transfer(5._sp, prcp, prr, parameters%cft(row, col), states%hft(row, col), qr)

                        qd = max(0._sp, prd + l)

                        qt = (qr + qd)

                        q(row, col) = qt*mesh%dx*mesh%dx*0.001_sp/setup%dt

                        !% =================================================================================================== %!
                        !%   Store simulated net rainfall on domain (optional)
                        !%   The net rainfall over a surface is a fictitious quantity that corresponds to
                        !%   the part of the rainfall water depth that actually causes runoff.
                        !% =================================================================================================== %!

                        if (setup%save_net_prcp_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_net_prcp_domain(k, t) = qt

                            else

                                output%net_prcp_domain(row, col, t) = qt

                            end if

                        end if

                        !% =================================================================================================== %!
                        !%   Store simulated discharge on domain (optional)
                        !% =================================================================================================== %!

                        if (setup%save_qsim_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_qsim_domain(k, t) = q(row, col)

                            else

                                output%qsim_domain(row, col, t) = q(row, col)

                            end if

                        end if

                    end if !% [ END IF ACTIVE CELL ]

                end if !% [ END IF PATH ]

            end do !% [ END DO SPACE ]

            !% =============================================================================================================== %!
            !%   Store simulated discharge at gauge
            !% =============================================================================================================== %!

            do g = 1, mesh%ng

                output%qsim(g, t) = q(mesh%gauge_pos(g, 1), mesh%gauge_pos(g, 2))

            end do

        end do !% [ END DO TIME ]

        ! Fortran coupled Gamma Model Here ?
        ! Compilation ?
        
        ! use mod_gamma_routing_setup
        ! use mod_gamma_routing_mesh
        ! use mod_gamma_routing_parameters
        ! use mod_gamma_routing_states
        ! use mod_gamma_routing_results
        ! use mod_gamma_function
        ! use mod_gamma_interface
        ! use mod_gamma_routing
        
        !INPUTS MUST BE READ OR SET SOMEWHERE :
        !In a new derved type setup%coupling ?
        ! In a new deived type parameter%coupling
        !dt
        !hydraulics_coefficient
        !spreading
        !vmin=0.1,vmax=10.,&
        !&elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=0.1,ponderation_regul=10000.0&
        !&,velocity_computation="qm3",varying_spread=.true.

!~ FROM HERE TO ..

!~        npdt=int(setup._ntime_step*setup.dt/dt)

!~        call routing_setup_self_initialisation(routing_setup,npdt=npdt,dt=900.,vmin=0.1,vmax=10.,&
!~        &elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=0.1,ponderation_regul=10000.0&
!~        &,velocity_computation="qm3",varying_spread=.true.)
        
!~        nb_nodes=np.sum(mesh.active_cell)
!~        nb_upstream_nodes=8 #properties of gridded mesh
!~        routing_mesh_self_initialisation(nb_nodes=nb_nodes,nb_upstream_nodes=nb_upstream_nodes) #initialise the gamma mesh
        
           !Here 3 functions written in python must be translated in Fortran
!~         VectorMesh1D=SmashMesh2DToVector(smash_model)
!~         nodes_linker=ComputeNodeLinker(smash_model)
!~         gauge_pos=VectorIndexesOfSmashGauges(smash_model)
        
!~         #Filling the Gamma mesh properties
!~         routing_mesh%upstream_to_downstream_nodes=np.argsort(VectorMesh1D["FlowAcc"])+1
!~         routing_mesh%nodes_linker=nodes_linker
!~         routing_mesh%controlled_nodes=gauge_pos+1
!~         routing_mesh%surface=(smash_model.mesh.dx)**2. / 1000.0**2.  #!km²
!~         routing_mesh%dx=smash_model.mesh.dx
        
!~         #Update and finalise the gamma mesh for the model_gamma object
!~         call mesh_update(routing_mesh)

!~         call routing_parameter_self_initialisation(routing_parameter=routing_parameter,routing_setup=routing_setup,&
!~          &routing_mesh=routing_mesh,&
!~          &hydraulics_coefficient=0.5,spreading=1.)
!~         call routing_state_self_initialisation(routing_setup,routing_mesh,routing_parameter,routing_states)
!~         call compute_gamma_parameters(routing_setup,routing_mesh,routing_states)
!~         call routing_results_self_initialisation(routing_setup,routing_mesh,routing_results)
!~         write(*,*) "routing_hydrogram_forward..."
!~         write(*,*) ""
!~         allocate(observations(routing_setup%npdt,routing_mesh%nb_nodes))
!~         observations=0.0
!~         call routing_states_reset(routing_states)

!~ .. TO HERE : DO NOT DIFFERENTIATE WITH TAPENADE

!~ HOW TO MANAGE THE DIFFERENTIATION WITH RESPECT TO ROUTING_PARAMETER ? => New derived type parameters%coupling ?

!~         call routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,inflows,observations,&
!~         &routing_states,routing_results,cost)

!A ajuster
!~     do g = 1, routing_results%controlled_nodes

!~             output%qsim(g, t) = routing_results%discharges(routing_results%controlled_nodes(g, 1), routing_results%controlled_nodes(g, 2))

!~     end do

    end subroutine gr_g_forward



    subroutine vic_a_forward(setup, mesh, input_data, parameters, states, output)

        implicit none

        !% =================================================================================================================== %!
        !%   Derived Type Variables (shared)
        !% =================================================================================================================== %!

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(in) :: parameters
        type(StatesDT), intent(inout) :: states
        type(OutputDT), intent(inout) :: output

        !% =================================================================================================================== %!
        !%   Local Variables (private)
        !% =================================================================================================================== %!

        real(sp), dimension(mesh%nrow, mesh%ncol) :: q
        real(sp) :: prcp, pet, runoff, qi, qb, qt, qup, qrout
        integer :: t, i, row, col, k, g

        !% =================================================================================================================== %!
        !%   Begin subroutine
        !% =================================================================================================================== %!

        do t = 1, setup%ntime_step !% [ DO TIME ]

            do i = 1, mesh%nrow*mesh%ncol !% [ DO SPACE ]

                !% =============================================================================================================== %!
                !%   Local Variables Initialisation for time step (t) and cell (i)
                !% =============================================================================================================== %!

                runoff = 0._sp
                qi = 0._sp
                qb = 0._sp
                qt = 0._sp
                qup = 0._sp
                qrout = 0._sp

                !% =========================================================================================================== %!
                !%   Cell indice (i) to Cell indices (row, col) following an increasing order of flow accumulation
                !% =========================================================================================================== %!

                if (mesh%path(1, i) .gt. 0 .and. mesh%path(2, i) .gt. 0) then !% [ IF PATH ]

                    row = mesh%path(1, i)
                    col = mesh%path(2, i)
                    if (setup%sparse_storage) k = mesh%rowcol_to_ind_sparse(row, col)

                    !% ======================================================================================================= %!
                    !%   Global/Local active cell
                    !% ======================================================================================================= %!

                    if (mesh%active_cell(row, col) .eq. 1 .and. mesh%local_active_cell(row, col) .eq. 1) then !% [ IF ACTIVE CELL ]

                        if (setup%sparse_storage) then

                            prcp = input_data%sparse_prcp(k, t)
                            pet = input_data%sparse_pet(k, t)

                        else

                            prcp = input_data%prcp(row, col, t)
                            pet = input_data%pet(row, col, t)

                        end if

                        if (prcp .ge. 0 .and. pet .ge. 0) then !% [ IF PRCP GAP ]

                            !% =============================================================================================== %!
                            !%   Infiltration module
                            !% =============================================================================================== %!

                            call vic_infiltration(prcp, parameters%cusl1(row, col), parameters%cusl2(row, col), &
                            & parameters%b(row, col), states%husl1(row, col), states%husl2(row, col), &
                            & runoff)

                            !% =============================================================================================== %!
                            !%   Vertical transfer module
                            !% =============================================================================================== %!

                            call vic_vertical_transfer(pet, parameters%cusl1(row, col), parameters%cusl2(row, col), &
                            & parameters%clsl(row, col), parameters%ks(row, col), states%husl1(row, col), &
                            & states%husl2(row, col), states%hlsl(row, col))

                        end if !% [ END IF PRCP GAP ]

                        !% =================================================================================================== %!
                        !%   Horizontal transfer module
                        !% =================================================================================================== %!

                        call vic_interflow(5._sp, parameters%cusl2(row, col), states%husl2(row, col), qi)

                        call vic_baseflow(parameters%clsl(row, col), parameters%ds(row, col), &
                        & parameters%dsm(row, col), parameters%ws(row, col), states%hlsl(row, col), qb)

                        qt = (runoff + qi + qb)

                        !% =================================================================================================== %!
                        !%   Routing module
                        !% =================================================================================================== %!

                        call upstream_discharge(setup%dt, mesh%dx, mesh%nrow,&
                        &  mesh%ncol, mesh%flwdir, mesh%flwacc, row, col, q, qup)

                        call linear_routing(setup%dt, qup, parameters%lr(row, col), states%hlr(row, col), qrout)

                        q(row, col) = (qt + qrout*real(mesh%flwacc(row, col) - 1))&
                                     & *mesh%dx*mesh%dx*0.001_sp/setup%dt

                        !% =================================================================================================== %!
                        !%   Store simulated net rainfall on domain (optional)
                        !%   The net rainfall over a surface is a fictitious quantity that corresponds to
                        !%   the part of the rainfall water depth that actually causes runoff.
                        !% =================================================================================================== %!

                        if (setup%save_net_prcp_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_net_prcp_domain(k, t) = qt

                            else

                                output%net_prcp_domain(row, col, t) = qt

                            end if

                        end if

                        !% =================================================================================================== %!
                        !%   Store simulated discharge on domain (optional)
                        !% =================================================================================================== %!

                        if (setup%save_qsim_domain) then

                            if (setup%sparse_storage) then

                                output%sparse_qsim_domain(k, t) = q(row, col)

                            else

                                output%qsim_domain(row, col, t) = q(row, col)

                            end if

                        end if

                    end if !% [ END IF ACTIVE CELL ]

                end if !% [ END IF PATH ]

            end do !% [ END DO SPACE ]

            !% =============================================================================================================== %!
            !%   Store simulated discharge at gauge
            !% =============================================================================================================== %!

            do g = 1, mesh%ng

                output%qsim(g, t) = q(mesh%gauge_pos(g, 1), mesh%gauge_pos(g, 2))

            end do

        end do !% [ END DO TIME ]

    end subroutine vic_a_forward

end module md_forward_structure
