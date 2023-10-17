!%      This module `mw_forward` encapsulates all SMASH routine.
!%      This module is wrapped

module mw_forward

    use md_constant, only: sp, dp
    use mwd_setup, only: SetupDT
    use mwd_mesh, only: MeshDT
    use mwd_input_data, only: Input_DataDT
    use mwd_parameters, only: ParametersDT_initialise, ParametersDT, Hyper_ParametersDT
    use mwd_states, only: StatesDT_initialise, StatesDT, Hyper_StatesDT
    use mwd_output, only: OutputDT_initialise, OutputDT
    use mwd_parameters_manipulation, only: normalize_parameters
    use mwd_states_manipulation, only: normalize_states

    implicit none

contains

    subroutine forward(setup, mesh, input_data, parameters, &
    & parameters_bgd, states, states_bgd, output, cost)

        !% Notes
        !% -----
        !%
        !% Wrapping base_forward

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters, parameters_bgd
        type(StatesDT), intent(inout) :: states, states_bgd
        type(OutputDT), intent(inout) :: output
        real(sp), intent(inout) :: cost

        call base_forward(setup, mesh, input_data, parameters, &
        & parameters_bgd, states, states_bgd, output, cost)

    end subroutine forward

    subroutine forward_b(setup, mesh, input_data, parameters, &
    & parameters_b, parameters_bgd, parameters_bgd_b, states, &
    & states_b, states_bgd, states_bgd_b, output, output_b, &
    & cost, cost_b)

        !% Notes
        !% -----
        !%
        !% Wrapping base_forward_b
        !% Gradient of useful results: cost
        !% with respect to varying inputs: *(parameters.ci) *(parameters.cp)
        !%                *(parameters.beta) *(parameters.cft) *(parameters.cst)
        !%                *(parameters.alpha) *(parameters.exc) *(parameters.b)
        !%                *(parameters.cusl1) *(parameters.cusl2) *(parameters.clsl)
        !%                *(parameters.ks) *(parameters.ds) *(parameters.dsm)
        !%                *(parameters.ws) *(parameters.lr) *(states.hi)
        !%                *(states.hp) *(states.hft) *(states.hst) *(states.husl1)
        !%                *(states.husl2) *(states.hlsl) *(states.hlr)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters, parameters_b, &
        & parameters_bgd, parameters_bgd_b
        type(StatesDT), intent(inout) :: states, states_b, &
        & states_bgd, states_bgd_b
        type(OutputDT), intent(inout) :: output, output_b
        real(sp), intent(inout) :: cost, cost_b

        call base_forward_b(setup, mesh, input_data, parameters, &
        & parameters_b, parameters_bgd, parameters_bgd_b, states, &
        & states_b, states_bgd, states_bgd_b, output, output_b, &
        & cost, cost_b)

    end subroutine forward_b


    subroutine wrapped_forward_b(setup, mesh, input_data, parameters, parameters_b, &
        &  states, states_b, output, cost)

        !% Notes
        !% -----
        !%
        !% Wrapping base_forward_b
        !% Gradient of useful results: cost
        !% with respect to varying inputs: *(parameters.ci) *(parameters.cp)
        !%                *(parameters.beta) *(parameters.cft) *(parameters.cst)
        !%                *(parameters.alpha) *(parameters.exc) *(parameters.b)
        !%                *(parameters.cusl1) *(parameters.cusl2) *(parameters.clsl)
        !%                *(parameters.ks) *(parameters.ds) *(parameters.dsm)
        !%                *(parameters.ws) *(parameters.lr) *(states.hi)
        !%                *(states.hp) *(states.hft) *(states.hst) *(states.husl1)
        !%                *(states.husl2) *(states.hlsl) *(states.hlr)

        implicit none

        type(SetupDT), intent(inout) :: setup
        type(MeshDT), intent(inout) :: mesh
        type(Input_DataDT), intent(inout) :: input_data
        type(ParametersDT), intent(inout) :: parameters
        type(ParametersDT), intent(inout) :: parameters_b
        type(StatesDT), intent(inout) :: states
        type(StatesDT), intent(inout) :: states_b
        type(OutputDT), intent(inout) :: output
        real(sp), intent(inout) :: cost 

        type(ParametersDT) :: parameters_bgd, parameters_bgd_b
        type(StatesDT) :: states_bgd, states_bgd_b
        type(OutputDT) :: output_b
        real(sp) :: cost_b
       
        
        call ParametersDT_initialise(parameters_bgd_b, mesh)
        call ParametersDT_initialise(parameters_bgd, mesh)

        call StatesDT_initialise(states_bgd_b, mesh)
        call StatesDT_initialise(states_bgd, mesh)

        call OutputDT_initialise(output_b, setup, mesh)

        call normalize_parameters(setup, mesh, parameters)
        call normalize_states(setup, mesh, states)
        
        ! Background is normalize
        parameters_bgd = parameters
        states_bgd = states
        
        cost_b = 1._sp
        cost = 0._sp
        
        call base_forward_b(setup, mesh, input_data, parameters, &
        & parameters_b, parameters_bgd, parameters_bgd_b, states, &
        & states_b, states_bgd, states_bgd_b, output, output_b, &
        & cost, cost_b)
            
            
    end subroutine wrapped_forward_b


    subroutine forward_b0(setup, mesh, input_data, parameters, &
    & parameters_b, parameters_bgd, parameters_bgd_b, states, &
    & states_b, states_bgd, states_bgd_b, output, output_b, &
    & cost)

        !% Notes
        !% -----
        !%
        !% Wrapping base_forward_b0 
        !% Gradient of useful results: *(output.qsim_domain)
        !% with respect to varying inputs: *(parameters.ci) *(parameters.cp)
        !%                *(parameters.beta) *(parameters.cft) *(parameters.cst)
        !%                *(parameters.alpha) *(parameters.exc) *(parameters.b)
        !%                *(parameters.cusl1) *(parameters.cusl2) *(parameters.clsl)
        !%                *(parameters.ks) *(parameters.ds) *(parameters.dsm)
        !%                *(parameters.ws) *(parameters.lr) *(states.hi)
        !%                *(states.hp) *(states.hft) *(states.hst) *(states.husl1)
        !%                *(states.husl2) *(states.hlsl) *(states.hlr)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters, parameters_b, &
        & parameters_bgd, parameters_bgd_b
        type(StatesDT), intent(inout) :: states, states_b, &
        & states_bgd, states_bgd_b
        type(OutputDT), intent(inout) :: output, output_b
        real(sp), intent(inout) :: cost

        call base_forward_b0(setup, mesh, input_data, parameters, &
        & parameters_b, parameters_bgd, parameters_bgd_b, states, &
        & states_b, states_bgd, states_bgd_b, output, output_b, &
        & cost)
        

    end subroutine forward_b0



    subroutine wrapped_forward_b0(setup, mesh, input_data, parameters, parameters_b, &
        &  states, states_b, output, cost)

        !% Notes
        !% -----
        !%
        !% Wrapping base_forward_b
        !% Gradient of useful results: cost
        !% with respect to varying inputs: *(parameters.ci) *(parameters.cp)
        !%                *(parameters.beta) *(parameters.cft) *(parameters.cst)
        !%                *(parameters.alpha) *(parameters.exc) *(parameters.b)
        !%                *(parameters.cusl1) *(parameters.cusl2) *(parameters.clsl)
        !%                *(parameters.ks) *(parameters.ds) *(parameters.dsm)
        !%                *(parameters.ws) *(parameters.lr) *(states.hi)
        !%                *(states.hp) *(states.hft) *(states.hst) *(states.husl1)
        !%                *(states.husl2) *(states.hlsl) *(states.hlr)

        implicit none

        type(SetupDT), intent(inout) :: setup
        type(MeshDT), intent(inout) :: mesh
        type(Input_DataDT), intent(inout) :: input_data
        type(ParametersDT), intent(inout) :: parameters
        type(ParametersDT), intent(inout) :: parameters_b
        type(StatesDT), intent(inout) :: states
        type(StatesDT), intent(inout) :: states_b
        type(OutputDT), intent(inout) :: output
        real(sp), intent(inout) :: cost

        type(ParametersDT) :: parameters_bgd, parameters_bgd_b
        type(StatesDT) :: states_bgd, states_bgd_b
        type(OutputDT) :: output_b
        
        
        call ParametersDT_initialise(parameters_bgd_b, mesh)
        call ParametersDT_initialise(parameters_bgd, mesh)

        call StatesDT_initialise(states_bgd_b, mesh)
        call StatesDT_initialise(states_bgd, mesh)

        call OutputDT_initialise(output_b, setup, mesh)

        call normalize_parameters(setup, mesh, parameters)
        call normalize_states(setup, mesh, states)
        
        ! Trigger the denormalization subroutine in forward
        setup%optimize%denormalize_forward = .true.
        
        ! Background is normalize
        parameters_bgd = parameters
        states_bgd = states
        
        cost = 0._sp
        output_b%qsim_domain=1._sp
        
        call base_forward_b0(setup, mesh, input_data, parameters, &
        & parameters_b, parameters_bgd, parameters_bgd_b, states, &
        & states_b, states_bgd, states_bgd_b, output, output_b, &
        & cost)
        
        ! Remove the denormalization subroutine in forward
        setup%optimize%denormalize_forward = .false.
        
        end subroutine wrapped_forward_b0


    subroutine forward_d(setup, mesh, input_data, parameters, &
    & parameters_d, parameters_bgd, parameters_bgd_d, states, &
    & states_d, states_bgd, states_bgd_d, output, output_d, &
    & cost, cost_d)

        !% Notes
        !% -----
        !%
        !% Wrapping base_forward_d

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters, parameters_d, &
        & parameters_bgd, parameters_bgd_d
        type(StatesDT), intent(inout) :: states, states_d, &
        & states_bgd, states_bgd_d
        type(OutputDT), intent(inout) :: output, output_d
        real(sp), intent(inout) :: cost, cost_d

        call base_forward_d(setup, mesh, input_data, parameters, &
        & parameters_d, parameters_bgd, parameters_bgd_d, states, &
        & states_d, states_bgd, states_bgd_d, output, output_d, &
        & cost, cost_d)

    end subroutine forward_d

    subroutine hyper_forward(setup, mesh, input_data, parameters, &
    & hyper_parameters, hyper_parameters_bgd, states, hyper_states, &
    & hyper_states_bgd, output, cost)

        !% Notes
        !% -----
        !%
        !% Wrapping base_hyper_forward

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters
        type(Hyper_ParametersDT), intent(inout) :: hyper_parameters, hyper_parameters_bgd
        type(StatesDT), intent(inout) :: states
        type(Hyper_StatesDT), intent(inout) :: hyper_states, hyper_states_bgd
        type(OutputDT), intent(inout) :: output
        real(sp), intent(inout) :: cost

        call base_hyper_forward(setup, mesh, input_data, parameters, hyper_parameters, &
        & hyper_parameters_bgd, states, hyper_states, hyper_states_bgd, output, cost)

    end subroutine hyper_forward

    subroutine hyper_forward_b(setup, mesh, input_data, &
    & parameters, parameters_b, hyper_parameters, hyper_parameters_b, &
    & hyper_parameters_bgd, states, states_b, hyper_states, hyper_states_b, &
    & hyper_states_bgd, output, output_b, cost, cost_b)

        !% Notes
        !% -----
        !%
        !% Wrapping base_hyper_forward_b

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters, parameters_b
        type(Hyper_ParametersDT), intent(inout) :: hyper_parameters, hyper_parameters_b, hyper_parameters_bgd
        type(StatesDT), intent(inout) :: states, states_b
        type(Hyper_StatesDT), intent(inout) :: hyper_states, hyper_states_b, hyper_states_bgd
        type(OutputDT), intent(inout) :: output, output_b
        real(sp), intent(inout) :: cost, cost_b

        call base_hyper_forward_b(setup, mesh, input_data, parameters, &
        & parameters_b, hyper_parameters, hyper_parameters_b, &
        & hyper_parameters_bgd, states, states_b, hyper_states, hyper_states_b, &
        & hyper_states_bgd, output, output_b, cost, cost_b)

    end subroutine hyper_forward_b

    subroutine hyper_forward_d(setup, mesh, input_data, &
    & parameters, parameters_d, hyper_parameters, hyper_parameters_d, hyper_parameters_bgd, &
    & states, states_d, hyper_states, hyper_states_d, hyper_states_bgd, &
    & output, output_d, cost, cost_d)

        !% Notes
        !% -----
        !%
        !% Wrapping base_hyper_forward_d

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters, parameters_d
        type(Hyper_ParametersDT), intent(inout) :: hyper_parameters, hyper_parameters_d, hyper_parameters_bgd
        type(StatesDT), intent(inout) :: states, states_d
        type(Hyper_StatesDT), intent(inout) :: hyper_states, hyper_states_d, hyper_states_bgd
        type(OutputDT), intent(inout) :: output, output_d
        real(sp), intent(inout) :: cost, cost_d

        call base_hyper_forward_d(setup, mesh, input_data, parameters, &
        & parameters_d, hyper_parameters, hyper_parameters_d, &
        & hyper_parameters_bgd, states, states_d, hyper_states, hyper_states_d, &
        & hyper_states_bgd, output, output_d, cost, cost_d)

    end subroutine hyper_forward_d

end module mw_forward
