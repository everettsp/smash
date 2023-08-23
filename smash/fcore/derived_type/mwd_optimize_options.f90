!%      (MWD) Module Wrapped and Differentiated.
!%
!%      Type
!%      ----
!%
!%      - Optimize_OptionsDT
!%
!%          ================================== =======================================
!%          `Variables`                        Description
!%          ================================== =======================================
!%          ``mapping``                       Control mapping name
!%          ``optimizer``                     Optimizer name
!%          ``control_tfm``                   Type of transformation applied to control
!%          ``opr_parameters``                Opr parameters to optimize
!%          ``opr_initial_states``            Opr initial states to optimize
!%          ``l_opr_parameters``              Opr parameters lower bound
!%          ``u_opr_parameters``              Opr parameters upper bound
!%          ``l_opr_initial_states``          Opr initial states lower bound
!%          ``u_opr_initial_states``          Opr initial states upper bound
!%          ``opr_parameters_descriptor``     Opr parameters descriptor to use
!%          ``opr_initial_states_descriptor`` Opr initial states descriptor use
!%          ``maxiter``                       Maximum number of iterations
!%          ``factr``                         LBFGSB cost function criterion
!%          ``pgtol``                         LBFGSB gradient criterion
!%          ================================== =======================================
!%
!%      Subroutine
!%      ----------
!%
!%      - Optimize_OptionsDT_initialise
!%      - Optimize_OptionsDT_copy

module mwd_optimize_options

    use md_constant !% only: sp, lchar
    use mwd_setup !% only: SetupDT

    implicit none

    type Optimize_OptionsDT

        character(lchar) :: mapping = "..." !$F90W char
        character(lchar) :: optimizer = "..." !$F90W char
        character(lchar) :: control_tfm = "..." !$F90W char

        integer, dimension(:), allocatable :: opr_parameters
        real(sp), dimension(:), allocatable :: l_opr_parameters
        real(sp), dimension(:), allocatable :: u_opr_parameters
        integer, dimension(:, :), allocatable :: opr_parameters_descriptor

        integer, dimension(:), allocatable :: opr_initial_states
        real(sp), dimension(:), allocatable :: l_opr_initial_states
        real(sp), dimension(:), allocatable :: u_opr_initial_states
        integer, dimension(:, :), allocatable :: opr_initial_states_descriptor

        integer :: maxiter = -99
        real(sp) :: factr = -99._sp
        real(sp) :: pgtol = -99._sp

    end type Optimize_OptionsDT

contains

    subroutine Optimize_OptionsDT_initialise(this, setup)

        implicit none

        type(Optimize_OptionsDT), intent(inout) :: this
        type(SetupDT), intent(in) :: setup

        allocate (this%opr_parameters(setup%nop))
        this%opr_parameters = -99

        allocate (this%l_opr_parameters(setup%nop))
        this%l_opr_parameters = -99._sp

        allocate (this%u_opr_parameters(setup%nop))
        this%u_opr_parameters = -99._sp

        allocate (this%opr_parameters_descriptor(setup%nd, setup%nop))
        this%opr_parameters_descriptor = -99

        allocate (this%opr_initial_states(setup%nos))
        this%opr_initial_states = -99

        allocate (this%l_opr_initial_states(setup%nos))
        this%l_opr_initial_states = -99._sp

        allocate (this%u_opr_initial_states(setup%nos))
        this%u_opr_initial_states = -99._sp

        allocate (this%opr_initial_states_descriptor(setup%nd, setup%nos))
        this%opr_initial_states_descriptor = -99

    end subroutine Optimize_OptionsDT_initialise

    subroutine Optimize_OptionsDT_copy(this, this_copy)

        implicit none

        type(Optimize_OptionsDT), intent(in) :: this
        type(Optimize_OptionsDT), intent(out) :: this_copy

        this_copy = this

    end subroutine Optimize_OptionsDT_copy

end module mwd_optimize_options