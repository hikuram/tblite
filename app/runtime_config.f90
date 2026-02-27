! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Runtime configuration loader for tblite CLI.
module tblite_runtime_config
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number
   use tblite_toml, only : toml_table, toml_array, toml_error, toml_parse, get_value, len
   use tblite_solvation, only : solvation_input, alpb_input, cds_input, shift_input, &
      & solution_state, born_kernel
   implicit none
   private

   public :: load_run_defaults_from_config

contains

subroutine load_run_defaults_from_config(filename, method, param, solvation, error)
   character(len=*), intent(in) :: filename
   character(len=:), allocatable, intent(inout) :: method
   character(len=:), allocatable, intent(inout) :: param
   type(solvation_input), allocatable, intent(inout) :: solvation
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), allocatable :: root
   type(toml_error), allocatable :: terr
   type(toml_table), pointer :: t_method, t_solv, t_floor
   character(len=:), allocatable :: s
   character(len=:), allocatable :: model, solvent, kernel_s, state_s, mode_s
   real(wp) :: eps, kappa_default
   logical :: one_sided
   integer :: kernel, sol_state, floor_mode
   integer, parameter :: max_z = 118
   real(wp), allocatable :: rmin(:), kappa(:)
   type(toml_array), pointer :: a_elem, a_rmin, a_kappa
   integer :: i, n
   character(len=:), allocatable :: sym

   open(file=filename, newunit=i, status='old', action='read')
   call toml_parse(root, i, terr)
   close(i)
   if (allocated(terr)) then
      call fatal_error(error, "Failed to parse config file '"//trim(filename)//"': "//trim(terr%message))
      return
   end if

   ! [method]
   call get_value(root, "method", t_method, null())
   if (associated(t_method)) then
      call get_value(t_method, "method", method)
      call get_value(t_method, "param", param)
   end if

   ! [solvation]
   call get_value(root, "solvation", t_solv, null())
   if (.not.associated(t_solv)) return

   call get_value(t_solv, "model", model)
   call get_value(t_solv, "solvent", solvent)
   call get_value(t_solv, "kernel", kernel_s)
   call get_value(t_solv, "state", state_s)

   ! Dielectric constant: allow explicit eps, otherwise try to use solvent name only.
   eps = -1.0_wp
   call get_value(t_solv, "eps", eps, -1.0_wp)

   kernel = born_kernel%p16
   if (allocated(kernel_s)) then
      select case(trim(kernel_s))
      case("p16")
         kernel = born_kernel%p16
      case("still")
         kernel = born_kernel%still
      case default
         call fatal_error(error, "Unknown born kernel '"//trim(kernel_s)//"' in config")
         return
      end select
   end if

   sol_state = solution_state%gsolv
   if (allocated(state_s)) then
      select case(trim(state_s))
      case("gsolv")
         sol_state = solution_state%gsolv
      case("bar1mol")
         sol_state = solution_state%bar1mol
      case("reference")
         sol_state = solution_state%reference
      case default
         call fatal_error(error, "Unknown solution state '"//trim(state_s)//"' in config")
         return
      end select
   end if

   ! [solvation.born_floor]
   floor_mode = 0
   one_sided = .true.
   kappa_default = 5.0_wp
   allocate(rmin(max_z), kappa(max_z))
   rmin = -1.0_wp
   kappa = -1.0_wp

   call get_value(t_solv, "born_floor", t_floor, null())
   if (associated(t_floor)) then
      call get_value(t_floor, "mode", mode_s)
      if (allocated(mode_s)) then
         select case(trim(mode_s))
         case("off")
            floor_mode = 0
         case("soft")
            floor_mode = 1
         case("hard")
            floor_mode = 2
         case default
            call fatal_error(error, "Unknown born_floor.mode '"//trim(mode_s)//"' in config")
            return
         end select
      end if

      call get_value(t_floor, "one_sided", one_sided, .true.)
      call get_value(t_floor, "kappa_default", kappa_default, 5.0_wp)

      call get_value(t_floor, "elements", a_elem, null())
      call get_value(t_floor, "rmin_bohr", a_rmin, null())
      call get_value(t_floor, "kappa", a_kappa, null())

      if (associated(a_elem) .or. associated(a_rmin) .or. associated(a_kappa)) then
         if (.not.(associated(a_elem) .and. associated(a_rmin))) then
            call fatal_error(error, "born_floor.elements and born_floor.rmin_bohr must be provided together")
            return
         end if
         if (len(a_elem) /= len(a_rmin)) then
            call fatal_error(error, "born_floor.elements and born_floor.rmin_bohr must have same length")
            return
         end if
         if (associated(a_kappa)) then
            if (len(a_kappa) /= len(a_elem)) then
               call fatal_error(error, "born_floor.kappa must have same length as elements if provided")
               return
            end if
         end if

         n = len(a_elem)
         do i = 1, n
            call get_value(a_elem, i, sym)
            if (.not.allocated(sym)) cycle
            rmin(to_number(sym)) = get_real_from_array(a_rmin, i, error)
            if (allocated(error)) return
            if (associated(a_kappa)) then
               kappa(to_number(sym)) = get_real_from_array(a_kappa, i, error)
               if (allocated(error)) return
            end if
         end do
      end if
   end if

   ! Allocate solvation input (ALPB/GBSA only for now)
   if (.not.allocated(model)) then
      call fatal_error(error, "solvation.model must be set in config")
      return
   end if

   if (.not.allocated(solvation)) allocate(solvation)

   select case(trim(model))
   case("alpb", "gbsa")
      if (.not.allocated(solvent)) then
         call fatal_error(error, "solvation.solvent must be set for alpb/gbsa in config")
         return
      end if
      ! Note: eps is required by constructor; for named solvent, caller can supply eps explicitly.
      if (eps <= 0.0_wp) then
         call fatal_error(error, "solvation.eps must be provided in config when using solvation.model='"//trim(model)//"'")
         return
      end if

      solvation%alpb = alpb_input(eps, solvent=solvent, kernel=kernel, alpb=(trim(model)=="alpb"))
      solvation%cds = cds_input(alpb=(trim(model)=="alpb"), solvent=solvent)
      solvation%shift = shift_input(alpb=(trim(model)=="alpb"), solvent=solvent, state=sol_state)

      solvation%alpb%floor_mode = floor_mode
      solvation%alpb%floor_one_sided = one_sided
      solvation%alpb%floor_kappa_default = kappa_default
      allocate(solvation%alpb%floor_rmin(max_z), solvation%alpb%floor_kappa(max_z))
      solvation%alpb%floor_rmin = rmin
      solvation%alpb%floor_kappa = kappa

   case default
      call fatal_error(error, "Unsupported solvation.model '"//trim(model)//"' in config (only alpb/gbsa supported)")
      return
   end select

end subroutine load_run_defaults_from_config

function get_real_from_array(arr, idx, error) result(val)
   type(toml_array), intent(inout) :: arr
   integer, intent(in) :: idx
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: val
   real(wp) :: tmp

   tmp = 0.0_wp
   call get_value(arr, idx, tmp)
   val = tmp
end function get_real_from_array

end module tblite_runtime_config
