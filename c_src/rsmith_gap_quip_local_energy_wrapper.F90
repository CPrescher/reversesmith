subroutine rsmith_quip_local_energy_wrapper(nlocal, nghost, atomic_numbers, lmptag, &
   inum, sum_num_neigh, ilist, quip_num_neigh, quip_neigh, lattice, &
   quip_potential, n_quip_potential, quip_x, quip_local_e) bind(c)

   use iso_c_binding, only: c_double, c_int
   use system_module
   use linearalgebra_module
   use connection_module
   use atoms_types_module
   use atoms_module
   use potential_module
   implicit none

   integer(kind=c_int), intent(in) :: nlocal, nghost, inum, sum_num_neigh
   integer(kind=c_int), intent(in), dimension(nlocal+nghost) :: atomic_numbers
   integer(kind=c_int), intent(in), dimension(nlocal+nghost) :: lmptag
   integer(kind=c_int), intent(in), dimension(inum) :: ilist
   integer(kind=c_int), intent(in), dimension(inum) :: quip_num_neigh
   integer(kind=c_int), intent(in), dimension(sum_num_neigh) :: quip_neigh
   real(kind=c_double), intent(in), dimension(3,3) :: lattice
   integer(kind=c_int), intent(in), dimension(n_quip_potential) :: quip_potential
   integer(kind=c_int), intent(in) :: n_quip_potential
   real(kind=c_double), intent(in), dimension(3,nlocal+nghost) :: quip_x
   real(kind=c_double), intent(out), dimension(nlocal+nghost) :: quip_local_e

   type local_quip_potential
      type(Potential), pointer :: pot
   end type local_quip_potential

   integer :: i, j, n, nn, ni, i_n1n, j_n2n, error
   integer, save :: nn_guess = 0
   type(atoms), save :: at
   type(local_quip_potential) :: transferred
   type(Potential), pointer :: pot
   logical, dimension(:), pointer :: local
   real(dp) :: r_ij
   real(dp), parameter :: small_number = 1.0e-10_dp

   if (n_quip_potential == 0) then
      call system_abort('rsmith_quip_local_energy_wrapper: potential not initialised')
   else
      transferred = transfer(quip_potential, transferred)
      pot => transferred%pot
   endif

   if (nlocal <= 0) then
      quip_local_e = 0.0_dp
      return
   endif

   if (.not. is_initialised(at) .or. at%N /= (nlocal+nghost) .or. &
       nn_guess < maxval(quip_num_neigh)) then
      call finalise(at)
      call initialise(at, nlocal+nghost, lattice)
      call set_cutoff(at, cutoff(pot))
      nn_guess = 2 * maxval(quip_num_neigh)
      call connection_fill(at%connect, at%N, at%N, nn_guess=nn_guess)
   endif

   at%Z = atomic_numbers
   at%pos = quip_x
   call add_property(at, 'local', .false., ptr=local, overwrite=.true., error=error)
   call add_property(at, 'lmptag', lmptag, overwrite=.true., error=error)

   do i = 1, at%N
      at%connect%neighbour1(i)%t%N = 0
      at%connect%neighbour2(i)%t%N = 0
   enddo

   nn = 0
   do ni = 1, inum
      i = ilist(ni) + 1
      local(i) = .true.
      i_n1n = 0
      do n = 1, quip_num_neigh(ni)
         nn = nn + 1
         j = quip_neigh(nn)
         if (i <= j) then
            i_n1n = i_n1n + 1
            at%connect%neighbour1(i)%t%N = at%connect%neighbour1(i)%t%N + 1
            at%connect%neighbour1(i)%t%int(1,i_n1n) = j
            at%connect%neighbour1(i)%t%int(2:4,i_n1n) = 0
            r_ij = norm(at%pos(:,i) - at%pos(:,j))
            if (r_ij < small_number) then
               call system_abort('rsmith_quip_local_energy_wrapper: atoms overlap')
            endif
            at%connect%neighbour1(i)%t%real(1,i_n1n) = r_ij
            at%connect%neighbour2(j)%t%N = at%connect%neighbour2(j)%t%N + 1
            j_n2n = at%connect%neighbour2(j)%t%N
            at%connect%neighbour2(j)%t%int(1,j_n2n) = i
            at%connect%neighbour2(j)%t%int(2,j_n2n) = i_n1n
         endif
      enddo
   enddo

   call calc(pot, at, local_energy=quip_local_e, &
      args_str='atom_mask_name=local pers_idces_name=lmptag lammps do_calc_connect=F')

end subroutine rsmith_quip_local_energy_wrapper
