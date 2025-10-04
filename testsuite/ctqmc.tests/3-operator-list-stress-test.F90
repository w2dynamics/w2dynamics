program test
   use iso_c_binding
   use MTrace
   use MRandom
   implicit none

   interface
      subroutine abort() bind(c)
      end subroutine abort
   end interface

   type(TTrace)         :: trace
   type(TOperPointer)   :: oper(2)
   type(TOper), pointer :: element
   integer              :: i, j, pos(2)
   integer, allocatable :: seed(:), prewinfcmp(:, :, :)
   integer, parameter   :: steps=75000, nbands=3, nslide=50
   real, parameter      :: pinsert=0.501, beta=100.0, tdmax=10.0
   real(kind(0.0d0))    :: throwaway
   logical              :: firstc, lastc, prewindow, inwindow
   logical, allocatable :: ffirstc(:, :), flastc(:, :)

   allocate(trace%ffirst(nbands, 2))
   allocate(trace%flast(nbands, 2))
   allocate(trace%fprewin(nbands, 2))
   allocate(trace%NOSOper(nbands, 2))
   trace%NOSOper = 0
   allocate(trace%NOSCAOper(nbands, 2, 2))
   trace%NOSCAOper = 0
   allocate(trace%PreWinFCount(nbands, 2, 2))
   trace%beta = beta
   trace%b_offdiag = .true.

   trace%window_position = 0
   trace%taudiff_max = tdmax
   trace%num_windows = max(floor(trace%beta / trace%taudiff_max) - 1, 1)
   trace%tauwin_min = (dble(trace%window_position)/(dble(trace%num_windows + 1))) * trace%beta
   trace%tauwin_max = (dble(trace%window_position + 2)/(dble(trace%num_windows + 1))) * trace%beta

   trace%NPairs = 0
   trace%PreWinFCount = 0
   trace%PreWinFCountValid = .true.

   trace%prewin => null()
   trace%postwin => null()
   trace%last_prewin => null()
   do i = 1, size(trace%fprewin, 1)
      do j = 1, size(trace%fprewin, 2)
         trace%fprewin(i, j)%p => null()
      end do
   end do

   call set_window(trace)

   call random_seed()
   call random_seed(size=i)
   allocate(seed(i))
   call random_seed(get=seed)
   do j = 2, i
      seed(1) = ieor(seed(1), seed(j))
   end do
   call sgrnd(seed(1))
   deallocate(seed)

   do i = 1, steps
      if (mod(i, nslide) == 0) call shift_window(trace)

      if (grnd() <= pinsert) then
         ! add a pair
         do j = 1, 2
            allocate(oper(j)%p)
            oper(j)%p%has_hyb=.true.
         end do

         call pair_OperAdd(trace, nbands, oper, 2, throwaway, force_diagonal=.false.)

         if (.not. process_OperAdd(trace, oper, pos, 2)) then
            do j = 1, 2
               deallocate(oper(j)%p)
            end do
            cycle
         end if

         call insert_Oper(trace, oper)

         call update_NPairs(trace, oper, 1)
      else
         ! remove a pair
         if (.not. gen_OperRemove(trace, oper, pos, 2, throwaway, with_hyb=.true.)) then
            cycle
         end if

         call remove_Oper(trace, oper)

         call update_NPairs(trace, oper, -1)

         do j = 1, 2
            deallocate(oper(j)%p)
         end do
      end if
   end do

   firstc = .false.
   lastc = .false.
   allocate(ffirstc(nbands, 2))
   ffirstc = .false.
   allocate(flastc(nbands, 2))
   flastc = .false.
   allocate(prewinfcmp(nbands, 2, 2))
   prewinfcmp = 0

   do i = 1, size(trace%fprewin, 1)
      do j = 1, size(trace%fprewin, 2)
         if (associated(trace%fprewin(i, j)%p)) then
            if (trace%fprewin(i, j)%p%tau >= trace%tauwin_min) call test_error(__LINE__)
            if (associated(trace%fprewin(i, j)%p%fnext)) then
               if (trace%fprewin(i, j)%p%fnext%tau < trace%tauwin_min) call test_error(__LINE__)
            end if
         end if
      end do
   end do
   if (associated(trace%prewin)) then
      if (trace%prewin%tau >= trace%tauwin_min) call test_error(__LINE__)
      if (associated(trace%prewin%next)) then
         if (trace%prewin%next%tau < trace%tauwin_min) call test_error(__LINE__)
      end if
   end if
   if (associated(trace%postwin)) then
      if (trace%postwin%tau <= trace%tauwin_max) call test_error(__LINE__)
      if (associated(trace%postwin%prev)) then
         if (trace%postwin%prev%tau > trace%tauwin_max) call test_error(__LINE__)
      end if
   end if

   if (associated(trace%prewin)) then
      prewindow = .true.
      inwindow = .false.
   else
      prewindow = .false.
      inwindow = .true.
   end if

   element => trace%first
   if (trace%noper > 0) then
      if (associated(element%prev)) call test_error(__LINE__)
   end if
   do i = 1, trace%noper
      if (associated(trace%postwin, element)) then
         if (prewindow) call test_error(__LINE__)
         if (inwindow) then
            inwindow = .false.
         else
            call test_error(__LINE__)
         end if
      end if

      if (prewindow) then
         if (element%tau >= trace%tauwin_min) call test_error(__LINE__)
         prewinfcmp(element%orbital, element%spin, element%ca)&
            = prewinfcmp(element%orbital, element%spin, element%ca) + 1
      else if (inwindow) then
         if (element%tau < trace%tauwin_min&
             .or. element%tau > trace%tauwin_max) call test_error(__LINE__)
      else
         if (element%tau <= trace%tauwin_max) call test_error(__LINE__)
      end if

      if (associated(trace%prewin, element)) then
         if (prewindow .and. .not. inwindow) then
            prewindow = .false.
            inwindow = .true.
         else
            call test_error(__LINE__)
         end if
      end if

      if (i /= 1) then
         if (.not. associated(element%prev)) call test_error(__LINE__)
      end if

      if (associated(element%prev)) then
         if (element%prev%tau >= element%tau) call test_error(__LINE__)
         if (.not. associated(element%prev%next, element)) call test_error(__LINE__)
      else
         if (.not. associated(trace%first, element)) call test_error(__LINE__)
         if (.not. firstc) then
            firstc = .true.
         else
            call test_error(__LINE__)
         end if
      end if

      if (associated(element%fprev)) then
         if (element%fprev%tau >= element%tau) call test_error(__LINE__)
         if (element%fprev%orbital /= element%orbital&
             .or. element%fprev%spin /= element%spin) call test_error(__LINE__)
         if (.not. associated(element%fprev%fnext, element)) call test_error(__LINE__)
      else
         if (.not. associated(trace%ffirst(element%orbital,&
                                           element%spin)%p,&
                              element)) call test_error(__LINE__)
         if (.not. ffirstc(element%orbital, element%spin)) then
            ffirstc(element%orbital, element%spin) = .true.
         else
            call test_error(__LINE__)
         end if
      end if

      if (associated(element%next)) then
         if (element%next%tau <= element%tau) call test_error(__LINE__)
         if (.not. associated(element%next%prev, element)) call test_error(__LINE__)
      else
         if (.not. associated(trace%last, element)) call test_error(__LINE__)
         if (.not. lastc) then
            lastc = .true.
         else
            call test_error(__LINE__)
         end if
      end if

      if (associated(element%fnext)) then
         if (element%fnext%tau <= element%tau) call test_error(__LINE__)
         if (element%fnext%orbital /= element%orbital&
             .or. element%fnext%spin /= element%spin) call test_error(__LINE__)
         if (.not. associated(element%fnext%fprev, element)) call test_error(__LINE__)
      else
         if (.not. associated(trace%flast(element%orbital,&
                                          element%spin)%p,&
                              element)) call test_error(__LINE__)
         if (.not. flastc(element%orbital, element%spin)) then
            flastc(element%orbital, element%spin) = .true.
         else
            call test_error(__LINE__)
         end if
      end if

      if (i == trace%noper) then
         if (.not. associated(trace%last, element)) call test_error(__LINE__)
         if (associated(element%next)) call test_error(__LINE__)
      else
         if (.not. associated(element%next)) call test_error(__LINE__)
      end if
      element => element%next
   end do

   if (.not. firstc) then
      if (associated(trace%first)) call test_error(__LINE__)
   end if
   if (.not. lastc) then
      if (associated(trace%last)) call test_error(__LINE__)
   end if

   call ensure_valid_PreWinFCount(trace)
   if (any(trace%PreWinFCount /= prewinfcmp)) call test_error(__LINE__)

   do i = 1, size(ffirstc, 1)
      do j = 1, size(ffirstc, 2)
         if (.not. ffirstc(i, j)) then
            if (associated(trace%ffirst(i, j)%p)) call test_error(__LINE__)
         end if
         if (.not. flastc(i, j)) then
            if (associated(trace%flast(i, j)%p)) call test_error(__LINE__)
         end if
      end do
   end do

   deallocate(ffirstc, flastc, prewinfcmp)

contains

   subroutine test_error(line)
      integer :: line

      write (*, "('Operator list stress test error in line ', I4)") line
      call abort()
   end subroutine test_error

end program test
