!===============================================================================
!! This module is for debugging; it prints vectors and 2D Matrices in a formatted way.
!! The syntax is: call ausgabe($Object,shape($Object),"$Object")
!! $Object can be a matrix or vector of
!! real, integer or KINDROC
!! and allocatable or not allocatable.
!! Written by Andreas Hausoel
module MAusgabe
!===============================================================================

use MParameters

interface ausgabe
   module procedure print_matrix_real
   module procedure print_matrix_int
   module procedure print_matrix_KINDR
   module procedure print_matrix_KINDR_file
   module procedure print_vector_real
   module procedure print_vector_int
   module procedure print_vector_KINDR
   module procedure print_vector_KINDR_file
end interface

contains

!===============================================================================
subroutine print_matrix_real(matrix,sh,string,stellen,bounds)
!===============================================================================
integer, intent(in), dimension(2) :: sh
integer :: iz41, is41, stellen
real, dimension(sh(1),sh(2)), intent(in) :: matrix
character(len=*) :: string 
character(len=60) :: format_string
character(len=15) :: s2
character(len=15) :: s1
logical :: bounds

write(s1,*) 7+stellen
write(s2,*) stellen

s1=trim(s1)
s2=trim(s2)

format_string= trim("(ES"//s1//"."//s2//",A2)")
! write(unit=*,fmt=*) format_string


write(*,"(A)",advance="no") trim(string)
write(*,"(A14)",advance="no") ",  Dimension:"
write(unit=*,fmt="(I2,'  x ',I2)") shape(matrix)
if (bounds) write(unit=*,fmt="(A9)",advance="no") "Bounds:  "
if (bounds) write(unit=*,fmt="(A3,I3,A4,I3,A8,I3,A4,I3)") "x  ", lbound(matrix,1),&
& "--> ", ubound(matrix,1), "     y  ", lbound(matrix,2), "--> ",ubound(matrix,2)

do iz41=lbound(matrix,1),ubound(matrix,1)
do is41=lbound(matrix,2),ubound(matrix,2)
write(*,format_string,advance="no") Matrix(iz41,is41), " "
!"(F6.3,A2)"
enddo
write(*,*)
enddo
write(*,*)
end subroutine print_matrix_real

!===============================================================================
subroutine print_matrix_int(matrix,sh,string,stellen,bounds)
!===============================================================================
integer, intent(in), dimension(2) :: sh
integer :: iz41, is41, stellen
integer, dimension(sh(1),sh(2)), intent(in) :: matrix
character(len=*) :: string 
character(len=60) :: format_string
character(len=15) :: s1
logical :: bounds

write(s1,*) stellen

s1=trim(s1)

format_string= trim("(I"//s1//",A2)")
! write(unit=*,fmt=*) format_string


write(*,"(A)",advance="no") trim(string)
write(*,"(A14)",advance="no") ",  Dimension:"
write(unit=*,fmt="(I2,'  x ',I2)") shape(matrix)
if (bounds) write(unit=*,fmt="(A9)",advance="no") "Bounds:  "
if (bounds) write(unit=*,fmt="(A3,I3,A4,I3,A8,I3,A4,I3)") "x  ", lbound(matrix,1),&
& "--> ", ubound(matrix,1), "     y  ", lbound(matrix,2), "--> ",ubound(matrix,2)

do iz41=lbound(matrix,1),ubound(matrix,1)
do is41=lbound(matrix,2),ubound(matrix,2)
write(*,format_string,advance="no") Matrix(iz41,is41), "  "
enddo
write(*,*)
enddo
write(*,*)
end subroutine print_matrix_int

!===============================================================================
subroutine print_matrix_KINDR(matrix,sh,string,stellen,bounds)
!===============================================================================
integer, intent(in), dimension(2) :: sh
integer :: iz41, is41, stellen
real(KINDR), dimension(sh(1),sh(2)), intent(in) :: matrix
character(len=*) :: string 
character(len=60) :: format_string
character(len=15) :: s2
character(len=15) :: s1
logical :: bounds

write(s1,*) 7+stellen
write(s2,*) stellen

s1=trim(s1)
s2=trim(s2)

format_string= trim("(ES"//s1//"."//s2//",A2)")
! write(unit=*,fmt=*) format_string


write(*,"(A)",advance="no") trim(string)
write(*,"(A14)",advance="no") ",  Dimension:"
write(unit=*,fmt="(I2,'  x ',I2,A5)",advance="yes") shape(matrix), "     "
if (bounds) write(unit=*,fmt="(A9)",advance="no") "Bounds:  "
if (bounds) write(unit=*,fmt="(A3,I3,A4,I3,A8,I3,A4,I3)") "x  ", lbound(matrix,1),&
& "--> ", ubound(matrix,1), "     y  ", lbound(matrix,2), "--> ",ubound(matrix,2)

do iz41=lbound(matrix,1),ubound(matrix,1)
do is41=lbound(matrix,2),ubound(matrix,2)
write(*,format_string,advance="no") Matrix(iz41,is41), " "
enddo
write(*,*)
enddo
write(*,*)
end subroutine print_matrix_KINDR



!===============================================================================
subroutine print_matrix_KINDR_file(matrix,sh,string,stellen,bounds,filenumber)
!===============================================================================
integer, intent(in), dimension(2) :: sh
integer :: iz41, is41, stellen, filenumber
real(KINDR), dimension(sh(1),sh(2)), intent(in) :: matrix
character(len=*) :: string 
character(len=60) :: format_string
character(len=15) :: s2
character(len=15) :: s1
logical :: bounds

write(s1,*) 7+stellen
write(s2,*) stellen

s1=trim(s1)
s2=trim(s2)

format_string= trim("(ES"//s1//"."//s2//",A2)")
! write(unit=*,fmt=*) format_string


write(filenumber,"(A)",advance="no") trim(string)
write(filenumber,"(A14)",advance="no") ",  Dimension:"
write(unit=filenumber,fmt="(I2,'  x ',I2,A5)",advance="yes") shape(matrix), "     "
if (bounds) write(unit=filenumber,fmt="(A9)",advance="no") "Bounds:  "
if (bounds) write(unit=filenumber,fmt="(A3,I3,A4,I3,A8,I3,A4,I3)") "x  ",&
& lbound(matrix,1), "--> ", ubound(matrix,1), "     y  ", lbound(matrix,2), "--> ",ubound(matrix,2)

do iz41=lbound(matrix,1),ubound(matrix,1)
do is41=lbound(matrix,2),ubound(matrix,2)
write(filenumber,format_string,advance="no") Matrix(iz41,is41), "     "
enddo
write(filenumber,*)
enddo
write(filenumber,*)
end subroutine print_matrix_KINDR_file

!===============================================================================
subroutine print_vector_real(vektor,sh,string,stellen,bounds)
!===============================================================================
integer, dimension(1), intent(in) :: sh
integer :: iii, stellen
real, dimension(sh(1)), intent(in) :: vektor
character(len=*) :: string 
character(len=60) :: format_string
character(len=15) :: s2
character(len=15) :: s1
logical :: bounds

write(s1,*) 7+stellen
write(s2,*) stellen

s1=trim(s1)
s2=trim(s2)

format_string= trim("(ES"//s1//"."//s2//",A2)")
! write(unit=*,fmt=*) format_string


write(*,"(A)",advance="no") trim(string)
write(*,"(A14)",advance="no") ",  Dimension:"
write(*,"(I3)") sh
if (bounds) write(unit=*,fmt="(A9)",advance="no") "Bound:  "
if (bounds) write(unit=*,fmt="(I3,A4,I3)") lbound(vektor), "--> ", ubound(vektor)

write(unit=*,fmt="(A1)", advance="no") "(" 
do iii=lbound(vektor,1),ubound(vektor,1)
write(*,format_string,advance="no") vektor(iii), " "
enddo
write(unit=*,fmt="(A1)") ")" 
write(*,*) ""
end subroutine print_vector_real

!===============================================================================
subroutine print_vector_int(vektor,sh,string,stellen,bounds)
!===============================================================================
integer, dimension(1), intent(in) :: sh
integer :: iii, stellen
integer, dimension(sh(1)), intent(in) :: vektor
character(len=*) :: string 
character(len=60) :: format_string
character(len=15) :: s1
logical :: bounds

write(s1,*) stellen

s1=trim(s1)

format_string= trim("(I"//s1//",A2)")
! write(unit=*,fmt=*) format_string


write(*,"(A)",advance="no") trim(string)
write(*,"(A14)",advance="no") ",  Dimension:"
write(*,"(I3)") sh
if (bounds) write(unit=*,fmt="(A9)",advance="no") "Bound:  "
if (bounds) write(unit=*,fmt="(I3,A4,I3)") lbound(vektor), "--> ", ubound(vektor)

write(unit=*,fmt="(A2)", advance="no") "( "  
do iii=lbound(vektor,1),ubound(vektor,1)
write(*,format_string,advance="no") vektor(iii), " "
enddo
write(unit=*,fmt="(A1)") ")" 
write(*,*) ""
end subroutine print_vector_int

!===============================================================================
subroutine print_vector_KINDR(vektor,sh,string,stellen,bounds)
!===============================================================================
integer, dimension(1), intent(in) :: sh
integer :: iii, stellen
real(KINDR), dimension(sh(1)), intent(in) :: vektor
character(len=*) :: string 
character(len=60) :: format_string
character(len=15) :: s2
character(len=15) :: s1
logical :: bounds

write(s1,*) 7+stellen
write(s2,*) stellen

s1=trim(s1)
s2=trim(s2)

format_string= trim("(ES"//s1//"."//s2//",A2)")
! write(unit=*,fmt=*) format_string


write(*,"(A)",advance="no") trim(string)
write(*,"(A14)",advance="no") ",  Dimension:"
write(*,"(I3)") sh
if (bounds) write(unit=*,fmt="(A9)",advance="no") "Bound:  "
if (bounds) write(unit=*,fmt="(I3,A4,I3)") lbound(vektor), "--> ", ubound(vektor)

write(unit=*,fmt="(A2)", advance="no") "( " 
do iii=lbound(vektor,1),ubound(vektor,1)
write(*,format_string,advance="no") vektor(iii), " "
enddo
write(unit=*,fmt="(A1)") ")" 
write(*,*) ""
end subroutine print_vector_KINDR


!===============================================================================
subroutine print_vector_KINDR_file(vektor,sh,string,stellen,bounds,filenumber)
!===============================================================================
integer, dimension(1), intent(in) :: sh
integer :: iii, stellen, filenumber
real(KINDR), dimension(sh(1)), intent(in) :: vektor
character(len=*) :: string 
character(len=60) :: format_string
character(len=15) :: s2
character(len=15) :: s1
logical :: bounds

write(s1,*) 7+stellen
write(s2,*) stellen

s1=trim(s1)
s2=trim(s2)

format_string= trim("(ES"//s1//"."//s2//",A2)")
! write(unit=*,fmt=*) format_string


write(filenumber,"(A)",advance="no") trim(string)
write(filenumber,"(A14)",advance="no") ",  Dimension:"
write(filenumber,"(I3)") sh
if (bounds) write(unit=*,fmt="(A9)",advance="no") "Bound:  "
if (bounds) write(unit=*,fmt="(I3,A4,I3)") lbound(vektor), "--> ", ubound(vektor)

write(unit=filenumber,fmt="(A2)", advance="no") "( " 
do iii=lbound(vektor,1),ubound(vektor,1)
write(filenumber,format_string,advance="no") vektor(iii), " "
enddo
write(unit=filenumber,fmt="(A1)") ")" 
write(filenumber,*) ""
end subroutine print_vector_KINDR_file


end module MAusgabe

#ifdef ausgabe_test

!===============================================================================
program Prog_Ausgabe
!===============================================================================

use MAusgabe

integer,parameter :: KINDR=kind(0.d0)
integer :: i,j
integer, parameter :: N=5
integer, parameter :: offset=3

! Test-Objects beginning with index 1
real, allocatable :: array1(:)
! real, dimension(N) :: array1
integer, dimension(N) :: array2
real(KINDR), dimension(N) :: array3
real, dimension(N,N) :: matrix1
integer, dimension(N,N) :: matrix2
real(KINDR), dimension(N,N) :: matrix3

! Test-Objects beginning with index 0
real, dimension(0:N-1) :: array11
integer, dimension(0:N-1) :: array22
real(KINDR), dimension(0:N-1) :: array33
real, dimension(0:N-1,0:N-1) :: matrix11
integer, dimension(0:N-1,0:N-1) :: matrix22
real(KINDR), dimension(0:N-1,0:N-1) :: matrix33

! Test-Objects beginning with index 0
real, dimension(offset:N+offset-1) :: array111
integer, dimension(offset:N+offset-1) :: array222
real(KINDR), dimension(offset:N+offset-1) :: array333
real, dimension(offset:N+offset-1,offset:N+offset-1) :: matrix111
integer, dimension(offset:N+offset-1,offset:N+offset-1) :: matrix222
real(KINDR), dimension(offset:N+offset-1,offset:N+offset-1) :: matrix333

allocate(array1(N))

!Set initial Values
do i = 1,N
array1(i)=i
array2(i)=i
array3(i)=i
enddo

do i = 1,N
do j = 1,N
matrix1(i,j)=i+j
matrix2(i,j)=i+j
matrix3(i,j)=i+j
enddo
enddo

!Call functions
call ausgabe(array1,shape(array1),"array1")
call ausgabe(array2,shape(array2),"array2")
call ausgabe(array3,shape(array3),"array3")

call ausgabe(matrix1,shape(matrix11),"matrix1")
call ausgabe(matrix2,shape(matrix22),"matrix2")
call ausgabe(matrix3,shape(matrix33),"matrix3")
!

!Set initial Values
do i = 0,N-1
array11(i)=i
array22(i)=i
array33(i)=i
enddo

do i = 0,N-1
do j = 0,N-1
matrix11(i,j)=i+j
matrix22(i,j)=i+j
matrix33(i,j)=i+j
enddo
enddo

!Call functions
call ausgabe(array11,shape(array11),"array11")
call ausgabe(array22,shape(array22),"array22")
call ausgabe(array33,shape(array33),"array33")

call ausgabe(matrix11,shape(matrix11),"matrix11")
call ausgabe(matrix22,shape(matrix22),"matrix22")
call ausgabe(matrix33,shape(matrix33),"matrix33")
! 
! 

!Set initial Values
do i = offset,N+offset-1
array111(i)=i
array222(i)=i
array333(i)=i
enddo

do i = offset,N+offset-1
do j = offset,N+offset-1
matrix111(i,j)=i+j
matrix222(i,j)=i+j
matrix333(i,j)=i+j
enddo
enddo

!Call functions
call ausgabe(array111,shape(array111),"array111")
call ausgabe(array222,shape(array222),"array222")
call ausgabe(array333,shape(array333),"array333")

call ausgabe(matrix111,shape(matrix111),"matrix111")
call ausgabe(matrix222,shape(matrix222),"matrix222")
call ausgabe(matrix333,shape(matrix333),"matrix333")

end program Prog_Ausgabe

#endif
