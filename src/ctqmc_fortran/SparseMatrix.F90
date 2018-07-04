!implementation of a csr-matrix and a matrix-vector multiplication
!matrix is created dynamically by subroutine init_SparseMatrix

!===============================================================================
module MSparseMatrix
!===============================================================================
use MParameters

   type :: TSparseMatrix
      integer                          :: dim1
      integer,allocatable              :: ptr(:),indices(:)
      real(KINDR),allocatable          :: dat(:)
   end type TSparseMatrix

contains
!===============================================================================
subroutine init_SparseMatrix(this,Mat,dim1,dim2)
!===============================================================================
   type(TSparseMatrix)                 :: this
!input
   integer                             :: dim1,dim2,i,j,NElem
   real(KINDR)                         :: Mat(dim1,dim2)
   
   this%dim1=dim1
   NElem=0
   do i=1,dim1
      do j=1,dim2
         if(.not.eq(Mat(i,j),0d0))NElem=NElem+1
      enddo
   enddo

   allocate(this%ptr(dim1+1)) 
   allocate(this%indices(NElem))
   allocate(this%dat(NElem))

   this%ptr(dim1+1)=NElem+1
   NElem=1
   do i=1,dim1
      this%ptr(i)=NElem
      do j=1,dim2
         if(.not.eq(Mat(i,j),0d0))then
            this%indices(NElem)=j
            this%dat(NElem)=Mat(i,j)
            NElem=NElem+1
         endif
      enddo
   enddo
end subroutine init_SparseMatrix

!===============================================================================
function MatVecProd(this,x)
!===============================================================================
   type(TSparseMatrix)                 :: this
!input
   real(KINDR),intent(in)              :: x(:)
!output
   real(KINDR)                         :: MatVecProd(this%dim1)
!local
   integer                             :: i,j

   MatVecProd=0d0
   do i=1,this%dim1
      do j=this%ptr(i),this%ptr(i+1)-1
         MatVecProd(i)=MatVecProd(i)+this%dat(j)*x(this%indices(j))
      enddo      
   enddo
end function MatVecProd

!===============================================================================
subroutine SMatVecProd(this,x,y)
!===============================================================================
   type(TSparseMatrix)                 :: this
!input
   real(KINDR),intent(in)              :: x(:)
!output
   real(KINDR)                         :: y(:)
!local
   integer                             :: i,j

   y=0d0
   do i=1,this%dim1
      do j=this%ptr(i),this%ptr(i+1)-1
         y(i)=y(i)+this%dat(j)*x(this%indices(j))
      enddo      
   enddo
end subroutine SMatVecProd

!===============================================================================
subroutine SMatVecProdAdd(this,x,y)
!===============================================================================
   type(TSparseMatrix)                 :: this
!input
   real(KINDR),intent(in)              :: x(:)
!output
   real(KINDR)                         :: y(:)
!local
   integer                             :: i,j

   do i=1,this%dim1
      do j=this%ptr(i),this%ptr(i+1)-1
         y(i)=y(i)+this%dat(j)*x(this%indices(j))
      enddo      
   enddo
end subroutine SMatVecProdAdd

!===============================================================================
subroutine print_SparseMatrix(this)
!===============================================================================
   type(TSparseMatrix)                 :: this
!local
   integer                             :: i,j

   do i=1,size(this%ptr)-1
      do j=this%ptr(i),this%ptr(i+1)-1
         write(stdout,*)"Row/Col/Data",i,"/",this%indices(j),"/",this%dat(j)
      enddo
   enddo
end subroutine print_SparseMatrix

!===============================================================================
subroutine dest_SparseMatrix(this)
!===============================================================================
   type(TSparseMatrix)                 :: this 
   if(allocated(this%ptr))&
      deallocate(this%ptr)
   if(allocated(this%indices))&
      deallocate(this%indices)
   if(allocated(this%dat))&
      deallocate(this%dat)
end subroutine dest_SparseMatrix

!===============================================================================
end module MSparseMatrix
!===============================================================================

#ifdef SparseMatrix_Test

!===============================================================================
program Prog_SparseMatrix
!===============================================================================
use MParameters
use MSparseMatrix
real(KINDR),dimension(2,4)                :: Mat=reshape((/1d0,0d0,3d0,0d0,5d0,0d0,4d0,2d0/),(/2,4/))
real(KINDR),dimension(4)                  :: x=(/1d0,0d0,1d0,1d0/)
real(KINDR),dimension(2)                  :: y
type(TSparseMatrix)                       :: SparseMat


write(stdout,*)"================================================================"
write(stdout,*)"           creating the intial sparse matrix"
call init_SparseMatrix(SparseMat,Mat,2,4)
write(stdout,*)"================================================================"
write(stdout,*)"           printing the matrix in sparse form"
write(stdout,*)"================================================================"
call print_SparseMatrix(SparseMat)
write(stdout,*)"================================================================"
write(stdout,*)"           calculating y=Ax"
y=MatVecProd(SparseMat,x)
write(stdout,*)"================================================================"
write(stdout,*)"           the result is:"
write(stdout,*)y
write(stdout,*)"================================================================"
write(stdout,*)"           calculating y=Ax with matmul"
write(stdout,*)"================================================================"
write(stdout,*)Matmul(Mat,x)
write(stdout,*)"================================================================"
write(stdout,*)"           destructing the sparse matrix"
call dest_SparseMatrix(SparseMat)
write(stdout,*)"================================================================"

end program Prog_SparseMatrix

#endif
