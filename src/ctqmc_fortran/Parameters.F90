!===============================================================================
module MParameters
!===============================================================================
   use iso_c_binding, only: c_int64_t, c_long_double

   integer,parameter                   :: KINDR=kind(0.d0),KINDC=kind(0.d0)
   integer,parameter                   :: PMAX=100
   integer,parameter                   :: stdin=5,stdout=6,stderr=0
   integer,parameter                   :: IntegerNotFound=-9999999
   integer,parameter                   :: idlength=100,vallength=500,linelength=5000
   real(KINDR),parameter               :: RealNotFound=-9999999d0
   real(KINDR)                         :: EPSLANC=1d-15
   real(KINDR)                         :: EPSDEG=1d-12


   type,private :: TParameter
      character(idlength)              :: ID=' '
      character(vallength)             :: Value=' '
      type(TParameter),pointer         :: next
   end type TParameter

   type(TParameter),pointer,private    :: first=>null()
   type(TParameter),pointer,private    :: last=>null()

   ! NWormSectors has to be defined in CTQMC.F90
   !integer, parameter :: NWormSectors = 15
   
   integer, parameter :: SectorZ = 1
   integer, parameter :: SectorG = 2
   integer, parameter :: SectorGSigma = 3
   integer, parameter :: SectorG4 = 4
   integer, parameter :: SectorH4 = 5
   integer, parameter :: SectorP2 = 6
   integer, parameter :: SectorP2pp = 7
   integer, parameter :: SectorP3 = 8
   integer, parameter :: SectorP3pp = 9
   integer, parameter :: SectorQQ = 10
   integer, parameter :: SectorQ4 = 11
   integer, parameter :: SectorNQQdag = 12
   integer, parameter :: SectorQQdd = 13
   integer, parameter :: SectorUcaca = 14
   integer, parameter :: SectorUccaa = 15
   integer, parameter :: SectorQUDdag = 16
   !number of worm operators in specific sector
   integer, parameter, dimension(2:16)   :: NOperWorm = (/2,4,4,6,4,4,4,4,6,12,8,8,4,4,4/)
contains
!===============================================================================
logical elemental function eq(real1,real2)
!===============================================================================
!input
   real(KINDR),intent(in)              :: real1,real2
   eq=abs(real1-real2).lt.1d-12
end function eq

!===============================================================================
subroutine set_Parameter(ID,val)
!===============================================================================
!input
   character(*),intent(in)             :: ID
   character(*),intent(in)             :: val
!local
   type(TParameter),pointer            :: tmp

   tmp=>first
   do while(index(tmp%ID,ID).ne.1)
      if(associated(tmp%next).eqv..false.)then
         write(stderr,*)ID," not found in parameter list!"
         stop
      endif
      tmp=>tmp%next
   enddo
   tmp%value(:)=" "
   tmp%value=val
end subroutine set_Parameter

!===============================================================================
integer function get_Integer_Parameter(ID)
!===============================================================================
!input
   character(*),intent(in)             :: ID
!local
   type(TParameter),pointer            :: tmp

   real :: rval
   tmp=>first
   do while(index(tmp%ID,ID).ne.1)
      if(associated(tmp%next).eqv..false.)then
         write(stderr,*)ID," not found in parameter list!"
         get_Integer_Parameter=IntegerNotFound
         stop
      endif
      tmp=>tmp%next
   enddo
   read(tmp%value,*,err=98)get_Integer_Parameter
   return

   ! Attempt with real to support scientific notation
98 read(tmp%value,*,err=99) rval
   if( (rval - int(rval)) > 2*epsilon(1.) ) goto 99
   get_Integer_Parameter = int(rval)
   return

99 write(stderr,*) ID," is not an integer"
   stop
end function get_Integer_Parameter

!===============================================================================
integer(c_int64_t) function get_LongInteger_Parameter(ID)
!===============================================================================
!input
   character(*),intent(in)             :: ID
!local
   type(TParameter),pointer            :: tmp
   real(kind=c_long_double)            :: rval

   tmp=>first
   do while(index(tmp%ID,ID).ne.1)
      if(associated(tmp%next).eqv..false.)then
         write(stderr,*)ID," not found in parameter list!"
         get_LongInteger_Parameter=IntegerNotFound
         stop
      endif
      tmp=>tmp%next
   enddo
   read(tmp%value,*,err=99)get_LongInteger_Parameter
   return

   ! Attempt with real to support scientific notation
99 read(tmp%value,*) rval
   if( (rval - int(rval, KIND=c_int64_t)) > 2*epsilon(1.) ) then
      write(stderr,*) ID," is not an integer"
      stop
   endif
   get_LongInteger_Parameter = int(rval,KIND=c_int64_t)

end function get_LongInteger_Parameter

!===============================================================================
real(KINDR) function get_Real_Parameter(ID)
!===============================================================================
!input
   character(*),intent(in)             :: ID
!local
   type(TParameter),pointer            :: tmp

   tmp=>first
   do while(index(tmp%ID,ID).ne.1)
      if(associated(tmp%next).eqv..false.)then
         write(stderr,*)ID," not found in parameter list!"
         get_Real_Parameter=RealNotFound
         stop
      endif
      tmp=>tmp%next
   enddo
   read(tmp%value, *, err=99) get_Real_Parameter
   return

99 write (*,*) 'Invalid parameter:', ID, '=', tmp%value
   stop 'get_Real_Parameter'
end function get_Real_Parameter


!===============================================================================
subroutine get_Real_List(ID,Val)
!===============================================================================
!input
   character(*),intent(in)          :: ID
!output
   real(KINDR),intent(out)          :: Val(:)
!local
   type(TParameter),pointer         :: tmp

   tmp=>first
   do while(index(tmp%ID,ID).ne.1)
      if(associated(tmp%next).eqv..false.)then
         write(stderr,*)ID," not found in parameter list!"
         stop
      endif
      tmp=>tmp%next
   enddo
   read(tmp%value,*)Val
end subroutine get_Real_List

!===============================================================================
subroutine get_Integer_List(ID,Val)
!===============================================================================
!input
   character(*),intent(in)          :: ID
!output
   integer,intent(out)              :: Val(:)
!local
   type(TParameter),pointer         :: tmp

   tmp=>first
   do while(index(tmp%ID,ID).ne.1)
      if(associated(tmp%next).eqv..false.)then
         write(stderr,*)ID," not found in parameter list!"
         stop
      endif
      tmp=>tmp%next
   enddo
   read(tmp%value,*)Val
end subroutine get_Integer_List

!===============================================================================
character(vallength) function get_String_Parameter(ID)
!===============================================================================
!input
   character(*),intent(in)     :: ID
!local
   type(TParameter),pointer            :: tmp

   tmp=>first
   do while(index(tmp%ID,ID).ne.1)
      if(associated(tmp%next).eqv..false.)then
         write(stderr,*)ID," not found in parameter list!"
         get_String_Parameter="not found"
         stop
      endif
      tmp=>tmp%next
   enddo
   get_String_Parameter=tmp%value
end function get_String_Parameter

!===============================================================================
subroutine push_Parameter(ID,Value)
!===============================================================================
!input
   character(*),intent(in)     :: ID,Value

   if(associated(first).eqv..false.)then
      allocate(first)
      last=>first
      first%next=>null()
   else
      allocate(last%next)
      last=>last%next
      last%next=>null()
   endif
   last%ID(:)=' '
   last%ID=ID
   last%Value(:)=' '
   last%Value=Value
end subroutine push_Parameter

!===============================================================================
subroutine read_ParameterFile(Filename)
!===============================================================================
!input
   character(*),intent(in)             :: Filename
!local
   integer                             :: stat
   character(idlength)                 :: id
   character(vallength)                :: val
   character(linelength)               :: line

   line(:)=' '
   id(:)=' '
   val(:)=' '
   open(100,file=Filename,status='old')
   stat=0
   do while(stat.ge.0)
      read(100,'(A)',iostat=stat)line
      if(index(line,"#").ne.1)then
         do while(index(line,char(09)).ne.0)
            line(index(line,char(09)):index(line,char(09)))=" "
         enddo
         if(index(line,"#").gt.0)line(index(line,"#"):)=" "
         if(index(line,"=").gt.0)then
            id=trim(adjustl(line(:index(line,"=")-1)))
            val=trim(adjustl(line(index(line,"=")+1:)))
            call push_Parameter(id,val)
         endif
      endif
   end do
   EPSLANC=get_Real_Parameter("EPSLANC")
   close(100)
end subroutine read_ParameterFile

!===============================================================================
subroutine read_ParameterString(ParameterString)
!===============================================================================
!input
   character(*),intent(in)             :: ParameterString
!local
   integer                             :: PosOld,PosNew
   character(idlength)                 :: id
   character(vallength)                :: val
   character(linelength)               :: line

   first=>null()
   last=>null()
   line(:)=' '
   id(:)=' '
   val(:)=' '
   PosOld=1
   PosNew=index(ParameterString(PosOld:),char(10))
   line=ParameterString(PosOld:PosNew)
   do while(PosNew.ne.0)
      line=ParameterString(PosOld:PosOld+PosNew-1)
      if(index(line,"#").ne.1)then
         do while(index(line,char(09)).ne.0)
            line(index(line,char(09)):index(line,char(09)))=" "
         enddo
         if(index(line,"#").gt.0)line(index(line,"#"):)=" "
         if(index(line,"=").gt.0)then
            id=trim(adjustl(line(:index(line,"=")-1)))
            val=trim(adjustl(line(index(line,"=")+1:)))
            do while(index(val,",").ne.0)
               val(index(val,","):index(val,","))=" "
            enddo
            call push_Parameter(id,val)
         endif
      endif
      PosOld=PosOld+PosNew
      PosNew=index(ParameterString(PosOld:),char(10))
   end do
   EPSLANC=get_Real_Parameter("EPSLANC")
end subroutine read_ParameterString

!===============================================================================
subroutine dest_Parameters()
!===============================================================================
!local
   type(TParameter),pointer            :: tmp1,tmp2

   tmp1=>first
   do while(associated(tmp1))
      tmp2=>tmp1
      tmp1=>tmp1%next
      deallocate(tmp2)
   enddo
   first=>null()
end subroutine dest_Parameters
!===============================================================================
end module MParameters
!===============================================================================

#ifdef Parameters_Test

!===============================================================================
program Prog_Parameters
!===============================================================================
use MParameters
character(1000)   :: ParaStr

write(stdout,*)"================================================================"
write(stdout,*)"            reading parameter file Parameters.in"
write(stdout,*)"================================================================"
call read_ParameterFile("Parameters.in")
write(stdout,*)"            writing the converted parameter values"
write(stdout,*)"================================================================"
write(stdout,*)"U ",get_Real_Parameter("U")
write(stdout,*)"J ",get_Real_Parameter("J")
write(stdout,*)"Hamiltonian: ",get_String_Parameter("Hamiltonian")
write(stdout,*)"NBands: ",get_Integer_Parameter("Nd")
write(stdout,*)"beta: ",get_Real_Parameter("beta")
write(stdout,*)"NCorr: ",get_Integer_Parameter("NCorr")
write(stdout,*)"================================================================"
write(stdout,*)"            setting parameter to a different value"
write(stdout,*)"================================================================"
call set_Parameter("U","1.073846")
write(stdout,*)"U ",get_Real_Parameter("U")
write(stdout,*)"================================================================"
call dest_Parameters()

end program Prog_Parameters

#endif
