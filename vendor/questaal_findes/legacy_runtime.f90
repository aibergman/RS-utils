! Small in-memory replacements for the Questaal application services used by
! FINDES. The search is invoked in mode=1, so file-writing and command-line
! services are intentionally inert; numerical geometry and array operations
! remain active.
subroutine info0(a,b,c,message)
  implicit none
  integer, intent(in) :: a,b,c
  character(*), intent(in) :: message
end subroutine
subroutine info2(a,b,c,message)
  implicit none
  integer, intent(in) :: a,b,c
  character(*), intent(in) :: message
end subroutine
subroutine info5(a,b,c,message)
  implicit none
  integer, intent(in) :: a,b,c
  character(*), intent(in) :: message
end subroutine
integer function nglob(name)
  implicit none
  character(*), intent(in) :: name
  nglob=6
end function
integer function isw(value)
  implicit none
  logical, intent(in) :: value
  isw=merge(1,0,value)
end function
logical function cmdopt(option,n,mode,buffer)
  implicit none
  character(*), intent(in) :: option
  integer, intent(in) :: n,mode
  character(*), intent(inout) :: buffer
  cmdopt=.false.
end function
logical function a2bin(text,value,n1,n2,separator,index,mode)
  implicit none
  character(*), intent(in) :: text,separator
  logical, intent(out) :: value
  integer, intent(in) :: n1,n2,index,mode
  value=.false.
  a2bin=.false.
end function
subroutine numsyv(index)
  implicit none
  integer, intent(out) :: index
  index=0
end subroutine
subroutine lodsyv(name,index,value,status)
  implicit none
  character(*), intent(in) :: name
  integer, intent(in) :: index
  real(8), intent(in) :: value
  integer, intent(out) :: status
  status=0
end subroutine
subroutine clrsyv(index)
  implicit none
  integer, intent(in) :: index
end subroutine
subroutine pshpr(value)
  implicit none
  integer, intent(in) :: value
end subroutine
subroutine poppr()
end subroutine
subroutine rx(message)
  implicit none
  character(*), intent(in) :: message
  error stop 'Questaal FINDES runtime error'
end subroutine
subroutine rx0(message)
  implicit none
  character(*), intent(in) :: message
end subroutine
subroutine shorps(n,plat,mode,source,destination)
  implicit none
  integer, intent(in) :: n,mode(3)
  real(8), intent(in) :: plat(3,3),source(3,* )
  real(8), intent(out) :: destination(3,* )
  integer :: i
  do i=1,n
    destination(:,i)=source(:,i)
  end do
end subroutine
subroutine ptr_ctrl(s_ctrl,mode,element,n1,n2,n3,array)
  use structures, only: str_ctrl
  implicit none
  type(str_ctrl), intent(inout) :: s_ctrl
  integer, intent(in) :: mode,n1,n2,n3,array(*)
  character(*), intent(in) :: element
  integer, pointer :: work(:)
  integer :: n
  n=max(0,n1)
  if (trim(element)=='ipc') then
    allocate(work(n)); work=0
    if (associated(s_ctrl%ipc)) work(:min(n,size(s_ctrl%ipc)))=s_ctrl%ipc(:min(n,size(s_ctrl%ipc)))
    if (iand(mode,4)/=0) work(:n)=array(:n)
    if (associated(s_ctrl%ipc)) deallocate(s_ctrl%ipc)
    s_ctrl%ipc=>work
  elseif (trim(element)=='ips') then
    allocate(work(n)); work=0
    if (associated(s_ctrl%ips)) work(:min(n,size(s_ctrl%ips)))=s_ctrl%ips(:min(n,size(s_ctrl%ips)))
    if (iand(mode,4)/=0) work(:n)=array(:n)
    if (associated(s_ctrl%ips)) deallocate(s_ctrl%ips)
    s_ctrl%ips=>work
  end if
end subroutine
subroutine ptr_lat(s_lat,mode,element,n1,n2,n3,array)
  use structures, only: str_lat
  implicit none
  type(str_lat), intent(inout) :: s_lat
  integer, intent(in) :: mode,n1,n2,n3
  character(*), intent(in) :: element
  real(8), intent(in) :: array(*)
  real(8), pointer :: work(:,:)
  integer :: nold
  if (trim(element)/='pos') return
  nold=0
  if (associated(s_lat%pos)) nold=size(s_lat%pos,2)
  allocate(work(n1,n2)); work=0.0d0
  if (nold>0) work(:,1:min(nold,n2))=s_lat%pos(:,1:min(nold,n2))
  if (iand(mode,4)/=0) work(:n1,:n2)=reshape(array(1:n1*n2),[n1,n2])
  if (associated(s_lat%pos)) deallocate(s_lat%pos)
  s_lat%pos=>work
end subroutine
subroutine daxpy(n,alpha,x,incx,y,incy)
  implicit none
  integer, intent(in) :: n,incx,incy
  real(8), intent(in) :: alpha,x(*)
  real(8), intent(inout) :: y(*)
  integer :: i
  do i=0,n-1
    y(1+i*incy)=y(1+i*incy)+alpha*x(1+i*incx)
  end do
end subroutine
subroutine dpcopy(x,y,n1,n2,scale)
  implicit none
  integer, intent(in) :: n1,n2
  real(8), intent(in) :: x(*),scale
  real(8), intent(out) :: y(*)
  integer :: i
  do i=n1,n2
    y(i)=scale*x(i)
  end do
end subroutine
subroutine dpzero(x,n)
  implicit none
  integer, intent(in) :: n
  real(8), intent(out) :: x(*)
  x(1:n)=0.0d0
end subroutine
subroutine icopy(n,x,incx,y,incy)
  implicit none
  integer, intent(in) :: n,incx,incy,x(*)
  integer, intent(out) :: y(*)
  integer :: i
  do i=0,n-1
    y(1+i*incy)=x(1+i*incx)
  end do
end subroutine
subroutine iinit(x,n)
  implicit none
  integer, intent(in) :: n
  integer, intent(out) :: x(*)
  x(1:n)=0
end subroutine
double precision function dasum(n,x,incx)
  implicit none
  integer, intent(in) :: n,incx
  real(8), intent(in) :: x(*)
  integer :: i
  dasum=0.0d0
  do i=0,n-1
    dasum=dasum+abs(x(1+i*incx))
  end do
end function
double precision function dnrm2(n,x,incx)
  implicit none
  integer, intent(in) :: n,incx
  real(8), intent(in) :: x(*)
  integer :: i
  dnrm2=sqrt(sum(x(1:1+(n-1)*incx:incx)**2))
end function
subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  implicit none
  character, intent(in) :: transa,transb
  integer, intent(in) :: m,n,k,lda,ldb,ldc
  real(8), intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
  real(8), intent(inout) :: c(ldc,*)
  integer :: i,j,l
  real(8) :: av,bv
  do j=1,n
    do i=1,m
      c(i,j)=beta*c(i,j)
      do l=1,k
        if (transa=='N' .or. transa=='n') then
          av=a(i,l)
        else
          av=a(l,i)
        end if
        if (transb=='N' .or. transb=='n') then
          bv=b(l,j)
        else
          bv=b(j,l)
        end if
        c(i,j)=c(i,j)+alpha*av*bv
      end do
    end do
  end do
end subroutine
subroutine dmpy(a,nca,nra,b,ncb,nrb,c,ncc,nrc,n,m,l)
  implicit none
  integer, intent(in) :: nca,nra,ncb,nrb,ncc,nrc,n,m,l
  real(8), intent(in) :: a(*),b(*)
  real(8), intent(out) :: c(*)
  integer :: i,j,k
  do i=0,n-1
    do j=0,m-1
      c(1+i*nrc+j*ncc)=0.0d0
      do k=0,l-1
        c(1+i*nrc+j*ncc)=c(1+i*nrc+j*ncc)+a(1+i*nra+k*nca)*b(1+k*nrb+j*ncb)
      end do
    end do
  end do
end subroutine
subroutine nrmliz(n,x,y)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: x(3,n)
  real(8), intent(out) :: y(3,n)
  real(8) :: norm
  integer :: i
  do i=1,n
    norm=sqrt(sum(x(:,i)*x(:,i)))
    if (norm>0.0d0) then
      y(:,i)=x(:,i)/norm
    else
      y(:,i)=0.0d0
    end if
  end do
end subroutine
double precision function d1mach(index)
  implicit none
  integer, intent(in) :: index
  select case(index)
  case(1); d1mach=tiny(1.0d0)
  case(2); d1mach=huge(1.0d0)
  case(3); d1mach=epsilon(1.0d0)/2.0d0
  case(4); d1mach=epsilon(1.0d0)
  case default; d1mach=log10(huge(1.0d0))
  end select
end function
double precision function danrm2(x)
  implicit none
  real(8), intent(in) :: x(3)
  danrm2=sum(x*x)
end function
subroutine dvheap(n,m,a,iperm,tol,mode)
  implicit none
  integer, intent(in) :: n,m,mode
  real(8), intent(inout) :: a(n,m)
  integer, intent(out) :: iperm(m)
  real(8), intent(in) :: tol
  integer :: i,j
  do i=1,m
    iperm(i)=i
  end do
  do i=2,m
    j=i
    do while(j>1)
      if (a(1,iperm(j))>=a(1,iperm(j-1))-tol) exit
      block
        integer :: t
        t=iperm(j); iperm(j)=iperm(j-1); iperm(j-1)=t
      end block
      j=j-1
    end do
  end do
end subroutine
subroutine d1mach_dummy()
end subroutine
subroutine fopn_dummy()
end subroutine
integer function fopn(name)
  implicit none
  character(*), intent(in) :: name
  fopn=99
end function
integer function iosite(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v)
  implicit none
  integer, intent(in) :: a,b,c
  character(*), intent(in) :: d,e
  integer, intent(in) :: f,g
  character(*), intent(in) :: h
  real(8), intent(in) :: i,j,k,l,m
  integer, intent(in) :: n
  real(8), intent(in) :: o,p,q,r,s,t,u,v
  iosite=0
end function
subroutine awrit2(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
  implicit none
  character(*), intent(in) :: a,b,c,d
  integer, intent(in) :: e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z
end subroutine
subroutine chcase(direction,n,string)
  implicit none
  integer, intent(in) :: direction,n
  character(*), intent(inout) :: string
end subroutine
integer function chrpos(string,n,char)
  implicit none
  integer, intent(in) :: n
  character(*), intent(in) :: string,char
  chrpos=index(string(:min(n,len(string))),trim(char))
end function
subroutine uniqclabl(iopt,nclass,mode,classes,nlist,labels,label)
  implicit none
  integer, intent(in) :: iopt,nclass,mode,nlist,classes(nlist)
  character(8), intent(inout) :: labels(nclass),label
  integer :: suffix,i
  character(8) :: base,candidate
  base=label; suffix=1
  do
    candidate=base
    if (suffix>1) write(candidate(2:8),'(i7)') suffix
    if (.not.any(labels(:max(0,nclass-1))==candidate)) exit
    suffix=suffix+1
  end do
  label=candidate
  labels(classes(1))=candidate
end subroutine
subroutine zslabl(iopt,label,z)
  implicit none
  integer, intent(in) :: iopt
  character(*), intent(inout) :: label
  integer, intent(inout) :: z
  if (iopt==-1) then
    if (z==0) label='E'
  end if
end subroutine
integer function a2vec(text,values)
  implicit none
  character(*), intent(in) :: text
  integer, intent(out) :: values(*)
  a2vec=0
end function
integer function idamx_stub(n,x,inc)
  implicit none
  integer, intent(in) :: n,inc
  real(8), intent(in) :: x(*)
  idamx_stub=1
end function
