!---- never change ------

#include "fintrf.h"      

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

   implicit none
   
interface
 subroutine gcd(vp,gshort,c,d,nn,ne,nf,nv,f)
    integer*8 vp(4,*)
    integer*8 nn,ne,nf,nv
    integer*8,allocatable,dimension(:,:) :: gshort,c,d,f
  end subroutine  
end interface

! mwPointer mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
!---- end never change ------
	  
	  
! mwSize stuff for mexing
      mwSize mo,no,siz
	  
! mwPointer stuff for mexing
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetPr
	  mwPointer vp_pr,g_pr,c_pr,d_pr,nn_pr,ne_pr,nf_pr,nv_pr, f_pr
      mwPointer m, n
      mwPointer mxGetM, mxGetN

!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
	  
! fortran subroutine arguments
	  real*8,allocatable,dimension(:,:) :: vp_r
	  real*8,allocatable,dimension(:,:) :: g_dbl,c_dbl,d_dbl,f_dbl
      integer*8,allocatable,dimension(:,:) :: vp,g,c,d,f
	  integer*8 nn, ne, nf, nv
!	  real*8 
	  character*80 msg
      logical debu
       
      debu = .false. ! .true. o .false. per attivare o disattivare il debug
	  if(debu) open(unit=66,file='log.txt',status='unknown')

!    Check to see input is numeric.
!	  do ii = 1,6
        if (mxIsNumeric(prhs(1)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:assemble_L_1d:NonNumeric', &
                                'Inputs must be numeric.')
        endif
!	  enddo
	  if(debu) write(66,*) 'num_check'
	  
	  
!     Check that input #1 is scalar matrix and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
     ! if(m .ne. 4 .or. n .ne. N_hexa*3) then
      if(m .ne. 4) then 
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must 4 x NV')
      endif	  
	  nv=n
      siz = m*n
      vp_pr = mxGetPr(prhs(1))
	  allocate(vp_r(4,n))
      allocate(vp(4,n))
      call mxCopyPtrToReal8(vp_pr, vp_r, siz) !  double precision to real
	  vp = int(vp_r) ! real to integer
	  if(debu) write(66,*) 'input 1'

	  	  

! call the computational subroutine.

       call gcd(vp,g,c,d,nn,ne,nf,nv,f)
       if(debu) write(66,*) 'called'

	  deallocate(vp,vp_r)
! Create a matrix for the return argument 1 (g)
      mo=2
      no=ne
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output 1 into a MATLAB array.
      g_pr = mxGetPr(plhs(1))
      siz=mo*no
	  allocate(g_dbl(mo,no))
	  g_dbl=real(g,8)
      call mxCopyReal8ToPtr(g_dbl, g_pr, siz)
      if(debu) write(66,*) 'G'
	  deallocate(g)
	  
! Create a matrix for the return argument 2 (c)
      mo=3
      no=nf
	  ComplexFlag = 0
	  if(debu) write(66,*) 'mxCreateDoubleMatrix'
      plhs(2) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output2 into a MATLAB array.
      c_pr = mxGetPr(plhs(2))
      siz=mo*no
	  if(debu) write(66,*) 'mxCopyReal8ToPtr, sizes C and nf:',size(c,1),size(c,2),nf
	  allocate(c_dbl(mo,no))
	  c_dbl=real(c,8)	  
      call mxCopyReal8ToPtr(c_dbl, c_pr, siz)
      if(debu) write(66,*) 'C'
	  deallocate(c)
	  
! Create a matrix for the return argument 3 (d)
      mo=4
      no=nv
	  ComplexFlag = 0
      plhs(3) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output 3 into a MATLAB array.
      d_pr = mxGetPr(plhs(3))
      siz=mo*no
	  allocate(d_dbl(mo,no))
	  d_dbl=real(d,8)		  
      call mxCopyReal8ToPtr(d_dbl, d_pr, siz)
      if(debu) write(66,*) 'D'
	  deallocate(d)
	  
! Create a matrix for the return argument 4 (f)
      mo=5
      no=nf
	  ComplexFlag = 0
	  if(debu) write(66,*) 'mxCreateDoubleMatrix'
      plhs(4) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output2 into a MATLAB array.
      f_pr = mxGetPr(plhs(4))
      siz=mo*no
	   if(debu) write(66,*) 'mxCopyReal8ToPtr, sizes F and nf:',size(f,1),size(f,2),nf
	  allocate(f_dbl(mo,no))
	  f_dbl=real(f,8)	   
      call mxCopyReal8ToPtr(f_dbl, f_pr, siz)
      if(debu) write(66,*) 'F'	  
      deallocate(f)

      if(debu) write(66,*) 'matrices to matlab'
      
      if(debu) write(66,*) 'closing log.txt'
	  
      
      if(debu) close(66)
      return
      end