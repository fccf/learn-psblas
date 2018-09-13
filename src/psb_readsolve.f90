program psb_read_solve
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  implicit none

  integer :: ictxt, np, iam, unit

  ! dense matrices
  real(psb_dpk_), allocatable, target ::  aux_b(:,:), d(:)
  real(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  real(psb_dpk_), pointer  :: b_col_glob(:)
  type(psb_d_vect_type)    :: b_col, x_col, r_col

  ! sparse matrices
  type(psb_dspmat_type) :: aux_a,a

  ! communications data structure
  type(psb_desc_type):: desc_a

  ! preconditioner data
  type(psb_dprec_type)  :: prec

  ! solver paramters
  integer(psb_ipk_) :: iter, itmax, ierr, itrace, ircode,&
       & methd, istopc, irst
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  real(psb_dpk_)   :: err, eps, cond

  ! other variables
  integer(psb_ipk_) :: info
  integer(psb_ipk_) :: m_problem, i
  integer(psb_ipk_), allocatable :: perm(:)
  character(:),allocatable :: mtrx_file, afmt, ptype, kmethd

  real(psb_dpk_) :: t1, t2, tprec
  real(psb_dpk_) :: r_amax, b_amax, scale,resmx,resmxp
  integer(psb_ipk_) :: iparm(20)

  mtrx_file = "pde2961.mtx"
  afmt = "CSR"
  ptype = "BJAC"
  kmethd = "BiCGSTAB"

  itmax  = 500
  itrace = -1
  istopc = 2
  irst   = 2
  eps    = 1.0e-7

  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  call psb_bcast(ictxt,mtrx_file)
  call psb_bcast(ictxt,kmethd)
  call psb_bcast(ictxt,ptype)
  call psb_bcast(ictxt,afmt)
  call psb_bcast(ictxt,itmax)
  call psb_bcast(ictxt,itrace)
  call psb_bcast(ictxt,istopc)
  call psb_bcast(ictxt,irst)
  call psb_bcast(ictxt,eps)

  call psb_barrier(ictxt)
  t1 = psb_wtime()

  if (iam == psb_root_) then
    !> read matrix
    call mm_mat_read(aux_a,info,filename=mtrx_file)
    if (info /= psb_success_) error stop "Unable to read matrix file :"//mtrx_file


    !> set RHS
    m_problem = aux_a%get_nrows()
    call psb_bcast(ictxt,m_problem)
    call psb_mat_renum(psb_mat_renum_identity_,aux_a,info,perm)

    call psb_realloc(m_problem,1,aux_b,ierr)
    if(ierr /=0) error stop "Unable to allocate RHS"

    !
    b_col_glob => aux_b(:,1)
    do i=1, m_problem
      b_col_glob(i) = done
    enddo
  else
    call psb_bcast(ictxt,m_problem) !<?
  endif


  !> Partition
  if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
  call psb_matdist(aux_a, a,  ictxt,desc_a,info,fmt=afmt,parts=part_block)

  !>
  call psb_scatter(b_col_glob,b_col,desc_a,info,root=psb_root_)
  call psb_geall(x_col,desc_a,info)
  call x_col%zero()
  call psb_geasb(x_col,desc_a,info)
  call psb_geall(r_col,desc_a,info)
  call r_col%zero()
  call psb_geasb(r_col,desc_a,info)

  t2 = psb_wtime() - t1


  call psb_amx(ictxt, t2)

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
    write(psb_out_unit,'(" ")')
  end if

  !> building the preconditioner
  call prec%init(ptype,info)
  t1 = psb_wtime()

  call prec%build(a,desc_a,info)

  tprec = psb_wtime()-t1
  if (info /= psb_success_) error stop "Unable to building the preconditioner"

  call psb_amx(ictxt,tprec)

  if(iam == psb_root_) then
    write(psb_out_unit,'("Preconditioner time: ",es12.5)')tprec
    write(psb_out_unit,'(" ")')
  end if

  !> krylov
  cond = dzero
  iparm = 0
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_krylov(kmethd,a,prec,b_col,x_col,eps,desc_a,info,&
  & itmax=itmax,iter=iter,err=err,itrace=itrace,&
  & istop=istopc,irst=irst,cond=cond)

  !> output info
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  call psb_amx(ictxt,t2)
  call psb_geaxpby(done,b_col,dzero,r_col,desc_a,info)
  call psb_spmm(-done,a,x_col,done,r_col,desc_a,info)
  resmx  = psb_genrm2(r_col,desc_a,info)
  resmxp = psb_geamax(r_col,desc_a,info)

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)
  if (iam == psb_root_) then
    call prec%descr()
    write(psb_out_unit,'("Matrix: ",a)')mtrx_file
    write(psb_out_unit,'("Computed solution on ",i8," processors")')np
    write(psb_out_unit,'("Iterations to convergence: ",i6)')iter
    write(psb_out_unit,'("Error estimate on exit   : ",es12.5)') err
    write(psb_out_unit,'("Time to buil prec.       : ",es12.5)')tprec
    write(psb_out_unit,'("Time to solve system     : ",es12.5)')t2
    write(psb_out_unit,'("Time per iteration       : ",es12.5)')t2/(iter)
    write(psb_out_unit,'("Total time               : ",es12.5)')t2+tprec
    write(psb_out_unit,'("Residual norm 2          : ",es12.5)')resmx
    write(psb_out_unit,'("Residual norm inf        : ",es12.5)')resmxp
    write(psb_out_unit,'("Condition number         : ",es12.5)')cond
    write(psb_out_unit,'("Total memory occupation for A:      ",i12)')amatsize
    write(psb_out_unit,'("Total memory occupation for PREC:   ",i12)')precsize
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize
    write(psb_out_unit,'("Storage type for DESC_A           : ",a)')&
    &  desc_a%get_fmt()
  end if

  !> ouput results
  call psb_gather(x_col_glob,x_col,desc_a,info,root=psb_root_)
  if (info == psb_success_) &
       & call psb_gather(r_col_glob,r_col,desc_a,info,root=psb_root_)
  if (info /= psb_success_) error stop "Unable to gather data"

  if (iam == psb_root_) then
    open(newunit=unit,file="results")

    write(psb_err_unit,'(" ")')
    write(psb_err_unit,'("Saving x on file")')
    write(unit,*) 'matrix: ',mtrx_file
    write(unit,*) 'computed solution on ',np,' processors.'
    write(unit,*) 'iterations to convergence: ',iter
    write(unit,*) 'error estimate (infinity norm) on exit:', &
         & ' ||r||/(||a||||x||+||b||) = ',err
    write(unit,'("Residual norm 2          : ",es12.5)')resmx
    write(unit,'("Residual norm inf        : ",es12.5)')resmxp
    write(unit,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
    do i=1,m_problem
      write(unit,998) i,x_col_glob(i),r_col_glob(i),b_col_glob(i)
    enddo

    close(unit)
  end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))

  !> free
  call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call prec%free(info)
  call psb_cdfree(desc_a,info)

  call psb_exit(ictxt)

end program psb_read_solve
