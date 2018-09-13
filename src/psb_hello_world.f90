program psb_hello_world
  use psb_base_mod
  implicit none

  integer :: icontxt, np, ip

  call psb_init(icontxt)

  call psb_info(icontxt,ip,np)

  write(*,'(a,i0,a,i0,a)') "hello world from processor ",ip, ' in ',np, ' processor !'


  call psb_exit(icontxt)


end program psb_hello_world
