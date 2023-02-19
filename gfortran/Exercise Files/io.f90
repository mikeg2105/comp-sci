!
! module handle input to and from files
!

module io


implicit none

public :: write_file, read_file

private

contains

!write data file
subroutine write_file(filename,x,y)
    character(len=20), intent(in) :: filename
    real, intent(in) :: x(:),y(:)  ! arrays will automatically take size of input array
    integer :: i,iu
    character(len=20) :: fmtstring ="(es12.6,1x,es12.6)"


    print*,filename

    print "(a,i5)", 'writing to'//trim(filename)//'on unit ',iu
    open(unit=iu,file=filename,status='replace',action='write')
    do i=1,min(size(x),size(y))
        write(iu,fmtstring) x(i),y(i)
    enddo
    close(unit=iu)

end subroutine write_file

!read data file
subroutine read_file(filename, x, y)
    character(len=20), intent(in) :: filename
    real, intent(out) :: x(:),y(:)  ! arrays will automatically take size of input array
    integer :: i,iu,ierr

    print "(a,i5)", 'reading from '//trim(filename)//'on unit ',iu

    open(unit=iu,file=filename,status='old',action='read')
    !do i=1,min(size(x),size(y))
    do while(ierr==0)
        read(iu,*,iostat=ierr) x(i),y(i)
    enddo
    close(unit=iu)



end subroutine read_file




end module io
