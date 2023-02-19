!
! module handle input to and from files
!

module io


implicit none

public :: writeoutfile, read_file

private

contains

!write data file
subroutine writeoutfile(filename,x,y,nx,nt)
    character(len=20), intent(in) :: filename
    double precision, intent(in) :: x(nx),y(nx,nt)  ! arrays will automatically take size of input array
    integer, intent(in) :: nt,nx
    integer :: i,iu,it
    character(len=20) :: fmtstring ="(es12.6,1x,es12.6)"


    print*,filename

    print "(a,i5)", 'writing to'//trim(filename)//'on unit ',iu
    open(unit=iu,file=filename,status='replace',action='write')
    write(iu,*) nx,nt
    write(iu,*) x
    do it=1,nt
        !do i=1,min(size(x),size(y))
            write(iu,*) y(:,it)
        !enddo
    enddo
    close(unit=iu)

end subroutine writeoutfile

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

!write a model atmosphere file
subroutine writeatfile(height, dens, press, temp,bruntvas, nitems)
    implicit none
    integer, intent(in) :: nitems
    real, intent(in) :: height(nitems), dens(nitems), press(nitems), temp(nitems),bruntvas(nitems)
    integer :: i,iu,ierr


    open(unit=iu, file='atmos.txt')

    !starting at point number 4 because lower points had -ve enrgy density
    do i=37,nitems
        write(iu,*) height(i),temp(i),dens(i), press(i), bruntvas(i)
200     format(F16.6,2X,F16.6,2X,F16.6,2X,F16.6)
    end do

    close(unit=9)

end subroutine writeatfile


!write a model atmosphere file
subroutine readatfile(height, dens, press, temp, nitems)
    implicit none
    integer, intent(in) :: nitems
    real, intent(inout) :: height(nitems), dens(nitems), press(nitems), temp(nitems)
    integer :: i,iu,ierr


    open(unit=iu,file='atmos.txt',status='old',action='read')

    i=1
    do while(ierr==0)
        read(iu,*,iostat=ierr) height(i),temp(i),dens(i), press(i)
        i=i+1
    enddo
    close(unit=iu)



end subroutine readatfile






end module io
