program ifanimal
implicit none
logical :: isacow,isadog,hastwohorns
integer, parameter :: nhorns = 2
isacow = .true.
isadog = .false.
if (isacow) then ! check if our animal is a cow
write(*,''(a)'',advance='no') my animal is a cow...'
if (nhorns==2) write(*,*) ' ...with two horns'
elseif (isadog) then
print*,' my animal is a dog. Woof.'
else
print ``(a)'','my animal is not a cow'
hastwohorns = (nhorns==2)
if (hastwohorns) print ``(a)'',' but it has two horns'
endif
end program ifanimal
