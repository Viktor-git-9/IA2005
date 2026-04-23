subroutine print_mem(label)
  character(len=*), intent(in) :: label
  integer :: unit, ios
  character(len=256) :: line
  open(newunit=unit, file='/proc/self/status', status='old', iostat=ios)
  if (ios /= 0) return
  do
    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    if (line(1:6) == 'VmRSS:' .or. line(1:5) == 'VmPeak' &
        .or. line(1:5) == 'VmSize') then
      write(*,'(A,A,A)') '[MEM] ', trim(label), ': '//trim(line)
    end if
  end do
  close(unit)
end subroutine print_mem