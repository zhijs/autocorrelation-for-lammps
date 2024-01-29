program polymer_analysis
  use omp_lib

  implicit none

  character(len=256) :: infile, outfile
  integer :: num_threads, t, num_molecules, i, j, k
  integer total
  integer :: id
  real :: x, y, z
  integer :: ierr, istat
  real, allocatable :: coords(:,:)
  integer, allocatable :: ids(:)
  integer, dimension(:,:), allocatable :: neighbor_list

  ! Read input arguments
  if (command_argument_count() < 1) then
      print *, "Usage: polymer_analysis <polymer_COM.txt>"
      stop
  endif
  call get_command_argument(1, infile)
  outfile = "neighbors.txt"

  ! Initialize OpenMP
  num_threads = omp_get_max_threads()
  call omp_set_num_threads(num_threads)

  ! Open output file for writing
  open(unit=10, file=outfile, status='replace')

  ! Open input file for reading
  open(unit=20, file=infile, status='old', action='read')
  read(20, *, iostat=istat)
  read(20, *, iostat=istat)
  read(20, *, iostat=istat)
  ! Loop over all timesteps
  do
      ! Read timestep
      read(20, *, iostat=istat) t, num_molecules
      if (istat /= 0) exit

      ! Allocate memory for coords and ids
      allocate(coords(num_molecules,3), ids(num_molecules))
      
      ! Read coords and ids
      do i = 1, num_molecules
          read(20, *, iostat=istat) ids(i), coords(i,1), coords(i,2), coords(i,3)
      end do

      ! Allocate memory for neighbor_list
      allocate(neighbor_list(num_molecules,num_molecules))

      ! Initialize neighbor_list
      neighbor_list = 0
      total=0
     
      !$OMP PARALLEL DO PRIVATE(i,j) reduction(+:total) SHARED(num_molecules,coords,neighbor_list)
      
      do i = 1, num_molecules-1
          do j = i+1, num_molecules
              if (i /= j) then
                  if (distance(coords(i,:), coords(j,:)) < 4.0) then
                      total=total+1
                      neighbor_list(i,j) = 1
                      neighbor_list(j,i) = 1
                  endif
              endif
          end do
      end do
      !$OMP END PARALLEL DO

      ! Write neighbor_list to file
      if (total > 0) then
          write(10, "(I6)") total
      endif
      do i = 1, num_molecules-1
          do j = i+1, num_molecules
              if (neighbor_list(i,j) == 1) then
                  write(10, "(I8, 2X, I6, 2X, I6)")t, ids(i), ids(j)                
              endif
          end do
      end do
      ! Free memory
      deallocate(coords, ids, neighbor_list)
  end do

  ! Close files
  close(10)
  close(20)

contains

  function distance(a, b) result(dist)
    real, dimension(:), intent(in) :: a, b
    real :: dist
    dist = sqrt(sum((a-b)**2))
  end function distance

end program polymer_analysis
