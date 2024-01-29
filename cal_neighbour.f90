program polymer_analysis
  use omp_lib

  implicit none

  character(len=256) :: infile_1, infile_2, outfile
  integer :: num_threads, t, num_molecules_1, num_molecules_2, i, j, k
  integer total
  integer :: id
  real :: x, y, z
  integer :: ierr, istat_1, istat_2
  integer :: image
  integer :: shift
  integer :: stop_t
  real :: box_size(3)
  real :: ix(3)
  real, allocatable :: coords_1(:,:), coords_2(:,:)
  integer, allocatable :: ids_1(:), ids_2(:)
  integer, allocatable :: ix_1(:,:)
  integer, dimension(:,:), allocatable :: neighbor_list

  total=30
  stop_t=100000000
  shift = (total*2+1)*(total*2+1)*(total*2+1)
  box_size(1)=100
  box_size(2)=100
  box_size(3)=100
  ! Read input arguments
  if (command_argument_count() < 2) then
      print *, "Usage: polymer_analysis <polymer_COM.txt>"
      stop
  endif
  call get_command_argument(1, infile_1)
  call get_command_argument(2, infile_2)
  outfile= "neighbors.txt"

  ! Initialize OpenMP
  num_threads = omp_get_max_threads()
  call omp_set_num_threads(num_threads)

  ! Open output file for writing
  open(unit=10, file=outfile, status='replace')

  ! Open input file for reading
  open(unit=30, file=infile_2, status='old', action='read')
  read(30, *, iostat=istat_2)
  read(30, *, iostat=istat_2)
  read(30, *, iostat=istat_2)
  
  open(unit=20, file=infile_1, status='old', action='read')
  read(20, *, iostat=istat_1)
  read(20, *, iostat=istat_1)
  read(20, *, iostat=istat_1)
  ! Loop over all timesteps
  do

      
      
      ! Read timestep
      read(20, *, iostat=istat_1) t, num_molecules_1
  
      ! Allocate memory for coords and ids
      allocate(coords_1(num_molecules_1,3), ids_1(num_molecules_1),ix_1(num_molecules_1,3))
      
      ! Read coords and ids
      do i = 1, num_molecules_1
          read(20, *, iostat=istat_1) ids_1(i), coords_1(i,1), coords_1(i,2), coords_1(i,3)
      end do


      read(30, *, iostat=istat_2) t, num_molecules_2

      allocate(coords_2(num_molecules_2,3), ids_2(num_molecules_2))
      
      do i = 1, num_molecules_2
          read(30, *, iostat=istat_2) ids_2(i), coords_2(i,1), coords_2(i,2), coords_2(i,3)
      end do
      
      ! Allocate memory for neighbor_list
      allocate(neighbor_list(num_molecules_1,num_molecules_2))
      
      if (istat_2 /= 0) exit
      if (istat_1 /= 0) exit

      ! Initialize neighbor_list
      neighbor_list = 0

      ! Calculate neighbor_list
      !$OMP PARALLEL DO PRIVATE(i,j,ix) SHARED(num_molecules_1,num_molecules_2,coords_1,coords_2,neighbor_list)
      
      do i = 1, num_molecules_1
          do j = 1, num_molecules_2
                   
                  if (distance(coords_1(i,:), coords_2(j,:),ix) < 6.0) then
                      neighbor_list(i,1) = 1
                      ix_1(i,:)= ix
                  endif
   
          end do
      end do
      !$OMP END PARALLEL DO

      ! Write neighbor_list to file
      k=0

      do i = 1, num_molecules_1
          do j = 1, 1
              if ((neighbor_list(i,j) == 1)) then
                  k = k + 1
              endif
          end do
      end do
    if (t<=stop_t) then
    write (10, "(I8)") k
    endif
      do i = 1, num_molecules_1
          do j = 1, 1
              if (neighbor_list(i,j) == 1) then
                  image = (ix_1(i,1)+total)+(ix_1(i,2)+total)*(2*total+1)+(ix_1(i,3)+total)*(2*total+1)*(2*total+1)+1
                  if ((abs(ix_1(i,1))<=total) .and. (abs(ix_1(i,2))<=total) .and. (abs(ix_1(i,3))<=total) .and. (t<=stop_t)) then
                  write(10, "(I9, 2X, I6, 2X, I6)")t, image, ids_1(i)+shift    
                  endif            
              endif
          end do
      end do
      ! Free memory
      deallocate(coords_1,coords_2,ids_1,ids_2, ix_1, neighbor_list)
  end do

  ! Close files
  close(10)
  close(20)

contains

  function distance(a, b, ix) result(dist)
    real, dimension(:), intent(in) :: a, b
    real, dimension(:), intent(out) :: ix
    real :: dist
    integer :: i
    real :: dx(3)
 
    dx = a-b
    ix = 0
    do i=1,3

       do while (abs(dx(i)) > box_size(i)/2.0)
           if (dx(i) <0.0) then 
            dx(i) = dx(i)+box_size(i)
            ix(i) = ix(i)-1
            else
            dx(i) = dx(i)-box_size(i)
            ix(i) = ix(i)+1
            end if
       end do
    end do
           

      
    dist = sqrt(sum(dx**2))
  end function distance

end program polymer_analysis
