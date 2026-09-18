! Subroutine to cut appropriate part from the full heterogeneity map
! depending on hypocenter location and dimensions of the non-renormalization domain.
! The output of this subroutine takes on the role of dcorg.
! Author: Viktor Essbach (essbach@ens.psl.fr), 02/2026
MODULE DcCutter
    IMPLICIT NONE
    CONTAINS

SUBROUTINE cut_from_full_Dc(dc_full, dcorg, hypocenter_x, hypocenter_y, dim_x, dim_y)
    IMPLICIT NONE
    REAL(8), DIMENSION(:, :), INTENT(IN) :: dc_full
    REAL(8), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: dcorg
    REAL(8), INTENT(IN) :: hypocenter_x, hypocenter_y
    INTEGER, INTENT(IN) :: dim_x, dim_y

    INTEGER :: i_start, i_end, j_start, j_end, renormDim_x, renormDim_y

    write(*,*) "Cutting the appropriate part from the full heterogeneity map..."
    write(*,*) "Hypocenter location (x,y):", hypocenter_x, hypocenter_y
    ! Calculate the starting and ending indices for cutting the subarray
    ! Setup to center around upper left corner of grid cell containing hypocenter
    ! Approach taken from renormalize_Dcmap.f90.

    renormDim_x = 4096  ! this needs to be adapted to the size of the full renormalization domain,
    renormDim_y = 4096  ! which with nmax=64, nscale=4, npower=3 is 4096x4096

    i_start = renormDim_x/2 - dim_x/2 + INT(hypocenter_x - real(renormDim_x+1,8)/2.0d0) + 1
    if(i_start < 1) i_start = i_start + renormDim_x
    if(i_start > renormDim_x) i_start = i_start - renormDim_x
    i_end = i_start + dim_x - 1

    j_start = renormDim_y/2 - dim_y/2 + INT(hypocenter_y - real(renormDim_y+1,8)/2.0d0) + 1
    if(j_start < 1) j_start = j_start + renormDim_y
    if(j_start > renormDim_y) j_start = j_start - renormDim_y
    j_end = j_start + dim_y - 1

    ! write(*,*) "Calculated cutting indices: i_start =", i_start, ", i_end =", i_end, &
    !            ", j_start =", j_start, ", j_end =", j_end

    ! Allocate the output array based on the calculated dimensions
    ALLOCATE(dcorg(dim_x, dim_y)) ! there used to be +1 here.

    ! Cut the appropriate part from the full heterogeneity map
    dcorg = dc_full(i_start:i_end, j_start:j_end)
    write(*,*) "Cutting done!"

END SUBROUTINE cut_from_full_Dc

SUBROUTINE cut_from_full_Dc_old(dc_full, dcorg, hypocenter_x, hypocenter_y, dim_x, dim_y)
    IMPLICIT NONE
    REAL(8), DIMENSION(:, :), INTENT(IN) :: dc_full
    REAL(8), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: dcorg
    REAL(8), INTENT(IN) :: hypocenter_x, hypocenter_y
    INTEGER, INTENT(IN) :: dim_x, dim_y

    INTEGER :: i_start, i_end, j_start, j_end

    write(*,*) "Cutting the appropriate part from the full heterogeneity map..."
    write(*,*) "Hypocenter location (x,y):", hypocenter_x, hypocenter_y
    ! Calculate the starting and ending indices for cutting the subarray
    ! Setup to center around upper left corner of grid cell containing hypocenter
    !i_start = MAX(1, FLOOR(hypocenter_x - dim_x / 2))
    !i_end = MIN(SIZE(dc_full, 1), FLOOR(hypocenter_x + dim_x / 2))
    !j_start = MAX(1, CEILING(hypocenter_y - dim_y / 2))
    !j_end = MIN(SIZE(dc_full, 2), CEILING(hypocenter_y + dim_y / 2))

    i_start = MAX(1, INT(hypocenter_x - dim_x / 2.0d0))
    ! this is to emulate the behaviour of renormalize_Dcmap.f90
    ! there, a wraparound is performed when the following value is smaller 0
    ! value 4096 is hard coded at the moment -> size of the full renormalization domain
    if ( INT(hypocenter_x - 4096 / 2.0d0) < 0 ) then
        write(*,*) "Adjusted i_start"
        i_start = i_start + 1
    end if
    i_end = MIN(SIZE(dc_full, 1), i_start + dim_x)

    j_start = MAX(1, CEILING(hypocenter_y - dim_y / 2.0d0))
    if ( INT(hypocenter_y - 4096 / 2.0d0) < 0 ) then
        write(*,*) "Adjusted j_start from ", j_start, " to ", j_start + 1
        j_start = j_start + 1
    end if
    j_end = MIN(SIZE(dc_full, 2), j_start + dim_y)

    write(*,*) "Calculated cutting indices: i_start =", i_start, ", i_end =", i_end, &
               ", j_start =", j_start, ", j_end =", j_end

    ! Allocate the output array based on the calculated dimensions
    ALLOCATE(dcorg(dim_x, dim_y)) ! there used to be +1 here.
    write(*,*) "Allocated dcorg with dimensions: ", SIZE(dcorg, 1), " x ", SIZE(dcorg, 2)

    ! Cut the appropriate part from the full heterogeneity map
    dcorg = dc_full(i_start:i_end-1, j_start:j_end-1)
    write(*,*) "Cutting done!"

END SUBROUTINE cut_from_full_Dc_old

END MODULE DcCutter