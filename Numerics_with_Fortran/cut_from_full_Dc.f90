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
    REAL, DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: dcorg
    REAL(8), INTENT(IN) :: hypocenter_x, hypocenter_y
    INTEGER, INTENT(IN) :: dim_x, dim_y

    INTEGER :: i_start, i_end, j_start, j_end

    write(*,*) "Cutting the appropriate part from the full heterogeneity map..."
    write(*,*) "Hypocenter location (x,y):", hypocenter_x, hypocenter_y
    write(*,*) "Dimensions of non-renormalization domain (x,y):", dim_x, dim_y
    ! Calculate the starting and ending indices for cutting the subarray
    i_start = MAX(1, NINT(hypocenter_x - dim_x / 2))
    i_end = MIN(SIZE(dc_full, 1), NINT(hypocenter_x + dim_x / 2))
    j_start = MAX(1, NINT(hypocenter_y - dim_y / 2))
    j_end = MIN(SIZE(dc_full, 2), NINT(hypocenter_y + dim_y / 2))
    write(*,*) "Made it here 1"

    ! Allocate the output array based on the calculated dimensions
    ALLOCATE(dcorg(i_end - i_start, j_end - j_start)) ! there used to be +1 here.
    write(*,*) "Made it here 2"

    ! Cut the appropriate part from the full heterogeneity map
    dcorg = dc_full(i_start:i_end, j_start:j_end)
    write(*,*) "Made it here 3"

END SUBROUTINE cut_from_full_Dc

END MODULE DcCutter