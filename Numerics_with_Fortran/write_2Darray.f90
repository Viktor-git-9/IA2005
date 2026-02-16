Module writeArrays
IMPLICIT NONE
CONTAINS

SUBROUTINE write_real_1DArray(data, filename)
    ! Subroutine write_real_1Darray
    ! Write the indices and elements of a real double precision 1D array to a .dat file.
    ! Careful, this is slow for large arrays. Use for debugging.
    ! Inputs: data (1D array of real double precision numbers), filename (character string to designate file)
    ! Outputs: writes data to file
    ! Author: Viktor Essbach (v.essbach@ens.psl.fr), 02/26
    IMPLICIT NONE
    INTEGER :: i, nRows
    DOUBLE PRECISION, intent(in) :: data(:)
    CHARACTER(len=*), intent(in) :: filename
    write(*,*) "Writing to file..."
    nRows = size(data)
    open(19, file=filename)
    do i = 1, nRows
        write(19, 100) i, data(i)
    enddo
    close(19)
    100 format(i5, 1x, f15.6)
    write(*,*) "File written!"

END SUBROUTINE

SUBROUTINE write_real_2DArray(data, filename)
    ! Subroutine write_real_2Darray
    ! Write the indices and elements of a real double precision 2D array to a .dat file.
    ! Careful, this is slow for large arrays. Use for debugging.
    ! Inputs: data (2D array of real double precision numbers), filename (character string to designate file)
    ! Outputs: writes data to file
    ! Author: Viktor Essbach (v.essbach@ens.psl.fr), 12/25
    IMPLICIT NONE
    INTEGER :: i, j, nRows, nCols
    DOUBLE PRECISION, intent(in) :: data(:,:)
    CHARACTER(len=*), intent(in) :: filename

    write(*,*) "Writing to file..."
    nRows = size(data,1)
    nCols = size(data,2)
    open(19, file=filename)
    do i = 1, nRows
        do j = 1, nCols
            write(19, 100) i, j, data(i,j)
        enddo
    enddo
    close(19)
    100 format(i5, 1x, i5, 1x, f15.6)
    write(*,*) "File written!"

END SUBROUTINE

SUBROUTINE write_cmplx_2DArray(data, filename)
    ! Subroutine write_cmplx_2Darray
    ! Write the indices and elements of a complex double precision 2D array to a .dat file.
    ! Careful, this is slow for large arrays. Use for debugging.
    ! Inputs: data (2D array of complex double precision numbers), filename (character string to designate file)
    ! Outputs: writes data to file
    ! Author: Viktor Essbach (v.essbach@ens.psl.fr), 12/25
    IMPLICIT NONE
    INTEGER :: i, j, nRows, nCols
    DOUBLE COMPLEX, intent(in) :: data(:,:)
    CHARACTER(len=*), intent(in) :: filename

    write(*,*) "Writing to file..."
    nRows = size(data,1)
    nCols = size(data,2)
    open(19, file=filename)
    do i = 1, nRows
        do j = 1, nCols
            write(19, 100) i, j, real(data(i,j)), aimag(data(i,j))
        enddo
    enddo
    close(19)
    100 format(i5, 1x, i5, 1x, f15.6, 1x, f15.6)
    write(*,*) "File written!"

END SUBROUTINE

SUBROUTINE write_real_2DArray_bin(data, filename)
    IMPLICIT NONE
    DOUBLE PRECISION, intent(in) :: data(:,:)
    CHARACTER(len=*), intent(in) :: filename

    write(*,*) "Writing real 2D to file ", filename
    write(*,*) "Dimensions of data array:", size(data,1), size(data,2), size(data)
    open(unit=19, file=filename, form="unformatted", access="stream")
    write(19) data
    close(19)
    write(*,*) "File written!"

END SUBROUTINE

SUBROUTINE write_cmplx_2DArray_bin(data, filename)
    IMPLICIT NONE
    DOUBLE COMPLEX, intent(in) :: data(:,:)
    CHARACTER(len=*), intent(in) :: filename

    write(*,*) "Writing cx 2D to file ", filename
    write(*,*) "Dimensions of data array:", size(data,1), size(data,2), size(data)
    open(unit=19, file=filename, form="unformatted", access="stream")
    write(19) data
    close(19)
    write(*,*) "File written!"

END SUBROUTINE

SUBROUTINE write_real_3DArray_bin(data, filename)
    IMPLICIT NONE
    DOUBLE PRECISION, intent(in) :: data(:,:,:)
    CHARACTER(len=*), intent(in) :: filename

    write(*,*) "Writing real 3D to file ", filename
    write(*,*) "Dimensions of data array:", size(data,1), size(data,2), size(data,3), size(data)
    open(unit=19, file=filename, form="unformatted", access="stream")
    write(19) data
    close(19)
    write(*,*) "File written!"

END SUBROUTINE

SUBROUTINE write_cmplx_3DArray_bin(data, filename)
    IMPLICIT NONE
    DOUBLE COMPLEX, intent(in) :: data(:,:,:)
    CHARACTER(len=*), intent(in) :: filename

    write(*,*) "Writing cx 3D to file ", filename
    write(*,*) "Dimensions of data array:", size(data,1), size(data,2), size(data,3), size(data)
    open(unit=19, file=filename, form="unformatted", access="stream")
    write(19) data
    close(19)
    write(*,*) "File written!"

END SUBROUTINE

END MODULE writeArrays
