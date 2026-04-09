PROGRAM Dcmap_generator
    use makeDCmap
    use precision_mod
    IMPLICIT NONE

    INTEGER nmax, nscale, nscale2, npower, npower2, ixmax, ndense, nhypo, stopInd, cutSection
    REAL(preci) ds, dc0, r0
    REAL(preci), DIMENSION(:), ALLOCATABLE :: x0, y0, x0_4_saving, y0_4_saving
    REAL(preci), DIMENSION(:, :), ALLOCATABLE :: dcorg, dc_4_saving
    character(len=256) :: filename_hetero, filename_x0, filename_y0, line, key
    character(len=50)  :: str_var1, str_var2
    character(len=1) :: eq
    integer :: ios
    CHARACTER(*), PARAMETER :: savePath2 = '/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/IA2005/Numerics_with_Fortran/heterogeneity/'

    ! --- read input ---
    open(unit=10, file="input.txt", status="old", action="read")

    do
        read(10, '(A)', iostat=ios) line
        if (ios /= 0) exit   ! end of file

        ! Skip empty or comment lines (optional but recommended)
        if (len_trim(line) == 0) cycle
        if (line(1:1) == "!") cycle

        ! Read the key (left-hand side)
        read(line, *) key

        select case (trim(key))

        case ("nmax")
            read(line, *) key, eq, nmax

        case ("nscale")
            read(line, *) key, eq, nscale

        case ("npower")
            read(line, *) key, eq, npower

        case ("ds")
            read(line, *) key, eq, ds

        case ("ndense")
            read(line, *) key, eq, ndense

        case ("dc0")
            read(line, *) key, eq, dc0

        case ("r0")
            read(line, *) key, eq, r0

        case ("cutSection")
            read(line, *) key, eq, cutSection

        case ("stopInd")
            read(line, *) key, eq, stopInd

        end select
    end do

    close(10)

    ! --- derived quantities ---
    ixmax  = nmax * nscale**npower
    dc0 = dc0 * ds
    r0  = r0  * ds

    nscale2 = nscale/2
    npower2 = npower*2
	nhypo = ndense*(nscale2*nscale2)**npower2

    cutSection = 1 ! whether to cut a section from the full het. map for saving
    stopInd = 255

    ALLOCATE( dcorg(0:ixmax, 0:ixmax) )
    ALLOCATE( dc_4_saving(0:stopInd, 0:stopInd) )
    ALLOCATE( x0(nhypo), y0(nhypo) )

    call make_fractal_DCmap(dcorg, x0, y0, nscale, npower, ndense, ixmax, dc0, r0)

    if (cutSection == 1) then ! there is also some trouble with this when setting nmax=64 - 02.04.2026
        write(*,*) "Cutting a section of the full heterogeneity map for saving..."
        dc_4_saving = dcorg(0:stopInd, 0:stopInd) ! cut appropriate section from dcorg
        x0_4_saving = pack(x0, x0 <= stopInd)
        y0_4_saving = pack(y0, y0 <= stopInd)
    else
        dc_4_saving = dcorg ! save the full map without cutting
        x0_4_saving = x0
        y0_4_saving = y0
    end if

    ! Convert variables to strings
    write(str_var1, '(G0)') ndense
    write(str_var2, '(G0)') nscale

    ! Build full filename
    filename_hetero = savePath2 // 'hetero_' // trim(str_var1) // '_' // trim(str_var2) // '.bin'
    filename_x0 = savePath2 // 'x0_' // trim(str_var1) // '_' // trim(str_var2) // '.bin'
    filename_y0 = savePath2 // 'y0_' // trim(str_var1) // '_' // trim(str_var2) // '.bin'

    open(12, file=filename_hetero, form="unformatted", access="stream")
    write(12) dc_4_saving
    close(12)

    open(12, file=filename_x0, form="unformatted", access="stream")
    write(12) x0_4_saving
    close(12)

    open(12, file=filename_y0, form="unformatted", access="stream")
    write(12) y0_4_saving
    close(12)

END