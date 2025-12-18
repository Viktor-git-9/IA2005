MODULE makeDCmap
IMPLICIT NONE
Contains

SUBROUTINE make_fractal_DCmap(dcorg, x0, y0, nscale, npower, ndense, ixmax, dc0, r0)
  ! Description goes here
  IMPLICIT NONE
  INTEGER :: nscale, npower, idum, iscale, ihypo, i, j, i1, j1, &
  nhypo, ndense, nscale2, npower2, nasp, ixmax
  REAL, intent(inout) :: dcorg(:,:)
  REAL(8), intent(out), allocatable :: x0(:), y0(:)
  REAL(8) :: dcmax, dc0, r0, r0dum, dcdum, xo, yo, rad

  interface
    REAL FUNCTION ran1(idhypo)
    INTEGER, intent(in) :: idhypo
    END FUNCTION ran1
  end interface


  idum = -411
  dcmax = dc0*nscale**(npower + 1)
	dcorg = real(dcmax)

  nscale2 = nscale/2
  npower2 = npower*2
	nhypo = ndense*(nscale2*nscale2)**npower2
  write(*,*) nhypo
  write(*,*) "Creating fractal asperity map..."
  ALLOCATE( x0(nhypo), y0(nhypo) )

        do iscale = 0,  npower2
          nasp = ndense*(nscale2*nscale2)**(npower2 - iscale)
          r0dum = r0*nscale2**iscale
          write(*,*) r0dum
          dcdum = dc0*nscale2**iscale
          do ihypo = 1, nasp
            xo = ran1(idum)*ixmax
            yo = ran1(idum)*ixmax
            if(iscale.eq.0 ) then
                x0(ihypo) = xo
                y0(ihypo) = yo
            endif

            do i = int(xo - r0dum)-1, int(xo + r0dum)+1
              i1 = i
              if(i1.lt.1) i1 = ixmax + i
              if(i1.gt.ixmax) i1 =  i - ixmax
              do j = int(yo - r0dum)-1, int(yo + r0dum)+1
                j1 = j
                if(j1.lt.1) j1 = ixmax + j
                if(j1.gt.ixmax) j1 =  j - ixmax

                rad = sqrt((i-xo)**2 + (j-yo)**2)
                if( rad.le.r0dum ) then
                  if( dcorg(i1,j1).gt.dcdum ) dcorg(i1,j1) = real(dcdum)
                endif
	      enddo
	    enddo

          enddo
        enddo
write(*,*) "Fractal asperity map ready."
END SUBROUTINE

SUBROUTINE make_homogeneous_DCmap(dcorg, x0, y0, ixmax, dc0, dcmax, r_asperity, ihypo)
  IMPLICIT NONE
  INTEGER :: i, j, ixmax, ihypo
  REAL(8) :: dcmax, dc0, r_asperity, rad
  REAL, intent(inout) :: dcorg(:,:)
  REAL(8), intent(out), allocatable :: x0(:), y0(:)

  ALLOCATE( x0(ixmax-2), y0(ixmax-2) )

  write(*,*) "Creating hom. asperity map."
  do i = 1, ixmax
    do j = 1, ixmax
      if(i.gt.0.and.i.lt.ixmax) x0(i) = i - 0.5
      if(j.gt.0.and.j.lt.ixmax) y0(j) = j - 0.5
    enddo
  enddo
  dcorg(:,:) = dc0
  do i = 0, ixmax
    do j = 0, ixmax
      rad = sqrt((i-x0(ihypo))**2 + (j-y0(ihypo))**2)
      !if(rad.lt.r_asperity) dcorg(i,j) = dc0 ! this can be used to insert a single circular asperity.
    enddo
  enddo
  write(*,*) "Hom. asperity map done!"

END SUBROUTINE

END MODULE makeDCmap

