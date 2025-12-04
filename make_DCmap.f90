MODULE makeDCmap
IMPLICIT NONE
Contains

SUBROUTINE make_fractal_DCmap(dcorg, x0, y0, nscale, npower, ndense, ixmax, dc0, r0)
    ! Description goes here
    IMPLICIT NONE
    INTEGER :: nscale, npower, idum, iscale, ihypo, i, j, i1, j1, &
    nhypo, ndense, nscale2, npower2, nasp, ixmax
    REAL :: dcorg(:,:)
    REAL(8) :: x0(:), y0(:)
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

        do iscale = 0,  npower2
          nasp = ndense*(nscale2*nscale2)**(npower2 - iscale)
          r0dum = r0*nscale2**iscale
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
      END SUBROUTINE
    END MODULE makeDCmap

