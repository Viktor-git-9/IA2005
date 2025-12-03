! Subroutine get_response_function.f
! Calls kernel functions to create fourier transformation of response function 
! according to specifications given by the user.
! Inputs:
! Outputs: 
! Author: Viktor Essbach (v.essbach@ens.psl.fr), 12/25
Module getResponse
IMPLICIT NONE
CONTAINS
    SUBROUTINE get_resp(p000, zker, itmx, ndata1, ndata2, nmax, offset, facbiem, stressType)
        IMPLICIT NONE
        INTEGER :: k, idata, i, j, ix, iy, itmx, ndata1, ndata2, nmax, stressType
        INTEGER :: ndata(2)
        DOUBLE COMPLEX :: zker(ndata1*ndata2, itmx), zresp(ndata1*ndata2)
        DOUBLE PRECISION :: p000, facbiem, piece1, offset
        interface
            REAL(8) FUNCTION ker31s(c1,c2,c3,t,w)
            REAL(8), intent(in) :: c1,c2,c3,t,w
            end function ker31s

            REAL(8) FUNCTION ker33s(c1,c2,c3,t,w)
            REAL(8), intent(in) :: c1,c2,c3,t,w
            end function ker33s
        end interface
        !
        ! RESPONSE FUNCTION (for of-plane shear stress calculation)
        !
        ndata(1) = ndata1
	    ndata(2) = ndata2
        select case(stressType)
        case (31)
                p000 = ker31s(0.0d0, 0.0d0, offset, 0.0d0, facbiem)
                do k = 1, itmx
                    do idata=1, ndata1*ndata2
                    zresp(idata) = cmplx(0.0d0, 0.0d0)
                    enddo
                    do i = 1-nmax, nmax-1
                        ix = i + 1
                        if(i.lt.0) ix = ix + ndata1
                        do j = 1-nmax, nmax-1
                            piece1 = ker31s(dble(i), dble(j), offset, dble(k), facbiem)
                            iy = j + 1
                            if(j.lt.0) iy = iy + ndata2
                            idata = ix + (iy-1)*ndata1
                            zresp(idata) = cmplx(piece1, 0.0d0)
                        enddo
                    enddo
                    call fourn(zresp, ndata, 2, 1)
                    do idata=1, ndata1*ndata2
                        zker(idata, k) = zresp(idata)
                    enddo
                enddo

            case (33)
                p000 = ker33s(0.0d0, 0.0d0, offset, 0.0d0, facbiem)
                do k = 1, itmx
                    do idata=1, ndata1*ndata2
                    zresp(idata) = cmplx(0.d0, 0.d0)
                    enddo
                    do i = 1-nmax, nmax-1
                        ix = i + 1
                        if(i.lt.0) ix = ix + ndata1
                        do j = 1-nmax, nmax-1
                            piece1 = ker33s(dble(i), dble(j), offset, dble(k), facbiem)
                            iy = j + 1
                            if(j.lt.0) iy = iy + ndata2
                            idata = ix + (iy-1)*ndata1
                            zresp(idata) = cmplx(piece1, 0.0d0)
                        enddo
                    enddo
                    call fourn(zresp, ndata, 2, 1)
                    do idata=1, ndata1*ndata2
                        zker(idata, k) = zresp(idata)
                    enddo
                enddo
            end select
            return

    END SUBROUTINE
END MODULE getResponse
                
