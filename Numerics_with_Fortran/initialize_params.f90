! Subroutine get_response_function.f
! Sets slip weakening friction law parameters for simulation
! according to specifications given by the user.
! Inputs:
! Outputs: 
! Author: Viktor Essbach (v.essbach@ens.psl.fr), 01/26

Module iniParams
IMPLICIT NONE
CONTAINS
    SUBROUTINE homogeneous_friction(w, tau0, tp, tr, dc, sigma, a, iv, irup, t0, tp0, &
         tr0, ns, ds, rad, nmax, iter, xhypo, yhypo, rini)
        IMPLICIT NONE
        INTEGER :: i, j, nmax, iter, ns
        INTEGER :: iv(nmax, nmax), irup(nmax, nmax)
        DOUBLE PRECISION :: t0, tp0, tr0, ds, rad, xhypo, yhypo, rini
        DOUBLE PRECISION :: w(nmax, nmax), tau0(nmax, nmax), tp(nmax, nmax), &
        tr(nmax, nmax), dc(nmax, nmax), sigma(nmax, nmax), a(nmax, nmax)

        do i = 1, nmax ! loop over all x-positions
            do j = 1, nmax ! loop over all y-positions
                w(i,j) = 0.0d0 ! time integrated slip, set to 0
                tau0(i,j) = t0 ! background stress<tau0
                tp(i,j) = tp0 ! peak strength parameter for cell, set to uniform tp0
                tr(i,j) = tr0 ! residual strength after failure, uniform tr0
                dc(i,j) = dc(i,j)/ns ! renormalized fracture energy (normalized per-scale). Ask Hideo for unit explanation!
                sigma(i,j) = tp(i,j) ! current shear strength?
                a(i,j) = (tp(i,j) - tr(i,j))/dc(i, j) ! slope of linear slip-weakening law
                iv(i,j) = 0 ! integer state variable for cell status during slip: 0 = not slipping, 1 = slipping, 2 = failed
	            irup(i,j) = -1 ! time index for first rupture occurance at given cell, to determine slip time. -1 = not slipped yet
        
                if (iter.eq.0) then ! first rupture initiation if program is at first scale stage!
                    rad = sqrt((i-xhypo)**2 + (j-yhypo)**2)*ds ! get distance to hypocenter for current cell
                        if ( rad.le.rini ) then ! if distance is smaller than initialization radius:
                            tp(i,j) = 0.0d0 ! set cell strength to zero
                            dc(i,j) = 0.0d0 ! set rupture energy to zero
                            a(i,j) = 0.0d0 ! set slope of slip weakening to zero
                            sigma(i,j) = tp(i,j) ! set shear strength to zero
                            iv(i,j) = 2 ! set state variable to "failed!" -> rupture now has initialized at this cell!
	                    endif
	            endif
            enddo
        enddo

    END SUBROUTINE

    SUBROUTINE long_asperity(w, tau0, tp, tr, dc, sigma, a, iv, irup, t0, tp0, &
         tr0, ns, ds, rad, nmax, iter, xhypo, yhypo, rini)
        IMPLICIT NONE
        INTEGER :: i, j, nmax, iter, ns
        INTEGER :: iv(nmax, nmax), irup(nmax, nmax)
        DOUBLE PRECISION :: t0, tp0, tr0, ds, rad, xhypo, yhypo, rini
        DOUBLE PRECISION :: w(nmax, nmax), tau0(nmax, nmax), tp(nmax, nmax), &
        tr(nmax, nmax), dc(nmax, nmax), sigma(nmax, nmax), a(nmax, nmax)

        do i = 1, nmax
            do j = 1, nmax
                tr(i, j) = tr0 ! assign residual stress, uniform everywhere
                if (abs((nmax+1)/2 - j) < (nmax - 1)/2) then ! assign critical stress...
                    tp(i,j) = 10.d0 ! high at boundaries...
                else
                    tp(i,j) = 0.98d0 ! and low inbetween.
                endif

                if (abs((nmax+1)/2 - j) > (nmax-14)/2) then ! assign initial stress
                    tau0(i,j) = 1.d0
                else
                    tau0(i,j) = 0.4d0
                endif

            enddo
        enddo
    END SUBROUTINE
END MODULE iniParams