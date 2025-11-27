! program stoc1-fft-IA05-distrib.f90 (26th October 2018).
! under GNU Lesse General Public License v3.0. on @GitHub.com
!
! Refer Ide and Aochi (JGR, 2005). doi:1O.1029/2004JB003591.
! See README file or related file.
! >ifort stoc1-fft-IA05-distrib.f90 ran1.f fourn-d.f kernel31s_05Avril.f (verified on Unix and Windows)
IMPLICIT NONE
! 
INTEGER nmax, ndata1, ndata2
INTEGER npower, nscale, ixmax, itmx
PARAMETER ( nmax = 64, nscale = 4, npower = 3, ixmax = nmax*nscale**npower )
PARAMETER ( itmx = 500 )
PARAMETER ( ndata1 = 4*nmax, ndata2 = 4*nmax )
REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: vel, vel2
REAL(8), DIMENSION(:, :), ALLOCATABLE :: tau0, tp, tr, &
	stress, sigma, w, a, tau, dc
REAL(8), DIMENSION(:), ALLOCATABLE :: x0, y0, smrate, smoment
REAL(8) :: pi, mu, const, facbiem, facfft, &
		tp0, tr0, dc0, t0, dsreal, dtreal, coef, &
		p000, ker31s, dtau, dsigma, alpha, &
		ds, dt, rad, r0, rini, xhypo, yhypo, &
		xo, yo, r0dum, dcdum, dcmax, dim, smo, mw, &
		t, piece1, piece1_offset, ans, offset
REAL, DIMENSION(:, :), ALLOCATABLE :: dcorg
REAL :: ran1
INTEGER, DIMENSION(:, :), ALLOCATABLE :: iv, irup
INTEGER :: i, j, k, l, m, n, ndir, idata, ix, iy, &
		kmin, iter, kmax, itmx1, icheck, icheck2, it, &
		ndense, nasp, idum, iscale, ihypo, isim, isim0, nhypo, &
		nmax2, ns, i0, j0, i1, j1, k1, nscale2, npower2
INTEGER :: ndata(2)
DOUBLE COMPLEX zdata(ndata1*ndata2), zresp(ndata1*ndata2), zresp_offset(ndata1*ndata2)
DOUBLE COMPLEX zvel(ndata1*ndata2, itmx)
DOUBLE COMPLEX zker(ndata1*ndata2, itmx), zker_offset(ndata1*ndata2, itmx)
DOUBLE COMPLEX zans(ndata1*ndata2)
EXTERNAL ker31s, ran1
CHARACTER*40 name2, name3, name4, name5, name6, name7, name8, name9, dir, param_file
CHARACTER*5  num, num2

! PARAMETER FILE
        param_file = "IA05.prm"
        open(11, file=param_file, status="old", err=99)
        read(11,*) isim0
        close(11)
        
!
! FIXED PARAMETER
!
	pi = acos(-1.0d0)
        facbiem = 2.0d0
	ds = 1.0d0
	dt = ds/facbiem
	kmin = itmx/nscale
! for FFT
	facfft = 1.0d0/(ndata1*ndata2)
	ndata(1) = ndata1
	ndata(2) = ndata2
!
! OUTPUT DIRECTORY
!
	dir = '.'
        ndir = index(dir, ' ')-1
!
! SCALE-INDEPENDENT PARAMETER
!
	mu = 32.40d0
	alpha = 6.0d0
	tp0 = 5.0d0
	tr0 = 0.0d0
	t0 = 3.0d0
	const = sqrt(3.0d0)/(4.0d0*pi)*mu
!
! SCALING PARAMETER
! unit (mu) = mu [GPa]/tb[MPa] dc0[m]/ds [km] = 1
! dc0 should be normalized by 0.001 ds.
!
! INITIAL SETTING
	dc0 = 0.250d0*ds
	r0 = 5.6250d0*ds
	rini = 3.75d0*ds
	ndense = 4
	dim = -2.0
	xhypo = (1+nmax)/2.
	yhypo = (1+nmax)/2.

ALLOCATE( vel(nmax, nmax, 0:itmx), vel2(nmax, nmax, 0:itmx) )
ALLOCATE(tau0(nmax, nmax),     tp(nmax, nmax),   dc(nmax, nmax), &
       stress(nmax, nmax),     tr(nmax, nmax),    a(nmax, nmax), &
     	sigma(nmax, nmax),      w(nmax, nmax), &
     	   iv(nmax, nmax),   irup(nmax, nmax),  tau(nmax, nmax) )
ALLOCATE( smrate(0:itmx), smoment(0:itmx) )
ALLOCATE( dcorg(ixmax, ixmax) )

!
! CREATING DC
! The following section creates the fractal asperity map.
! to be fixed (-411) for obtaining the same result as case1-1
        idum = -411
        dcmax = dc0*nscale**(npower + 1)
	dcorg = real(dcmax)

        nscale2 = nscale/2
        npower2 = npower*2
	nhypo = ndense*(nscale2*nscale2)**npower2
ALLOCATE( x0(nhypo), y0(nhypo) )

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
	name7 = dir(1:ndir)//'/hetero.org'
	ns = ixmax/256
	if(ns.lt.1) ns = 1
	open(71, file=name7)
        do i=1, ixmax, ns
          do j=1, ixmax, ns
            write(71,'(2i7, f10.3)') i, j, dcorg(i,j)
          enddo
        enddo
        close(71)

!
! RESPONSE FUNCTION (for in-plane shear stress calculation)
!
p000 = ker31s(0.0d0, 0.0d0, 0.0d0, 0.0d0, facbiem) ! get reference value of the Kernel at the origin
do k = 1, itmx ! loop over each time step
  do idata=1, ndata1*ndata2 ! set complex response array zresp (1D array storing all 2D grid points)
    zresp(idata) = cmplx(0.0d0, 0.0d0) ! reset to 0 for each time step
  enddo
  do i = 1-nmax, nmax-1 ! loop over first spatial dimension
    ix = i + 1 ! periodic index creation (ask Hideo!)
    if(i.lt.0)  ix = ix + ndata1 ! for negative values of i, shift into positive index values
    do j = 1-nmax, nmax-1 ! inner loop over second spatial dimension
      piece1 = ker31s(dble(i), dble(j), 0.d0, dble(k), facbiem) ! evaluate kernel at location (i, j, 0), time step k
      iy = j + 1
      if(j.lt.0) iy = iy + ndata2 ! periodic indexing for j (index in y-direction)
      idata = ix + (iy-1)*ndata1 ! map 2D index (ix, iy) to 1D vector (to fill zresp)
      zresp(idata) = cmplx(piece1, 0.0d0) ! fill zresp with complex numbers with imaginary part 0 -> builds spatial kernel at time step k
    enddo
  enddo
  call fourn(zresp, ndata, 2, 1) ! call 2D fourier transform routine for zresp
  do idata=1, ndata1*ndata2
    zker(idata, k) = zresp(idata) ! store fourier transformed Green's function for all time steps!
  enddo
enddo

!
! RESPONSE FUNCTION (for of-plane shear stress calculation)
!
offset = 10.d0 ! z-coordinate of off-plane measurement plane
do k = 1, itmx
  do idata=1, ndata1*ndata2
    zresp_offset(idata) = cmplx(0.0d0, 0.0d0)
  enddo
  do i = 1-nmax, nmax-1
    ix = i + 1
    if(i.lt.0) ix = ix + ndata1
    do j = 1-nmax, nmax-1
      piece1_offset = ker31s(dble(i), dble(j), offset, dble(k), facbiem)
      iy = j + 1
      if(j.lt.0) iy = iy + ndata2
      idata = ix + (iy-1)*ndata1
      zresp_offset(idata) = cmplx(piece1_offset, 0.0d0)
    enddo
  enddo
  call fourn(zresp, ndata, 2, 1)
  do idata=1, ndata1*ndata2
    zker_offset(idata, k) = zresp_offset(idata)
  enddo
enddo



name8 = dir(1:ndir)//'/kernel.dat'

open(18, file=name8)
do idata = 1, ndata1*ndata2
  write(18, '(100g15.5)') zker(idata,:)
enddo
close(18)

name9 = dir(1:ndir)//'/kernel_offset.dat'

open(19, file=name9)
do idata = 1, ndata1*ndata2
  write(19, '(100g15.5)') zker_offset(idata,:)
enddo
close(19)

name2 = dir(1:ndir)//'/hoge2.dat'

open(12, file=name2)
write(12, '(5i10)') nmax, nscale, npower, ixmax, itmx
write(12, '(5f10.3)') mu, alpha, tp0, tr0, t0
write(12, '(4f10.3)') dc0, dcmax, r0, rini
write(12, '(4i10)') ndense, nscale2, npower2, nhypo
close(12)

!!
!! ITERATION OF HYPOCENTER LOCATION
!!
!! isim0 red by a parameter file
do isim = isim0, isim0
  ihypo = isim
!! test 
  if(isim.eq.0) ihypo = 1
!! two scenarios of Mw3.8 for Aochi & Burnol (2018)
  if(isim.eq.1) ihypo = 537
  if(isim.eq.2) ihypo = 546
! for checking large events
  if(isim.eq.1) ihypo = 806
  if(isim.eq.2) ihypo = 7987 
  if(isim.eq.3) ihypo = 9141
  if(isim.eq.4) ihypo = 10426
  if(isim.eq.5) ihypo = 12746 
  if(isim.eq.6) ihypo = 13375

  write(num, '(i5.5)') ihypo 
  name6 = dir(1:ndir)//'/output'//num(1:5)//'i.dat'
  open(16, file=name6)
  write(16, '(i10, 2f10.3)') ihypo, x0(ihypo), y0(ihypo)
  write(16, '(f10.3)') rini
  close(16)

!!	
!! ITERATION OF STAGE
!!
  STAGE: do iter = 0, npower ! loop over different scales
    kmax = itmx
    ns = nscale**iter ! scale factor determines how many original grid cells are aggregated into one cell. nscale = 4 by default

! RENORMALIZATION
    nmax2 = nmax*nscale**iter ! nmax: number of grid cells in coarse grid, nmax2: "physical" size of region that coarse grid spans, in elementary grid cells
    ! remember: nmax2 = nmax * ns
    do i = 1, nmax ! loop over all grid cells in coarse grid, in x direction
      i0 = ixmax/2 - nmax2/2 + ns*(i-1) + 1 ! find x-index of upper left corner of the coarse cell with index (i,j), in terms of elementary grid cells
      do j = 1, nmax
        j0 = ixmax/2 - nmax2/2 + ns*(j-1) + 1 ! find y-index of upper left corner of the coarse cell with index (i,j), in terms of elementary grid cells
        dc(i,j) = 0.0d0 ! initialize coarse grid value for fracture energy
        do i1 = 1, ns
          do j1 = 1, ns
            ix = i0 + i1 - 1 + int(x0(ihypo) - (ixmax+1)/2.) ! shift center of coarse graining step to hypocenter location, get proper x-index in terms of elementary cells
            if(ix.lt.1) ix = ixmax + ix ! periodic boundary conditions?
            if(ix.gt.ixmax) ix =  ix - ixmax
            iy = j0 + j1 - 1 + int(y0(ihypo) - (ixmax+1)/2.) ! shift center of coarse graining step to hypocenter location, get proper x-index in terms of elementary cells
            if(iy.lt.1) iy = ixmax + iy ! periodic boundary conditions?
            if(iy.gt.ixmax) iy =  iy - ixmax

            dc(i,j) = dc(i,j) + dble(dcorg(ix, iy)) ! now, sum up contributions from all elementary cells in coarse grid cell and...
          enddo
        enddo
        dc(i,j) = dc(i,j)/(ns**2) ! average them using the count of elementary cells in said grid cell!
      enddo
    enddo

! INITIALIZATION OF PARAMETERS
    vel = 0.0d0 ! vel(i,j,k): slip velocity at every coarse grid-cell and time set to 0
    smrate = 0.0d0 ! smrate(k): moment release rate at every time set to 0

    do i = 1, nmax ! loop over all x-positions
      do j = 1, nmax ! loop over all y-positions
        w(i,j) = 0.0d0 ! time integrated slip, set to 0
        tau0(i,j) = t0 ! background stress<ftau0
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

    if( iter.ne.0 ) then ! after first scale stage:
      kmax = itmx
      if( iter.eq.npower) kmax = itmx
      kmin = itmx1/nscale ! this has to do with continuing the calculation at a certain time step after changing scales, I think

!    do k = 1, itmx
      do k = 1, nscale*(itmx1/nscale) ! maybe transfer results from previous finer scale to coarser scale? ask Hideo
        n = (k-1)/nscale + 1
        do i = 1, nmax
          l = (i-1)/nscale + 1 + (nscale-1)*nmax/(2*nscale) ! mapping from (i,j,k) indices to (l,m,n) indices
          do j = 1, nmax
            m = (j-1)/nscale + 1 + (nscale-1)*nmax/(2*nscale)
            if( vel2(i,j,k).ne.0. )then
              vel(l,m,n) = vel(l,m,n) + &
		vel2(i,j,k)/(nscale*nscale*nscale) ! transfer slip velocities temporarily saved in vel2 (see below) to new indices
            endif
          enddo
        enddo
      enddo
    endif

    icheck = 0 ! later used to check if rupture reached domain boundaries
    itmx1 = kmax ! effective time limit for current iteration stage
    smoment = 0.0d0 ! seismic moment will be accumulated in smoment(k)

! ITERATION OF TIME
    TIME: do k = 1, kmax ! this is the main time loop, central piece of this code!!!
      if( mod(k, 20).eq.0 ) write(6,'(a20, 3i5)') "SIMULATION", ihypo, iter, k ! write hypocenter location, scale stage, and time step to terminal
      icheck2 = 0 ! tracks if any cell is slipping at current time-step

      if ( k.ne.1 ) then ! if the first time step has passed...
        do idata=1, ndata1*ndata2 ! ... prepare to compute convolution with velocity history via FFT
          zdata(idata) = cmplx(0.0d0, 0.0d0) ! velocity history will be stored here
        enddo
        do j = 1, nmax
          do i = 1, nmax
            idata = i + (j-1)*ndata1 ! transform (i,j) to single index for 1D vector
            zdata(idata) = cmplx(dble(vel(i,j,k-1)), 0.0d0) ! fill zdata with slip velocities from previous time step
          enddo
        enddo
        call fourn(zdata, ndata, 2, 1) ! now, run FFT with slip velo history
        do idata=1, ndata1*ndata2
          zvel(idata, k-1) = zdata(idata) ! zvel contains FFT results for all time steps
        enddo

        do idata=1, ndata1*ndata2
          zans(idata) = (0.0d0, 0.0d0) ! set zans to 0 
          do n=1, k-1 ! this is the actual convolution, carried out over all previous time steps -> over the full slip velo history
            zans(idata) = zans(idata) +  &
		zker(idata, k-n)*zvel(idata, n) ! multiply kernel with slip velo histories in frequency domain and sum up -> convolution in space domain
          enddo
        enddo
        call fourn(zans, ndata, 2, -1) ! FFT backtransform to obtain stress change caused by past slip action
      endif

      do i=1, nmax
        do j=1, nmax
          idata = i + (j-1)*ndata1 ! map (i,j) to 1D vector array again

          if( k.ne.1 ) then ! again for all time steps other than the first...
            ans = facfft*dble(zans(idata)) ! get the result of the FFT backtransform by multiplying with a normalization factor.
            dtau = tau0(i,j) + const*ans ! now obtain driving shear stress on cell at (i,j) at time step k by adding the elastic stress perturbation to the background stress
          else ! remember: const = sqrt(3.0d0)/(4.0d0*pi)*mu is a factor to get the stress units right
  	    dtau = tau0(i,j) ! in first time step, dtau is simply the background stress
          endif
          dsigma = sigma(i,j) ! current shear strength at cell (i,j), to be compared with dtau

	  if( iter.ne.0.and.k.le.kmin ) then ! special case exception for non-zero iterations
	    if( vel(i,j,k).ne.0. ) then ! if velo at (i,j,k) is not 0...
	      iv(i,j) = 1 ! assign the cell as "slipping"
	    else
	      iv(i,j) = 0 ! else, assign it as "locked", except if dtau (driving stress) exceeds dsigma (shear strength)
	      if( dtau.gt.dsigma ) iv(i,j) = 2 ! in that case: mark as failed/nucleated
	    endif
	  else
            if( iv(i,j).eq.0 ) then ! this is the main velocity update, core of the code! if cell is assigned as locked...
              vel(i,j,k) = 0.0d0 ! ... set its slip velo to 0.
              if( dtau.gt.dsigma ) iv(i,j) = 2 ! but if driving stress exceeds strength, set it as failing
            else
              iv(i,j) = 1 ! if cell is assigned as slipping...
	      if( (const*p000 + a(i,j)*dt).ge.0. ) then ! ... and this criterion holds (ask Hideo. Cell reached residual strength?),
	        vel(i,j,k) = (tr(i,j) - dtau)/(const*p000) ! assign slip velocity according to velo = (residual stress - driving stress)/unit factor (?)
	      else
                vel(i,j,k)  = (dsigma - dtau)/(const*p000 + a(i,j)*dt) ! if the criterion does not hold, assign velo according to slightly altered law
                if((w(i,j)+vel(i,j,k)*dt).gt.dc(i,j)) then
                  vel(i,j,k) = ( tr(i,j) - dtau)/(const*p000)
                endif
	      endif
              if(vel(i,j,k).lt.0.) then ! if slip velo becomes negative, set it to 0 and assign the cell as locked
                vel(i,j,k)  = 0.0d0
                iv(i, j) = 0
                if( dtau.gt.dsigma ) iv(i,j) = 2 ! ..., but if dtau > dsigma reassign it as failing
              endif
            endif
	  endif
          ! now, update stress, slip, etc
          stress(i,j) = dtau + const*p000*vel(i,j,k) ! instantaneous stress? Ask Hideo
          tau(i,j) = dtau ! store driving stress for later use
          w(i,j) = w(i,j) + vel(i,j,k)*dt ! Euler time integration of slip to obtain w at cell
	  if( iter.ne.npower ) vel2(i,j,k) = vel(i,j,k) ! if the final scale stage has not been reached, store velocity field in vel2
          if( w(i,j).gt.dc(i,j) ) then ! if accumulated slip exceeds fracture energy...
            sigma(i,j) = tr(i,j) ! ..., sigma becomes residual stress
          else
            sigma(i,j) = tp(i,j) - a(i,j)*w(i,j) ! if not, it is updated according to the linear slip weakening with slide from tp to tr with slope a.
          endif
	  if( irup(i,j).lt.0.and.vel(i,j,k).ne.0.) irup(i,j) = k ! track which cells have been reached by rupture and assign rupture time(index)

          if((i.eq.1.or.i.eq.nmax).and.vel(i,j,k).ne.0.) icheck = 1 ! if cell at edge of scale stage is reached by rupture, set icheck to 1
          if((j.eq.1.or.j.eq.nmax).and.vel(i,j,k).ne.0.) icheck = 1
	  if( vel(i,j,k).ne.0.) icheck2 = 1 ! if there is a rupture at all, set icheck2 to 1 (if icheck2 = 0, rupture has died out)

	  smrate(k) = smrate(k) + vel(i,j,k)*alpha ! update moment release rate by summing up slip velos from all cells, scaled by alpha (ask Hideo)
          smoment(k) = smoment(k) + w(i,j)*ns ! update total released moment with time integrated slip, scaled by ns
        enddo
      enddo

!      if( iter.le.npower ) then ! if final scale stage has not been reached, write output files every time step (took this out)
!	write(num, '(i5.5)') ihypo
!        write(num2,'(i4.4)') iter*1000 + k
!        name3 = dir(1:ndir)//'/'//num(1:5)//'_step'//num2(1:4)//'.dat'
!        open(13, file=name3)
	   
!        do i = 1, nmax
!          do j = 1, nmax
!            write(13, '(f7.1, 1x, f7.1, 1x, 3f13.8)') &
!     		real(i), real(j), vel(i,j,k)*alpha, w(i,j)*ns, stress(i,j)
!          enddo
!        enddo
!        close(13)
!      endif

      if( k.eq.itmx.or.(iter.ne.npower.and.icheck.eq.1).or.icheck2.eq.0) then ! if maximum iterations are reached, or rupture has reached boundary of largest scale, or rupture has died out...
        write(num,'(i5.5)') ihypo ! write final results to output files
        write(num2,'(i1.1)') iter
        name3 = dir(1:ndir)//'/output'//num(1:5)//num2(1:1)//'.dat'
        open(13, file=name3)

	smo = 0.0d0
        do i = 1, nmax
          do j = 1, nmax
            write(13, 105) i, j, dc(i,j)*ns, irup(i,j), w(i,j)*ns
	    if(w(i,j).ne.0.) smo = smo + w(i,j)*ns
 105	  format(2i5, 1x, f9.3, i10, 1x, f9.3)
          enddo
        enddo
        close(13)

        coef = (0.4)**3*(ds*ns)**2*mu*10.0**9
	dsreal = 4.d0*ns*ds
	dtreal = dsreal/(2.*alpha*1000.)

	smo = 4.*(smo/1000.)*dsreal**2*mu*10.0**9
	mw = (log10(smo)-9.1)/1.5 ! calculate the magnitude of the event
        name4 = dir(1:ndir)//'/output'//num(1:5)//'f.dat'
        open(14, file=name4)
	write(14, '(3i5)') ihypo, iter, k
	write(14, '(e10.4, 2f10.3, f10.4)') smo, mw, dsreal, dtreal
        do it=1, k
          write(14, '(i5,f12.4,2e12.4)') it, it*dtreal, smoment(it)*coef, &
                  (smoment(it)-smoment(it-1))*coef/dtreal
        enddo
	close(14)

	name5 = dir(1:ndir)//'/moment'//num(1:5)//num2(1:1)//'.dat'
	open(15, file=name5)
	smo = 0.0d0
	do k1 = 1, k
	  smo = smo + smrate(k1)*dtreal*dsreal**2*mu*10.0**9
	  mw =  (log10(smo)-9.1)/1.5
	  write(15, '(f15.6, e15.5, e15.5, f15.5)') &
		k1*dtreal, smrate(k1)*dsreal**2*mu*10.0**9, smo, mw
	enddo
	close(15)

      endif

      if(iter.ne.npower.and.icheck.eq.1) then
        itmx1 = k - 1
        write(6,*) "returning at ...", itmx1
        exit TIME
      endif

      if(icheck2.eq.0) then
	write(6,*) "rupture terminated", iter, k
	exit STAGE
      endif

    enddo TIME
  enddo STAGE
enddo

 99     continue
        write(*,*) "END OF SIMULATION"

END

