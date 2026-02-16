MODULE renormalizeDcmap
    IMPLICIT NONE
    Contains

    SUBROUTINE renormalize_dcmap(dc, dcorg, x0, y0, ihypo, nmax, nmax2, ixmax, ns)
        IMPLICIT NONE
        REAL(8), intent(out) :: dc(:,:)
        REAL(8), intent(in) :: x0(:), y0(:)
        REAL(8), intent(in) :: dcorg(:,:)
        INTEGER, intent(in) :: ihypo, nmax, nmax2, ixmax, ns
        INTEGER :: i, j, i0, j0, i1, j1, ix, iy

         do i = 1, nmax ! loop over all grid cells in coarse grid, in x direction
            i0 = ixmax/2 - nmax2/2 + ns*(i-1) + 1 ! find x-index of upper left corner of the coarse cell with index (i,j), in terms of elementary grid cells
            do j = 1, nmax
               j0 = ixmax/2 - nmax2/2 + ns*(j-1) + 1 ! find y-index of upper left corner of the coarse cell with index (i,j), in terms of elementary grid cells
               dc(i,j) = 0.0d0 ! initialize coarse grid value for fracture energy
               do i1 = 1, ns
                  do j1 = 1, ns
                     ix = i0 + i1 - 1 + int(x0(ihypo) - real(ixmax+1,8)/2.0d0) ! shift center of coarse graining step to hypocenter location, get proper x-index in terms of elementary cells
                     if(ix.lt.1) ix = ixmax + ix ! periodic boundary conditions?
                     if(ix.gt.ixmax) ix =  ix - ixmax
                     iy = j0 + j1 - 1 + int(y0(ihypo) - real(ixmax+1,8)/2.0d0) ! shift center of coarse graining step to hypocenter location, get proper x-index in terms of elementary cells
                     if(iy.lt.1) iy = ixmax + iy ! periodic boundary conditions?
                     if(iy.gt.ixmax) iy =  iy - ixmax

                     dc(i,j) = dc(i,j) + dble(dcorg(ix, iy)) ! now, sum up contributions from all elementary cells in coarse grid cell and...
                  enddo
               enddo
               dc(i,j) = dc(i,j)/(ns**2) ! average them using the count of elementary cells in said grid cell!
            enddo
         enddo

    END SUBROUTINE renormalize_dcmap
END MODULE renormalizeDcmap