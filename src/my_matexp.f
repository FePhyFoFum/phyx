      subroutine wrapalldmexpv(n,m,t, v, w, tol, anorm, wsp, lwsp, iwsp,
     .                       liwsp, itrace, iflag,ia,ja,a,nz,res )
      implicit none

      integer,intent(inout) :: m, lwsp, liwsp, itrace, iflag
      integer,intent(inout) :: iwsp(liwsp)
      double precision,intent(inout) :: t, tol, anorm, v(n), w(n)
      double precision,intent(inout) :: wsp(lwsp)
      integer,intent(in) :: nz,n,ia(nz),ja(nz)
      double precision,intent(in) :: a(nz)
      double precision, intent(out) :: res(n*n)
      
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )
      intrinsic ABS
      integer i,j
      
      do i = 1,n
         wsp(i) = ZERO
      enddo
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      enddo
      anorm = wsp(1)
      do i = 2,n
         if ( anorm.lt.wsp(i) ) anorm =  wsp(i)
      enddo
      do i = 1,n
         do j = 1,n
            v(j) = ZERO
         enddo
         v(i) = ONE
         call myDMEXPV( n, m, t,v,w, tol, anorm,wsp,lwsp,
     .             iwsp,liwsp, itrace,iflag,ia,ja,a,nz )
         do j = 1,n
            res(((i-1)*n)+j) = w(j)
         enddo
      enddo
      end subroutine wrapalldmexpv
      
      subroutine wrapsingledmexpv(n,m,t,v,w,tol,anorm,wsp,lwsp,iwsp,
     .                       liwsp, itrace, iflag,ia,ja,a,nz,res )
      implicit none

      integer,intent(inout) :: m, lwsp, liwsp, itrace, iflag
      integer,intent(inout) :: iwsp(liwsp)
      double precision,intent(inout) :: t, tol, anorm, v(n), w(n)
      double precision,intent(inout) :: wsp(lwsp)
      integer,intent(in) :: nz,n,ia(nz),ja(nz)
      double precision,intent(in) :: a(nz)
      double precision, intent(out) :: res(n)
      
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )
      intrinsic ABS
      integer i,j
      
      do i = 1,n
         wsp(i) = ZERO
      enddo
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      enddo
      anorm = wsp(1)
      do i = 2,n
         if ( anorm.lt.wsp(i) ) anorm =  wsp(i)
      enddo
      call myDMEXPV( n, m, t,v,w, tol, anorm,wsp,lwsp,
     .         iwsp,liwsp, itrace,iflag,ia,ja,a,nz )
      do j = 1,n
         res(j) = w(j)
      enddo
      end subroutine wrapsingledmexpv

      subroutine wrapdgpadm(ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,
     .                        ns,iflag )
      implicit none

      integer,intent(inout) :: ideg,m,ldh,lwsp,iexph,ns,iflag,ipiv(m)
      double precision,intent(inout) :: t,H(ldh,m),wsp(lwsp)
      integer i,j
*      print 9001,( (H(i,j), j=1,ldh), i=1,m )
 9001 format( 4(1X,D11.2) )
 9000 format( 4(1X,I2.1) )
      call DGPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)
*      print 9001,( (wsp(iexph+(j-1)*m+i-1), j=1,m), i=1,m )
*      print 9000,( ((iexph+(j-1)*m+i-1), j=1,m), i=1,m )
      end subroutine wrapdgpadm     


