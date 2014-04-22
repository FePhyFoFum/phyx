*---  Utility timer routine (seconds).
*---  uncomment the appropriate one and comment the others 

*---  SUN & SGI --
      double precision function clock()
      real*4 etime, tm(2)
      clock = etime( tm )
      end

*---  IBM & CRAY --
*      double precision function clock()
*      real*8 timef
*      clock = 1000.0d0 * timef()
*      end

*---  others ??


*---  if in trouble, use this to get out of the hook!
*      double precision function clock()
*      clock = 0.0d0
*      end
