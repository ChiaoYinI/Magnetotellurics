      logical inan
      integer k
      real w
      w=sqrt(k-1.5)
      write(*,*)w
      write(*,*)inan(w)
      stop
      end
c Lana Erofeeva, 2004
      logical function inan(w)
      real w
      integer k 
      character*80 c 
      inan=.false.
      write(c,*)w
      do k=1,80
       if(c(k:k).ne.' ')go to 1
      enddo
      return ! blanks only
1     c=c(k:80)
      write(*,*)'c(1:3):',c(1:3)
      if(c(1:3).eq.'NaN')inan=.true.
      return
      end
