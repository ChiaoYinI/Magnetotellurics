c
c*****************************************
c
      subroutine edfwts(x,nch,nf,p1,p2,h,w)

      complex h(3),x(nch,nf)
      real edf,p1,p2,w(nf)

         do 10 i = 1,nf
         edf = x(1,i)*conjg(x(1,i))*h(1) + x(2,i)*conjg(x(2,i))*h(3)
     &             +2.*real(conjg(x(2,i))*x(1,i)*h(2))

         if(edf.gt.p2) then
            w(i) = 0.
         else if(edf.gt.p1) then
            w(i) = sqrt(p1/edf)
         else
            w(i) = 1.0
         end if
10       continue

      return
      end
c
c*****************************************
c
      subroutine edfwtsRR(x,nch,iref,nf,p1,p2,p3,h,w)

ccc   like edfwts, but weights depend also on similarity
ccc   of signal power between local and remote H sites
ccc   p3 = max deviation of ratio of amplitudes for reference channel
      complex h(3),x(nch,nf)
      real edf,edfref,p1,p2,w(nf)
      integer iref

      i1 = iref
      i2 = iref+1
      do i = 1,nf
         edf = x(1,i)*conjg(x(1,i))*h(1) + x(2,i)*conjg(x(2,i))*h(3)
     &             +2.*real(conjg(x(2,i))*x(1,i)*h(2))
         edfref = x(i1,i)*conjg(x(i1,i))*h(1) +
     &            x(i2,i)*conjg(x(i2,i))*h(3)
     &             +2.*real(conjg(x(i2,i))*x(i1,i)*h(2))

         if(edf.gt.p2) then
            w(i) = 0.
         else if(edf.gt.p1) then
            w(i) = sqrt(p1/edf)
         else
            w(i) = 1.0
         end if

         if(edfref.gt.p2) then
            w(i) = 0.
         else if(edf.gt.p1) then
            w(i) = w(i)*sqrt(p1/edf)
         end if

         if(edfref/edf .gt. p3 .or. edf/edfref .gt. p3 ) w(i) = 0.
       
      enddo

      return
      end
