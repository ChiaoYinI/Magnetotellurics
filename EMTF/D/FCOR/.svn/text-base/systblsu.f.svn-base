c___________________________________________________
c
        subroutine systblsu(dr,nch,nfil,afparam,iftype,sc,cda,
     *                      cfsTF,sensdir)
        
c       computes table of transfer functions for subroutine fctf
c       for uncorrecting transfer functions back to measurement units

ccc         corrections include
c       3) correction of units of fourier coefficients (to counts/(sqrt(hz)))
c       4) analogue low pass (anti-alias) filters of instruments
c       5) analogue high pass filters of instruments (if any)
c        5') system response table
c       6) scale factors for individual channels to convert form counts
c                to physically meaningful units; for magnetic channels
c                this will be nt for electric mv/km; transformation of
c                measured channels to a "nice" coordinate system is not
c                done with this 

c       =========> afparam(nfilmax,nch) analogue filter/system response parameters
c>>>>>>>>>> 28 Feb, 1991 changed to character*80 array
c       =========> iftype(nfilmax,nch) indicator variable for formula for 
c                    filter/response
c       =========> sc(nch) = scale factors to convert counts to physical
c                       units
c       =========> cda clock offset in seconds (+ for fast)
c      
c       <========  systbl(nch,nsystbl) = conmplex table output
      
      include 'fcor.inc'
      parameter (pi=3.141592654)

      real dr,sc(nch),freq(0:nsystbl),a0
      complex systbl(nchmx,nsystbl),afcor,temp,tc,zeros(100),poles(100)
      integer iftype(nfilmax,nch),nfil(nch),nzeros,npoles,dum
      character*256 afparam(nfilmax,nch)
      character*40 cfile
      character*80 cfsTF,cfrsp,sensdir
      integer nsens

      do k=1,80
       if(sensdir(k:k).eq.' ')then
         nsens=k-1
         go to 100
       endif
      enddo
100   continue
      if(nch.gt.nchmx) then
         write(0,*) 'ERROR IN systblsu !!!!'
         write(0,*) 'NCH = ',nch,'  NCHMX = ',nchmx
         write(0,*) 'STOPPING'
         stop
      endif
c>>>> correct units of fc's & for fixed clock offset
      freq(0) = steplog/(2.*dr)
      do j = 1,nsystbl
         freq(j) = freq(j-1)/steplog
c         tc = cexp(cmplx(0.,freq(j)*2.*pi*cda))
         do i = 1,nch
CCCCC   THIS (mult by sqrt dr) SHOULD NOT BE HERE !!!!! 
c            systbl(i,j) = sc(i)*sqrt(dr)
            systbl(i,j) = sc(i)
         enddo
      enddo
CCCCC   THIS SHOULD NOT BE HERE !!!!! 

c>>>>>   corrections for analogue instrument filters/system response,
      do 40 ich = 1,nch
         do 35 k = 1,nfil(ich)
c         read info for each filter out of character array, set
c               up tables for table look up of system response
         j = iftype(k,ich)
         write(0,*) 'ich,k,j = ',ich,k,j
         write(0,*) 'afparam = ',afparam(k,ich)
         if((j.eq.3).or.(j.eq.4).or.(j.eq.6)) then
            read(afparam(k,ich),*) gain,t0,alpha
         else if((j.eq.2).or.(j.eq.5)) then
            read(afparam(k,ich),*) gain,t0
         else if(j.eq.7) then
            gain = 1.0
            read(afparam(k,ich),*) cfile,imode,igain
            cfrsp = sensdir(1:nsens)//'/'//cfile
            call resptbl1(cfrsp,j,imode,igain,freq,temp,0)
         else if(j.eq.8) then
            gain = 1.0
            read(afparam(k,ich),*) cfile
            cfrsp = sensdir(1:nsens)//'/'//cfile
            call resptbl1(cfrsp,j,imode,igain,freq,temp,0)
         else if(j.eq.10) then
ccc         poles and zeros
            read(afparam(k,ich),*) nzeros,npoles,a0
            read(afparam(k,ich),*) dum,dum,dum,
     &            (zeros(kz),kz=1,nzeros),(poles(kp),kp=1,npoles)
cBooker:     filter to correct timing offsets in individual channels
         elseif (j.eq.101) then
ccc         read time shift (assumed to be in seconds)
            read(afparam(k,ich),*) dt
         end if

         do i = 1,nsystbl
            period = 1./freq(i)
c        corect for filters, system response
            if(j.le.6) then
               temp = afcor(j,t0,period,alpha)
            else if((j.eq.7).or.(j.eq.8)) then
               cfrsp = sensdir(1:nsens)//'/'//cfile
               call resptbl1(cfrsp,j,imode,igain,freq(i),temp,1)
            else if(j.eq.9) then
               temp = 1.0
               gain = 1.0
            else if(j.eq.10) then
               call pz_rsp(nzeros,npoles,a0,zeros,poles,freq(i),temp)
               gain = 1
            elseif (j.eq. 101) then
cBooker:       filter to correct timing offsets in individual channels
cEgbert: 	apparent phase error in sign assumed by Booker
               gain = 1.
               temp = cmplx(0.0,-2*pi*freq(i)*dt)
               temp = cexp(temp)
            end if
            temp = temp*gain
            systbl(ich,i) = systbl(ich,i)/temp
         enddo ! do i = 1,nsystbl

35       continue 
40    continue
                       
       open(unit=93,file=cfsTF)
       write(93,'(2i5)')  nch,nsystbl
       do i = 1,nsystbl
          write(93,'(e12.4)') freq(i)
          write(93,'(10e12.4)') (systbl(j,i),j=1,nch)
       enddo
      close(93)
      return
      end
