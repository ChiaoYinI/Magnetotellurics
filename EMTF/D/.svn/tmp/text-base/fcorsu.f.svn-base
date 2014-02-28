c___________________________________________________
c
      subroutine fcorsu(dr,nfil,afparam,iftype,sc,rnrmt,cda,lfd)
        
      include 'iounits.inc'
      include 'decimate.inc'

c       computes table of transfer functions for subroutine rnorm
c       the transfer functions are used for correcting for filter
c       responses;  six things are corrected for
c       1) effect of first differencing
c       2) digital low pass filtering prior to decimation
c       3) correction of units of fourier coefficients (to counts/(sqrt(hz)))
c       4) analogue low pass (anti-alias) filters of instruments
c       5) analogue high pass filters of instruments (if any)
c        5') system response table
c       6) scale factors for individual channels to convert form counts
c                to physically meaningful units; for magnetic channels
c                this will be nt for electric mv/km; transformation of
c                measured channels to a "nice" coordinate system is not
c                done with this 

c       =========> fc(nd,nfcmx) = filter coefficients for digital filtering
c                       prior to decimation
c       =========> nfc(nd) = ith element is number of filter coffeficients 
c                       for decimating to level i
c       =========> npts(nd) = array of window lengths
c       =========> dr(nd) = array of sampling rates 
c       =========> afparam(nfilmax,nch) analogue filter/system response parameters
c>>>>>>>>>> 28 Feb, 1991 changed to character*80 array
c       =========> iftype(nfilmax,nch) indicator variable for formula for 
c                    filter/response
c       =========> sc(nch) = scale factors to convert counts to physical
c                       units
c       =========> cda clock offset in seconds (+ for fast)
c      
c       =========> lfd  logical array; true if corresponding
c                    channel/decimation level is to be first
c                     differenced
c       <========  rnrmt(nch,nd,npts/2) = conmplex table output
    
      parameter (pi=3.141592654)
      
      real dr(nd),sc(nch) 
      complex rnrmt(nch,nd,*),afcor,temp,tc
      integer iftype(nfilmax,nch),nfil(nch),nzeros,npoles,kz,kp
      logical lfd(nch,nd)
      real dt,t
      character*256 afparam(nfilmax,nch)

      integer id,i,npts2
      complex poles(100),zeros(100),w
      real*4 a0, dum
 
ccc   corrections for first difference if appropriate
ccc   correct units of fc's
      do id=1,nd
         npts2=nwin(id)/2
         pi2n=2.*pi/nwin(id)
         t = sqrt(dr(id))
         do ich = 1,nch
            w = cmplx(1.,0.)
            do i = 1,npts2
               if(lfd(ich,id))then
                  w=cmplx(0.0,pi2n*i)
                  w=1./(exp(w)-(1.0,0.0))
               end if
               rnrmt(ich,id,i)=w*sc(ich)*t
            enddo
         enddo
      enddo
        
ccc   corrections for low pass decimation filters
ccc   New version: 2-23-98   corrects for effect of low pass
ccc   filters on all decimation levels
      do id = 2,nd
         do jd = id,nd
            npts2 = nwin(jd)/2
            pi2n =2*pi*dr(id-1)/(dr(jd)*nwin(jd))
            do i = 1,npts2
               g = fc(id,1)
               do j=2,nfc(id)
                  t = pi2n*i*(j-1)
                  g=g+2.*fc(id,j)*cos(t)
               enddo
               do ich=1,nch
                  rnrmt(ich,jd,i) = rnrmt(ich,jd,i)/g
               enddo
            enddo
         enddo
      enddo
      
ccc  possible correction of fixed clock offset (seldom done here
ccc  anymore
      do id=1,nd
         do i = 1,npts2
            freq = i*(1./(dr(id)*nwin(id)))
            tc = cexp(cmplx(0.,freq*2.*pi*cda))
            do ich=1,nch
               rnrmt(ich,id,i)=tc*rnrmt(ich,id,i)
            enddo
         enddo
      enddo

ccc   corrections for analogue instrument filters/system response
      do ich = 1,nch
         do k = 1,nfil(ich)

c         read info for each filter out of character array, set
c               up tables for table look up of system response
            j = iftype(k,ich)
            if((j.eq.3).or.(j.eq.4).or.(j.eq.6)) then
               read(afparam(k,ich),*) gain,t0,alpha
            else if((j.eq.2).or.(j.eq.5)) then
               read(afparam(k,ich),*) gain,t0
               
ccc   EMI MT-1, BF-4 coils or EFSC
            else if(j.eq.7) then
               gain = 1.0
               read(afparam(k,ich),*) cfrsp,imode,igain
               cfrsp = 'sensors/'//cfrsp
               call resptbl(j,imode,igain,freq,temp,0)
            else if(j.eq.8) then
               gain = 1.0
               read(afparam(k,ich),*) cfrsp
               cfrsp = 'sensors/'//cfrsp
               call resptbl(j,imode,igain,freq,temp,0)
            else if(j.eq.10) then
ccc            POLES AND ZEROS representation of TF
               read(afparam(k,ich),*) nzeros,npoles,a0
               read(afparam(k,ich),*) dum,dum,dum, 
     &            (zeros(kz),kz=1,nzeros),(poles(kp),kp=1,npoles)

C-PDAS
C       31/7/96 read individual calibration files for each PAB, MFS05
            else if ((j.eq.17) .or. (j.eq.18)) then
               read (afparam(k,ich),*) cfrsp
               cfrsp = 'sensors/'//cfrsp
               call resptbl(j,imode,igain,freq,temp,0)
C     -PDAS
            elseif (j.eq.98 .or. j.eq.99) then
ccc            generic table
               read(afparam(k,ich),*) cfrsp
               cfrsp = 'sensors/'//cfrsp(1:iclong(cfrsp,80))
               call resptbl(j,imode,igain,freq,temp,0)
            elseif (j.eq.101) then
cBooker:       filter to correct timing offsets in individual channels
ccc            read time shift (assumed to be in seconds)
               read(afparam(k,ich),*) dt
cBooker:       filter to correct timing offsets in individual channels
            end if


ccc         this inner loop computes the actual filter response as a
ccc            function of frequency
            do id = 1,nd
            npts2 = nwin(id)/2
               do i = 1,npts2
                  freq = i*(1./(dr(id)*nwin(id)))
                  period = 1./freq
c               corect for filters, system response

C-PDAS
c       if it is a PDAS system iftype = j = 11,12
c       all parameters of the two filter types PD and MX are in 
c       function AFCOR, there is no need to read these from the char string 
c       afparam
                  if ((j.eq.11).or.(j.eq.12)) then
                     temp = afcor(j,t0,period,alpha)
                     gain = 1.0
                  end if
C-PDAS end
C LMT
c                 if it is an LMT system iftype = j = 13,14,15 or 16
c                 all parameters of the two filters are in FUNCTION AFCOR
                  if (((j.eq.13).or.(j.eq.14)).or.
     &                 ((j.eq.15).or.(j.eq.16))) then
                     temp = afcor(j,t0,period,alpha)
                     gain = 1.0
                  end if

                  if(j.le.6) then
                     temp = afcor(j,t0,period,alpha)
                  else if((j.eq.7).or.(j.eq.8)) then
                     call resptbl(j,imode,igain,freq,temp,1)
                  else if(j.eq.10) then
ccc                  Poles and Zeros
                     call pz_rsp(nzeros,npoles,a0,zeros,poles,freq,
     &                     temp)
                     gain = 1
                     
C                -PDAS  new calibration data files
                  else if ((j.eq.17) .or. (j.eq.18)) then
                     call resptbl(j,imode,igain,freq,temp,1)
                     gain = 1.0
C-PDAS

cvan    iftype = j = 97 is uncallibrated coil (i omega response)
                  else if (j.eq.97) then
                      temp = cmplx(0.0,freq)
                      gain = 1 

cnew    iftype = j = 98 is a generic amplitude/phase table
cnew    iftype = j = 99 is a generic real/imaginary table
                     elseif (j.eq.98 .or. j.eq.99) then
                     call resptbl(j,imode,igain,freq,temp,1)
                     gain = 1.0

		  elseif (j.eq. 101) then
cBooker:             filter to correct timing offsets in individual channels
                     gain = 1.
                     temp = cmplx(0.0,-2*pi*freq*dt)
                     temp = cexp(temp)

                  else if(j.eq.9) then
                     temp = 1.0
                     gain = 1.0
                  end if
                  temp = temp*gain
                  rnrmt(ich,id,i) = rnrmt(ich,id,i)/temp
               enddo
            enddo
         enddo
      enddo
 
      return
      end
