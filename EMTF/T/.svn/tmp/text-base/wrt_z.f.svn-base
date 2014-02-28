      subroutine wrt_z(z,nbt,fb,rdf,ldf,stname,orient,
     +   nchstk,ns,lrr,nsmx,zt,ibandlim,cdirout,level,
     +   stcor,decl,chid,samprate,nd,chead,ntmx,l_SDM,bw,idf)

cnew	append a DOS-style extension to the output file name
cnew	.zss for single site processing (lrref = .false.)
cnew	.zrr for remote reference (lrref = .true.)



      complex z(nsmx,*),zt(nsmx,*)
      real fb(nbt),rdf(ntmx,nbt),orient(2,*)
      real stcor(2,*),decl(*),samprate(nd),dw
      parameter (vcorps1  = .17)
      integer ldf(nbt),ibandlim(2,*),level(nbt),idf(nbt)
      logical lrr
      character*6 chid(*)
      character*40 stname
      character*40 cdirout
      character*80 cfile,chead
      character*4 c_ext
      real cfac(100)
      integer NJ,NI

      logical l_SDM

      if(lrr) then
         nch = nchstk-2
         c_ext = '.zrr'
      else
         nch = nchstk
         c_ext = '.zss'
      endif

      ll = iclong(stname,40)
      if (stname(ll-4:ll) .eq. '.zss' .or.
     &     stname(ll-4:ll) .eq. '.zrr') then
         c_ext = ''
      endif 

      m = iclong(cdirout,40)
      if (m.gt.0) then
         cfile = cdirout(1:m)//'/'//stname(1:ll)//c_ext
      else
         cfile = stname
      endif
      open(unit=3,file=cfile)

cnew  	if requested by the -S command line option write the spectral density
cnew	matrix to a file (extension .sdm)
      if (l_SDM) call sdmwrite(nbt,ns,stname,cdirout,
     &   nsmx,z,fb,rdf,nchstk,chid,orient,bw,ntmx,stcor,decl)
     
cnew      print*,'nbt ',nbt,' ns ',ns,' nsmx ',nsmx
      do i = 1,nbt
         do j = 1,ns
            zt(j,i) = z(j,i)
         enddo
      enddo
	  

c...  write header information 
      write(3,*) '**** IMPEDANCE IN MEASUREMENT COORDINATES ****'
      write(3,*) '********** WITH FULL ERROR COVARINCE**********'
      write(3,'(a80)') chead
      write(3,100) stname
      write(3,105) stcor(1,1),stcor(2,1),decl(1)
      write(3,110) nch,nbt
      write(3,*) 'orientations and tilts of each channel '
      do k=1,nch
         write(3,115) k,orient(1,k),orient(2,k),stname(1:3),chid(k)
      enddo
      write(3,*)

      NJ = 2
      NI = nch - NJ
      do ib = 1,nbt
         if (lrr) then
            call trlsrr(zt(1,ib),nch,ldf(ib))
         else
            call tranls(zt(1,ib),nch,ldf(ib))
         endif
      enddo

      do ib = 1,nbt
c...     correcting SR
         do i = 1,NI
            k = ((i+2)*(i+3))/2
            zt(k,ib) = cmplx(real(zt(k,ib)),0.0)
ccc         OLD correction for E(psi') term in RME asymptotic errors
ccc            cfac(i) = (2.*ldf(ib))/rdf(i,ib)
ccc         NEW: extra correction for redescending influence curve:
ccc         dw is fraction of downweighted df
            dw = (2*ldf(ib) - rdf(i,ib))/(2*ldf(ib))
ccc		write(0,*) ldf(ib), ib, rdf(i,ib),dw
ccc         OLD should be equivalent to
ccc         cfac(i) =  1/(1-dw)
CCC         NEW CORRECTION ALLOWS FOR NEGATIVE PART OF psi' curve
ccc              (approximately)
            cfac(i) = 1/(1-2*dw)
ccc         test Apr14
	    cfac(i) = 4*ldf(ib)*ldf(ib)/(idf(ib)*2*(2*rdf(i,ib)-2*idf(ib)))
         enddo
         do i = 1,NI
            do j = 1,i
               k = ((i+2)*(i+1))/2 + 2 + j
               zt(k,ib) = zt(k,ib)*cfac(i)*cfac(j)
            enddo
         enddo

c...     first, correcting SS for correlation of adjacent freq.
         cor = ibandlim(2,ib)-ibandlim(1,ib) + 1.
         cor = 2.*vcorps1*(cor-1.)/cor
         do k = 1,3
            zt(k,ib) = zt(k,ib)*(1.+cor)
         enddo

         write(3,120) fb(ib),level(ib),ibandlim(1,ib),ibandlim(2,ib)
         write(3,125) ldf(ib),1./samprate(level(ib))
c...     TF
         write(3,*) 'Transfer Functions'
         do i = 1,NI
            kk = ((i+2)*(i+1))/2 + 1
            write(3,140) (zt(k,ib),k=kk,kk+1)
         enddo
c...     SIGMA S (inverse coherent signal power matrix)
         write(3,*) 'Inverse Coherent Signal Power Matrix'
         write(3,140) (zt(k,ib),k=1,1)
         write(3,140) (zt(k,ib),k=2,3)
c...     SIGMA R (residual covariance matrix)
         write(3,*) 'Residual Covariance'
         do i = 1,NI
            kk = ((i+2)*(i+3))/2 + 1 - i
            write(3,140) (zt(k,ib),k=kk,kk+i-1)
         enddo
      enddo ! ib
      close(3)
      
      print*
      print*,' Output written to ',cfile(1:iclong(cfile,80))
      print*

100   format('station    :',a20)
105   format('coordinate ',f9.3,1x,f9.3,1x,' declination ',f8.2)
110   format('number of channels ',i3,2x,' number of frequencies ',i4)
115   format(i5,1x,f8.2,1x,f8.2,1x,a3,2x,a6)
120   format('period : ',f12.5,3x,' decimation level ',i3,3x,
     +       ' freq. band from ',i4,' to ',i4)
125   format('number of data point ',i6,' sampling freq. ',f7.3,' Hz')
140   format(16e12.4)
      return
      end ! of wrt_z.f
