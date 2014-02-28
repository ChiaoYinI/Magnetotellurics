ccc   NEW VERSION .... August, 1994  (modified repeatedly through Feb. 1996)
ccc   souped up, but fairly general version of multiple station EM
ccc   data analysis program.  Robust in several ways: robust covariance;
ccc   automatic noise scale estimates; timing errors fixed;
ccc   uses robust regression on principal components to "fix" local outliers
ccc   in all channels; can analyze extra sub-arrays (in principal -
ccc   this probably needs more testing); can be set to try separating plane-wave
ccc   from gradient; can be set to output residuals; can output principal
ccc   component time series, array TFs etc., etc., etc.
ccc   Does multiple remote apparent resistivity and phases ... 
c 
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>    parameters to set   <<<<<<<<<<<<<<<<<<<<<

c       array size parameters are set in iosoze
      include 'iosize.inc'
ccc   robust itaration maxima : defaults are set here, have _d appended
ccc     to variable names
ccc   itmax_sdm = maximum # of iterations used for robust estiamte of SDM
ccc                (using rotationaly invariant Huber estimate)
ccc   itmax_ln = maximim # of iterations used for robust regression for
ccc           each inner loop estimate of local noise
ccc   itmax_cln = maximum # of outer loops for local noise estimation/
ccc               individual channel outlier cleanup
ccc   itmax_rrr = maximum # of iterations used for robust regression for
ccc                final TF estimate (inner loop ... if coherent noise
ccc                    downweighting is used)
ccc   itmax_cn = maximum # of iterations used for coherent noise
ccc           downweighting in outer loop for final TF estimate
ccc             LEAVE THIS DEFAULT SET TO 0 !!!

      parameter (itmax_sdm_d = 40,
     &           itmax_ln_d=20,
     &           itmax_cln_d=3,
     &           itmax_rrr_d = 40,
     &           itmax_cn_d=0,
     &           nusetMax=50)

ccc   r0_d sets default value for outlier cutoff in Huber estimate
ccc   for robust regression with complex data

      parameter ( r0_d = 1.4)
cc

c=========================================================================

      parameter (nxmx = nchmx*nstamx*ntfmax,niwrk=nstamx*nsetmx,
     &      ncwrk=nstamx*nchmx*nsetmx,npmx=ntmx,numx=2*ntmx*npmx,
     &        ns1mx=(nchmx*(nchmx+1)/2),nsigmx=5,nxpcmx=nsigmx*ntfmax,
     &        p_d = .1,ntmx1 = ntmx+1)
      integer ibandlim(2,nbmax),isuse(2,ndmax,nusetMax),iband(2),ilxt,
     & ilx,idl(nbmax),ih_tot(ntmx1),nch(nstamx),ntape(nstamx),
     &iterrb(3,2), nd(nstamx),npts(ndmax,nstamx),na(1000),iwrk(niwrk),
     &iux(nstamx),ih(ntmx),ivunits(10),iotfn,ldf(ngrpmx,nbmax),
     & nt(50),nf(50),npoint(50),i1(ntmx),isuse_a(2,ndmax,nusetMax),
     & ie_tot(ntmx), n_grp,ih_grp(ntmx1),nch_grp(ntmx1),ioux,ioupw(2)
     & ,ngrp_tf,igrp_tf(ntmx+1,ngrpmx),sta_grp(nchemx,ngrpmx),id_list,
     &   nsets(nstamx),lar_id,iclong,brUnit

      integer nuset, dec_fac, nuset_a

      real w1(2,nsmx),dr(ndmax,nstamx),ev(ntmx),orient(2,ntmx),
     &  decl(nstamx),period(nbmax),decl_tot(nstamx),machep,sigpwr(2),
     &  stc_tot(2,nstamx),wrk(ntmx,ntmx,10),u(numx),tpower(nbmax),
     &  var(ntmx),rho(2,nstamx,nbmax),vartmp(ntmx),sndiag(ntmx,nbmax),
     &  eval(ntmx,nbmax),eval2(ntmx,nbmax),sig_te(nstamx),
     &  phi(2,nstamx,nbmax),sdiag(ntmx,nbmax),te(nstamx),noisepwr(2),
     &  stkm(2,nstamx),stc(2,nstamx),u2(numx),ev2(ntmx),cc(ntmx,nbmax)
     & ,ug(numx),u2g(numx),wt(ntfmax),dr1(nstamx),pwrmin
      complex cpwg(6,20),v(ntmx,5),xpc(nxpcmx),zt(nchemx,2)
      real pwg(11,20),dpwg(ntmx),ss(6),damp(5),varpwg(6),uvar(ntmx),
     &   z_err(2,3,nstamx,nbmax),rhoa(12,nstamx,nbmax)
     &   ,rhoa2(12,nstamx,nbmax)

      logical lprint,lplt,lex,lpw,lstcc,lx(nstamx,ntfmax)
     &   ,lterr,lch(ntmx),luse(ntfmax),la(nstamx,50),lgrad,
     &    ltfn,lres,lomit(ntfmax),lpwout,luse_all,lraw_out,
     &    l_CC,l_CM,l_PC,l_Pw,l_ref_proj,l_GRP,lrot,l_set_list,
     &   e_and_h,lBadRec,l_NoDeadChannels
      complex s(nsmx),xr(nxmx),xc(nxmx),xp(nxmx),
     &   z(nchemx,2,nstamx,nbmax),
     &   sig_s(2,2,nstamx,nbmax),sig_e(nchemx,nchemx,nstamx,nbmax)
     &   ,y(ntmx2),a(2),c(2),b(ntmx2),cwrk(ncwrk),w(ntmx,7)
     &   , y2(ntmx),s2(nsmx),pw(ntmx2),cnoise(ntmx,3)

      character*40 stnames(nstamx,ntpmax)
      character*20 cgrp(ngrpmx),arrayid
      character*10 cevec,cvec
      character*4 sta(nstamx),stat(nstamx),cdum
      character*2 cb(nbmax)
      character*40 cbandf
      character*82 arg
      character*80 comment,cfile,chead,badRecFile
      character*60 cf_array,cfile_grp
      character*6 chid(ntmx),chid_grp(2+nchemx,ngrpmx)
      character*1 cjob,cjobwt,option_grp
      logical inan

c============> logical variables which control program options
      do ista = 1,nstamx
         ljunk(ista) = .false.
      enddo
      l_NoDeadChannels = .true.
      lex = .false.
      lBadRec = .false.
c       lex is true to process sub-arrays
ccc     if lex = .true. need to declare logical array lx(nstamx,ntfmx)
ccc              (See above)

      lpw = .true.
c        lpw is true to transform two leading evecs so that average
c          horizontal fields are uniform linearly polarized

      lgrad = .false.
ccc       lgrad is .true. if source gradients should be fit for plane-wave
ccc         source estimates
ccc    (these last two are more or less mutually exclusive; change with
ccc      commanda line option -g)


      lterr = .false.
c         lterr is true to try fixing timing errors; program looks
c         for a file named *.TER , where * is the "array id"
c         specified in the cf_array file;  if this file is found,
c         time offsets for each station in seconds are read in and
c         used for the correction.  If the file is not found and
c         if lterr = .true., the program estimates timing errors by looking
c         at the relative phase variation of horizontal mag fields at high
c         frequencies;  (might want to adjust the frequency band parameters
c           used for this purpose .... see nterrb, iterrb(.,.) below).
c         When the program estimates timing errors, the results are output
c         in file *.TER
c         

      ltfn = .false.
ccc     ltfn is .true. to output TFs used for local noise estimation
ccc     (change default with command line option)

      luse_all = .false.
ccc      luse_all is true to omit data set ranges specified in
ccc      cf_array file only on first coherent noise downweighting
ccc      iteration; on subsequent iterations coherent noise/signal
ccc      ratios are used to decide which data should be omitted
ccc      (change default with command line option)

      lraw_out = .false.
ccc      lraw_out is true to output all components of raw data in PC_***** file
ccc   use command line option -r to change
      l_PC = .false.
ccc   use command line option -P to change
      l_CC = .false.
ccc   use command line option -c to change
      l_CM = .false.
ccc   use command line option -C to change
      l_Pw = .true.
      ioupw(1) = 55
      ioupw(2) = -1
ccc   use command line option -L to change
      l_ref_proj = .false.
ccc   use command line option  -s to change
      l_GRP = .false.
ccc   lrot is true to rotate magnetic components into fixed geomag/geographic
ccc   coordinate system (depends on what system the channel orientations
ccc   are given in)
      lrot = .false.
ccc   l_set_list is true to (only) make a list of available "set numbers"
ccc    for a fixed decimation level
      l_set_list = .false.

ccc   lres is .true. to output residuals (ouch)
      lres = .false.

ccc   option_grp tells what grouping of data channels to use for 
ccc   computing incoherent noise variance estimates
ccc   Use command line option -GX to change where X is T, S, or A
ccc   Options for defining groups are:
ccc   T ==> all components of a type at a single site (default)
ccc   S ==> all components at a single site
ccc   A ==> each component by itself
      option_grp = 'T'

      cjob = 'R'
      cjobwt = 'C'
      cvec =  'VECTOR #'
      cevec = 'SIG VEC #'
      ioux = 0

ccc   parse command line for options ...
      narg = iargc()
      r0 = r0_d
      itmax_rrr = itmax_rrr_d
      itmax_cn = itmax_cn_d
      itmax_ln = itmax_ln_d
      itmax_sdm = itmax_sdm_d  
      itmax_cln = itmax_cln_d
      p = p_d
      pwrmin = 0.0
      cf_array = 'array.cfg'
      do k = 1,narg
         call getarg(k,arg)

         if(arg(1:2).eq.'-n') then
ccc         turn off robust features:
            r0 = 100.
            itmax_rrr = 0
            itmax_ln = 0
            itmax_sdm = 0
            itmax_cln = 0
            itmax_cn = 0

         else if(arg(1:2) .eq. '-m') then
ccc         turn off all robust features EXCEPT SDM stack
            itmax_rrr = 0
            itmax_ln = 0
            itmax_cln = 0
            itmax_sdm = itmax_sdm_d

         else if(arg(1:2).eq.'-f') then
ccc         change default array file
            cf_array = arg(3:62)

         else if(arg(1:2).eq.'-T') then
ccc         turn on automatic timing error correction
            lterr = .true.

         else if(arg(1:2).eq.'-t') then
ccc         output TF files
            ltfn = .true.

         else if(arg(1:2).eq.'-u') then
ccc         assume signal and coherent noise are uncorrelated for
ccc         computing weights
            cjobwt = 'U'

         else if(arg(1:2).eq.'-N') then
ccc         dont transform eigenvectors in M_ file
            lpw = .false.

         else if(arg(1:2).eq.'-g') then
ccc         try to sort out plane wave/gradient sources geometrically
            lgrad = .true.
            lpw = .false.

         else if(arg(1:2).eq.'-c') then
            l_CC = .true.
            cjob = arg(3:3)
            if((cjob.ne.'H').and.(cjob.ne.'N')) then
ccc            change default output for canonical covariance
ccc            H = normalize E to nT using avg. impedance
ccc            R = canonical coherence (default)
ccc            N = noise units
               cjob = 'R'
            endif

         else if(arg(1:2).eq.'-i') then
ccc         change default max # of iterations used for robust regression 
ccc            for each inner loop estimate of local noise
            read(arg(3:60),*) itmax_ln

         else if(arg(1:2).eq.'-o') then
ccc         change default max # of outer loops for local noise estimation/
ccc               individual channel outlier cleanup
            read(arg(3:60),*) itmax_cln

         else if(arg(1:2).eq.'-I') then
ccc         change default max # of iterations used for robust regression 
ccc                for final TF estimate (inner loop ... if coherent noise
            read(arg(3:60),*) itmax_rrr

         else if(arg(1:2).eq.'-w') then
ccc         change default number of iterations for coherent noise
ccc         downweighting
            read(arg(3:60),*) itmax_cn

         else if(arg(1:2).eq.'-C') then
ccc         output CM file
            l_CM = .true.

         else if(arg(1:2).eq.'-L') then
ccc         don't output Pw file 
            l_Pw = .false.

         else if(arg(1:2).eq.'-P') then
ccc         output PC file
            l_PC = .true.

         else if(arg(1:2).eq.'-M') then
ccc         eliminate points with low estimated signal power from
ccc         fit; set cutoff level with the argument
            read(arg(3:60),*) pwrmin
            write(0,*) 'pwrmin = ',pwrmin

         else if(arg(1:2).eq.'-p') then
ccc         change default cutoff for coherent noise downweighting
            read(arg(3:60),*) p

         else if(arg(1:2).eq.'-a') then
ccc         omit data set ranges specified in
ccc         cf_array file only on first coherent noise downweighting
ccc         iteration; on subsequent iterations coherent noise/signal
ccc         ratios are used to decide which data should be omitted
            luse_all = .true.

         else if(arg(1:2).eq.'-r') then
            l_PC = .true.
            lraw_out = .true.

         else if(arg(1:2).eq.'-R') then
            l_ref_proj = .true.
            read(arg(3:62),*) irefsta

         else if(arg(1:2).eq.'-G') then
ccc         change channel grouping option for incoherent
ccc         noise computation
            option_grp = arg(3:3)

         else if(arg(1:2).eq.'--') then
            call usage(itmax_ln_d,itmax_cln_d,itmax_rrr_d)

         else if(arg(1:2).eq.'-S') then
ccc         just make a list of available sets and quite
ccc         use decimation level id_list for this ...
            read(arg(3:62),*) id_list
            l_set_list = .true.

         else if(arg(1:2).eq.'-s') then
            l_GRP = .true.
            cfile_grp = arg(3:62)

         else if(arg(1:2).eq.'-z') then
            lrot = .true.

         else if(arg(1:2).eq.'-b') then
            lBadRec = .true.
            read(arg(3:82),'(a80)') badRecFile

         endif
      enddo

c       get program control instructions from user:
c       open files for io, set up bands, make
c       directory array for direct access of frequency files
        
ccc   various initializations
      call mmt_init(machep,maxint,isuse,npts,ndmax,nstamx)
      nuset = 1

ccc     get info from array file (default is named cf_array)
       open(unit=1,file = cf_array,status='old')
ccc    read in number of stations from array file
       read(1,*) nsta_tot

ccc   get further information from array file, and then open output files,
ccc     do various set-up calculations etc.
      call setup(nsta_tot,ntape,stnames,arrayid,nd,nch,cbandf,
     &     stc_tot,decl_tot,orient,chid,dr,npts,var,isuse,nuset,sta,
     &     itmax_ln,nt_tot,ih_tot,ie_tot,ns,theta)
      write(0,*) 'multmtrn version of 29 Nov., 2006, second try'
      write(0,*) 'array ID = ',arrayid
      write(0,*) 'Done with setup'

      if(lBadRec) then
ccc      set isuse_a using info from bad record file
         brUnit = 88
         fac = .75
         call badToUse(badRecFile,brUnit,ndMax,nuset_a,isuse_a,
     &      npts,fac)
      else
         nuset_a = 0
      endif
      if (nuset_a .eq. 0) then
         nuset_a = nuset
         do i = 1,2
            do j = 1,ndmax
               isuse_a(i,j,1) = isuse(i,j,1)
            enddo
         enddo
      endif

      lar_id = iclong(arrayid,20)

      if(l_set_list) then
ccc      just make a list of available set numbers ...
ccc      first make directory for fourier coefficient files; open files
         call mkfdir(nsta_tot,ntape)
         if(.not.lfop) then
            call fop(nsta_tot,nch,ntape)
         end if
         call mk_list(nd,id_list,nsta_tot,ntape,nch,nsets,iwrk,cwrk)
         cfile = arrayid(1:lar_id)//'.SET'
         call wrt_sets(nsta_tot,id_list,sta,nsets,iwrk,cfile)
         stop
      endif

      if(l_GRP) then
ccc     read in channel grouping file to redefine chanel groupings
ccc     for TF estimation/z-file output
         open(unit=96,file=cfile_grp)
ccc      first line of non-standard grouping file is # of channel groups
         read(96,*) ngrp_tf
ccc      for each channel group ...
         do k = 1,ngrp_tf
ccc         first a name (for output file roots)
            read(96,'(a20)') cgrp(k)
            write(0,*) cgrp(k)
ccc         second the # of components in the group
            read(96,*) igrp_tf(1,k)
ccc         finally the list of component numbers 
            read(96,*) (igrp_tf(l,k),l=2,igrp_tf(1,k)+1)
            do l = 1,igrp_tf(1,k)
               do ista = 1,nsta_tot
                  if((igrp_tf(l+1,k).ge.ih_tot(ista)).and.
     &               (igrp_tf(l+1,k).lt.ih_tot(ista+1))) then
                     sta_grp(l,k) = ista
                     ll = igrp_tf(l+1,k)
                     chid_grp(l,k) = chid(ll)
                  endif
               enddo
            enddo
         enddo
         close(96)
      else
ccc      use default grouping determined from FC files
         ngrp_tf = nsta_tot
         ll = 0
         do k = 1,nsta_tot
            igrp_tf(1,k) = nch(k)
            cgrp(k) = sta(k)
            do l = 2,nch(k)+1
               ll = ll + 1
               sta_grp(l-1,k) = k
               igrp_tf(l,k) = ih_tot(k)+l-2
               chid_grp(l-1,k) = chid(ll)
            enddo
         enddo
      endif      

ccc   see if there is a "location fix file ... put this in to deal with
ccc   an error in the locations given in the FC files (there is in MSMT data)
      open(unit=88,file='STCOR',status='old',err=1)
         do ista = 1,nsta_tot
            read(88,*) stc_tot(1,ista),stc_tot(2,ista),decl_tot(ista)
         enddo
      close(88)
1     continue
      lprint = .true.
      nskip = 1
      
ccc.... set up for frequency band averaging
      call bset(nbt,ibandlim,idl,period,ndmax,nbmax,ierr,
     &     dr,npts,cbandf)

   
c     GE        Compute the sets to use for processing
c     (following (sort of) ME  modifications in tranmt)
         nduset = idl(nbt)
         do i = 1, nuset
            do k = 2, nduset
            dec_fac = nint(dr(k,1)/dr(k-1,1))
            isuse(1,k,i) = int(isuse(1,k-1,i)/dec_fac)
            isuse(2,k,i) = int(isuse(2,k-1,i)/dec_fac)
            if (mod(isuse(1,k-1,i),dec_fac).gt.0)
     &           isuse(1,k,i)=isuse(1,k,i)+1
            if (mod(isuse(2,k-1,i),dec_fac).gt.0)
     &           isuse(2,k,i)=isuse(2,k,i)+1
            enddo
         enddo
         do i = 1, nuset_a
            do k = 2, nduset
            dec_fac = nint(dr(k,1)/dr(k-1,1))
            isuse_a(1,k,i) = int(isuse_a(1,k-1,i)/dec_fac)
            isuse_a(2,k,i) = int(isuse_a(2,k-1,i)/dec_fac)
            if (mod(isuse_a(1,k-1,i),dec_fac).gt.0)
     &           isuse_a(1,k,i)=isuse_a(1,k,i)+1
            if (mod(isuse_a(2,k-1,i),dec_fac).gt.0)
     &           isuse_a(2,k,i)=isuse_a(2,k,i)+1
            enddo
         enddo
ccc   GE

      if(l_CM) then
ccc      file for correlation matrix
         iocm = 72
         open(unit=iocm,file=arrayid(1:lar_id)//'.CM',status='unknown')
      else
         iocm = 0
      endif
        
      if(l_PC) then
         ioupc = 52
         cfile = arrayid(1:lar_id)//'.PC'
         if(lraw_out) then
            call pc_outinit(ioupc,cfile,nsta_tot,nt_tot,nt_tot,nbt,
     &        stc_tot,decl_tot,orient,ih_tot)
         else
            call pc_outinit(ioupc,cfile,nsta_tot,nt_tot,nsigmx,nbt,
     &        stc_tot,decl_tot,orient,ih_tot)
         endif
      endif

      if(l_Pw) then
ccc      file for "Remote Reference Plane Wave TF"
         cfile = arrayid(1:lar_id)//'.Pw'
         open(unit=ioupw(1),file=cfile,status='unknown',
     &          form='unformatted')
         write(ioupw(1)) nt_tot,nsta_tot,2,nbt
         do ista = 1,nsta_tot
            write(ioupw(1)) nch(ista),ih_tot(ista),
     &      stc_tot(1,ista),stc_tot(2,ista),decl_tot(ista),sta(ista)
         enddo
         ll = 0
         do ista = 1,nsta_tot
            do l = 1,nch(ista)
               ll = ll + 1
               write(ioupw(1)) orient(1,ll),orient(2,ll),chid(ll),
     &            sta(ista)
            enddo
         enddo
      endif

ccc...make directory for fourier coefficient files; open files
      call mkfdir(nsta_tot,ntape)
      if(.not.lfop) then
         call fop(nsta_tot,nch,ntape)
      end if

ccc...setup output of signal/noise results (main + extra arrays)
ccc   for array combining program
      call init_xtr(nsta_tot,nt_tot,ntape,nch,nd,ih_tot,lx,isuse_a,
     &     nuset_a,xr,wrk,iwrk,cwrk,nxu,la,ivunits,arrayid,sta,ixs,ixf,
     &     stc_tot,chid,decl_tot,orient,nbt,nsigmx)

      if(.not.lex) nxu = 1      

ccc...files for TF outputs
      if(ltfn) then
ccc      output TFs used for local noise estimation
         iotfn = 73
         open(unit=iotfn,file=arrayid(1:lar_id)//'.TF',
     &        status='unknown')
         write(iotfn,'(3i4)') nt_tot,nt_tot,nbt
      else
         iotfn = 0
      endif

ccc.....setup complete; ready to start making SDMs
      if(lterr) then
         open(unit=57,file=arrayid(1:lar_id)//'.TER',
     &        status='old',err=950)
         do ista = 1,nsta_tot
            read(57,'(a4,2x,2e12.4)') cdum,te(ista),sig_te(ista)
            if(inan(te(ista)))te(ista)=0.
            if(inan(sig_te(ista)))sig_te(ista)=1e9
            if(abs(te(ista)).lt.(2.*sig_te(ista))) te(ista) = 0.
c            write(0,*)'OLD:',cdum,te(ista),sig_te(ista)
         enddo
         close(57)
         go to 960
950      continue
         nterrb = 1
ccc      OLD standard
         iterrb(1,1) = npts(1,1)/8
         iterrb(2,1) = (npts(1,1)/2)*.80
         iterrb(3,1) = 4
ccc      modified for NIMS 8 hz: go from band 16 (1 sec) --> 4 (4 sec), steps of 2
         iterrb(1,1) = 2
         iterrb(2,1) = 8
         iterrb(3,1) = 1

         twin = npts(1,1)*dr(1,1)
         call timerr(nsta_tot,nt_tot,ntape,nch,nd,orient,ih_tot,ie_tot,
     &      twin,xr,xc,lx,isuse_a,nuset,iterrb,nterrb,te,sig_te,u,wrk,
     &       cwrk,iwrk,s)
         open(unit=57,file=arrayid(1:lar_id)//'.TER',status='new')
         do ista=1,nsta_tot
            if(inan(te(ista)))te(ista)=0.
            if(inan(sig_te(ista)))sig_te(ista)=1e9
            if(abs(te(ista)).lt.(2.*sig_te(ista))) te(ista) = 0.
            write(57,'(a4,2x,2e12.4)') sta(ista),te(ista),sig_te(ista)
c            write(0,*)'NEW:',sta(ista),te(ista),sig_te(ista)
        enddo
         close(57)
      end if
960   continue

ccc   loop over NBT frequencies
      if(nbt .eq. 0 ) stop
      do ib=1,nbt
         iband(1)=ibandlim(1,ib)
         iband(2)=ibandlim(2,ib)
         id=idl(ib)
         ib1=iband(1)
         ib2=iband(2)
         idec = npts(id,1)*dr(id,1)/(npts(1,1)*dr(1,1))
ccc      garbage to screen to entertain user
         write(0,1010) ib,period(ib)
1010     format('______band #',i3,'; period=',g8.3,' s ______________',
     &       '_______________')
         write(0,1020) 1./dr(id,1),npts(id,1),id,iband
1020     format('Samp. Rate=  ',g8.3,'Hz; Window= ',i4,' Pts.;',
     &    ' Dec. Level= ',i1,'; FCs ',i3,'-',i3) 
ccc              end of garbage (for now)

ccc          get FCs for band ::: XX is data matrix (all data, including
ccc                 FCs with missing stations)
         call mkrec(nd,id,iband,nsta_tot,ntape,nch,nt_tot,isuse_a,
     &      nuset_a,nf_tot,xp,.true.,lx,ixs,ixf,cwrk,iwrk)

c         write(0,*) 'nf_tot = ',nf_tot
c         write(0,*) 'te',te

         if(lterr) then
ccc         adjust for timing errors
            call terrfix(xp,nt_tot,nf_tot,nsta_tot,ih_tot,period(ib),te)
         end if

ccc      Put data into format appropriate for rbstreg .... sort by
ccc      sub-arrays, generate pointers and array sizes for all sub-arrays,
ccc      reorder so that each channel is a column in the data matrix.
ccc      Result is put in array XR; XP can now be used for cleaning data
         call data_sort(xp,nf_tot,nt_tot,nsta_tot,ih_tot,lx,la,nxu,
     &        xr,nf,nt,npoint,luse,ixs,ixf)

c         write(0,*) 'nxu = ',nxu
c         write(0,*) 'nf,nt,npoint',(nf(k),nt(k),npoint(k),k=1,nxu)

ccc        loop over sub-arrays (array for all stations is first ...)
          do ix = 1,nxu
ccc          set up various component identifiers for sub-array ix
             call mask_ch(la,ix,lch,nsta_tot,nsta,ih_tot,ih,
     &            nt_tot,stc_tot,stc,decl_tot,decl)

ccc         for ix = 1 (full array do various coherence, noise calculations
             if(ix.eq.1) then
                if(lres) then
                   ioux = 59
                   cfile = 'RES_'//cb(ib)
                   call openx(ioux,cfile)
                   nfb = iband(2)-iband(1)+1
                   nseg = nf(1)/nfb
                   do i = 1,nt_tot
                      vartmp(i) = 1.0
                   enddo
                   comment = 'Raw data'
                   call wrtx(ioux,xr,nt_tot,nseg,nfb,vartmp,comment)
                endif
                call mkgrp(option_grp,ih_tot,ie_tot,nsta,nt(1),
     &                     ih_grp,n_grp,nch_grp)
                call n_rbst(xr,xc,xp,nf(1),nt(1),n_grp,ih_grp,nch_grp,
     &           sndiag(1,ib),s,period(ib),itmax_sdm,itmax_ln,
     &           itmax_cln,r0,ioux,iocm,iotfn)
ccc     Now S contains non-iterative SDM based on unweighted,
ccc     "cleaned" data.  Extract diagonal elements for comparison with
ccc     "robust" noise estimates.
                call diag(s,nt(1),sdiag(1,ib))
                do k = 1,nt(1)
                  if(sdiag(k,ib) .eq. 0) then
                     l_NoDeadChannels = .false.
                     write(0,*) 'Channel ',k ,' is dead'
                  endif
                enddo
ccc             now iteratively stack cleaned data
                call rbstk2(xc,nt(1),nf(1),s,xp,itmax_sdm)
                call diag(s,nt(1),sdiag(1,ib))
                if(l_CC) then
                   call lrcov(cjob,s,sndiag(1,ib),nt(1),nsta_tot,
     &               ih_tot,u,cc(1,ib),wrk)
                endif
             else
ccc             just form SDM
                call rbstk2(xr(npoint(ix)),nt(ix),nf(ix),s,xc,
     &            itmax_sdm)
             end if

ccc          find eigenvectors/values of weighted, cleaned SDM
             ii = 0
             do i=1,nt_tot
                if(lch(i)) then
                   ii = ii + 1
                   vartmp(ii) = sndiag(i,ib)
                endif
             enddo
             icon = 1
             call rsp(s,nt(ix),nsigmx,u,ev,icon,wrk,vartmp)
             nsig =  n_eig_sig(nf(ix),ev,nt(ix),nsigmx)
             if(lgrad) nsig = 5
             call filtpc(nt(ix),nsig,u,ev,u2,ev2,wrk,ih,nsta)

ccc          output leading evec, evals, variances used
             irec = ib+1
             ns = (nt(ix)*(nt(ix)+1))/2
ccc          now writing out S instead of leanding evecs,evals
             call wrt_s(ivunits(ix),irec,s,ns,vartmp,nt(ix),
     &         nf(ix),period(ib))
            
             if(ix.eq.1) then
ccc             output projection of FCs onto NSIG leading PCs
                nfb = iband(2)-iband(1)+1
                if(lraw_out) then
                   call pc_out(ioupc,xr,xpc,nt(1),nf(1),nfb,id,iband,
     &               ixs,u,ev,nsigmx,vartmp,period(ib),lraw_out)
                else if(l_PC) then
                   call pc_out(ioupc,xc,xpc,nt(1),nf(1),nfb,id,iband,
     &               ixs,u,ev,nsigmx,vartmp,period(ib),lraw_out)
                endif
ccc             and save eigenvalues ev and ev2
                tpower(ib) = 0.
                do i=1,nt(ix)
                   if(ev2(i).le.machep) ev2(i) = machep
                   tpower(ib) = tpower(ib) + ev2(i)
                   if(ev(i).le.machep) ev(i) = machep
                   eval(i,ib) = ev(i)
                   eval2(i,ib) = ev2(i)
                enddo
             endif

ccc          rotate columns of U, U2 into fixed geomagnetic coordinate system
ccc          first save into arrays UG, U2G
             nt2 = nsigmx*nt(ix)
             call ccopy(nt2,u,1,ug,1)
             if(lgrad) call ccopy(nt2,u2,1,u2g,1)
             if(lrot) then
                e_and_h = .true.
                call cc_geog(ug,nt(ix),nsig,nsta,ih_tot,ie_tot,
     &             orient,theta,e_and_h)
cXX                irot = 1
ccc             theta gives angle for rotation  of coordinate system
cXX                call geogcor(ug,nt(ix),nsig,nsta,ih,orient,decl,
cXX     &            theta,irot)
                if(lgrad) then
                   e_and_h = .true.
                   call cc_geog(u2g,nt(ix),nsig,nsta,ih_tot,ie_tot,
     &               orient,theta,e_and_h)
cXX                call geogcor(u2g,nt(ix),nsig,nsta,ih,orient,decl,
cXX     &            theta,irot)
                endif
             endif

ccc          Here different options for estimating plane wave response
ccc             space need to be considered
ccc          To start with ... just use eigenvectors
ccc          Need signal power, noise power
ccc          Plane wave response with average = normal
             if(lpw) then
ccc              plane-wave response the easy way
                call mkbih(b,nt(ix),ih,ie_tot,nsta,-1)
                c(1) = 1.
                c(2) = 0.
                call anfld(ug,nt(ix),nt(ix),b,c,a,y)
                if(lgrad) call anfld(u2g,nt(ix),nt(ix),b,c,a,y2)
                c(1) = 0.
                c(2) = 1.
                call anfld(ug,nt(ix),nt(ix),b,c,a,y(nt(ix)+1))
                if(lgrad) 
     &              call anfld(u2g,nt(ix),nt(ix),b,c,a,y2(nt(ix)+1))
                call mvc(y,ug,2*nt(ix))
                if(lgrad) call mvc(y2,u2g,2*nt(ix))
             end if

cccccccccccccccccccccc  gradient  ccccccccccccccccccccccccccccccccccccc
             if(lgrad) then
ccc             make plane wave and gradient vectors for sub-array
ccc             first center and convert station coordinates to km
                call llkm(stc,nsta,stkm,ier)
ccc               check to see if all stations had station coordinates
                if(ier.eq.-1) then
                   write(0,*) ' stopping; station coordinates missing'
                   stop
                endif
ccc             now make idealized uniform, gradient vectors (H only... so far)
                nvec = nsigmx
                call mkw(stkm,nt(ix),ih,nsta,w,theta,wrk)
                call pwgfit(u2g,w,cpwg,nt(ix),nvec,nsta,ih,pwg)
                evar = 1.
                call mkuvar(ev,nt(ix),pwg,nvec,evar,uvar,nf(ix))
                call pwgslv(cpwg,uvar,u2g,v,nt(ix),nvec)
ccc             v now contains plane-wave gradient separation
             endif
cccccccccccccccccccccc  gradient  ccccccccccccccccccccccccccccccccccccc

ccc          now move simple plane wave estimate into u2g
ccc          call mvc(y2,u2g,2*nt(ix))

             if(ix.eq.1 .and. l_NoDeadChannels) then
ccc            various other calculations, outputs for main array
           
ccc             Robust remote reference .... use cleaned data from all
ccc             channels with var(i) .gt. 0 to construct weighted PC
ccc             reference fields; do robust remote ref, cleaning
ccc             outliers for both local H and E ....
                n1 = 0
                do i=1,nt_tot
                   vartmp(i) = var(i)*sndiag(i,ib)
                   if(var(i).gt.0.) then
                      n1 = n1 + 1
                      i1(n1) = i
                   endif
                enddo

ccc             loop for iterating rr_rbst with weights determined
ccc             by estimated signal/coherent noise ratio
                iter = 0
                do while(iter.le.itmax_cn)
                  iter = iter+1
                  if(iter.gt.1) then
                    call cn_wt(cjobwt,p,nt(1),nf(1),npc,nfb,xpc,u,ev,pw
     &                ,v,wt,xp,sndiag(1,ib))
                    do i = 1,nf(1)
                      lomit(i)   = wt(i) .eq. 0.0
                    enddo
                  else
                    call init_lomit(lomit,ixs,nf(1),isuse,nuset,
     &                                       id,ndmax)
                    lpwout = .false.
                    npc = nsigmx
                  endif
                  if(iter.eq.itmax_cn+1) lpwout = l_Pw
                  if(l_ref_proj) then
                     call mkbih(b,nt(ix),ih,ie_tot,nsta,irefsta)
                  endif
c                  write(0,*) 'ntape before rr_rbst',(ntape(l),l=1,6)
                  call rr_rbst(xr,xc,lomit,vartmp,nt(1),nf(1),nsta,ih,
     &              ie_tot,n1,i1,xp,itmax_sdm,itmax_rrr,orient,decl,
     &              theta,ngrp_tf,igrp_tf,
     &              lpwout,ioupw,pw,z(1,1,1,ib),sig_s(1,1,1,ib),
     &              sig_e(1,1,1,ib),ldf(1,ib),l_ref_proj,b,period(ib),
     &                 pwrmin)
c                  write(0,*) 'ntape after rr_rbst',(ntape(l),l=1,6)
                enddo
ccc   make units of E same as H for output eigenvectors
                do ivec=1,nsig
                   iu1= nt(ix)*2*(ivec-1)+1
                   call eunits(ug(iu1),nt(ix),ih,ie_tot,nsta,period(ib))
                   if(lgrad)
     &                call eunits(u2g(iu1),nt(ix),ih,ie_tot,nsta,
     &                       period(ib))
                enddo  ! ivec
ccc    scale "secondary" eigenvectors to make H of order 1
                if(lpw) then
                   isec1 = 3
                else
                   isec1 = 1
                endif
                do ivec=isec1,nsig
                   iu1= nt(ix)*2*(ivec-1)+1
                   call hnorm(ug(iu1),nt(ix),ih_tot,ie_tot,nsta)
                   if(lgrad) call hnorm(u2g(iu1),nt(ix),ih,ie_tot,nsta)
                enddo  ! ivec
ccc      print out eigenvalues/vectors etc. for one band
ccc         (main array only)
                iounit = 3
                call prteig(ev,ev2,ug,period(ib),nt(ix),nf(ix),sta,
     &             iounit,arrayid,id,iband,nsig,nsta,ih_tot,.true.,
     &             stc_tot,theta,.true.,cvec)
                if(lgrad) then
                   call prteig(ev,ev2,v,period(ib),nt(ix),nf(ix),sta,
     &             iounit,arrayid,id,iband,nsigmx,nsta,ih_tot,.true.,
     &             stc_tot,theta,.false.,cevec)
                end if
             end if     !!!!!  if(ix.eq.1)
          enddo   ! ix   sub-array loop
       enddo      !ib    frequency band loop

      close(3)
      do ix = 1,nxu
         close(ivunits(ix))
      enddo

      close(2)
      if(l_PC) close(ioupc)
      if(l_Pw) close(ioupw(1))

ccc   output signal and noise powers
      close(2)
      open(unit=2,file=arrayid(1:lar_id)//'.SN',status='unknown')
      write(2,*) nt_tot,nbt
      write(2,700) (ldf(1,ib),ib=1,nbt)
ccc   FIRST:  Eigenvalues in "Noise units"
      write(2,*) 'Eigenvalues in noise units'
      do ib = 1,nbt
         write(2,710) period(ib),(eval(k,ib),k=1,nt_tot)
      enddo
ccc   SECOND:  SIGNAL POWER
      write(2,*) 'SIGNAL POWER'
      do ib = 1,nbt
         write(2,710) period(ib),(sdiag(k,ib),k=1,nt_tot)
      enddo
ccc   THIRD INCOHERENT NOISE POWER
      write(2,*) 'INCOHERENT NOISE POWER'
      do ib = 1,nbt
         write(2,710) period(ib),(sndiag(k,ib),k=1,nt_tot)
      enddo
ccc   FOURTH: RELATIVE INCOHERENT NOISE POWER
      write(2,*) 'RELATIVE INCOHERENT NOISE POWER'
      do ib = 1,nbt
         write(2,710) period(ib),(sndiag(k,ib)/sdiag(k,ib)
     &           ,k=1,nt_tot)
      enddo
ccc   FIFTH:   RELATIVE EIGENVALUES
      write(2,*) 'Relative Eigenvalues : ORDER CHANGED!!!!'
      do ib = 1,nbt
         write(2,710) period(ib),(eval2(k,ib)/tpower(ib),k=1,nsigmx)
      enddo
ccc   SIXTH:  Eigenvalues in nT  (only dominant NSIGMX)
      write(2,*) 'Eigenvalues in nT'
      do ib = 1,nbt
         write(2,710) period(ib),(eval2(k,ib),k=1,nsigmx)
      enddo

      if(l_CC) then
ccc   output canonical covariances
         open(unit=2,file=arrayid(1:lar_id)//'.CC',status='unknown')
         write(2,*) nt_tot,nbt
         do ib = 1,nbt
            write(2,710) period(ib),(cc(k,ib),k=1,nt_tot)
         enddo
      endif

700   format(20i6)
710   format(16e10.3)
720   format(e10.3,15f7.4)

      if(l_NoDeadChannels) then
ccc   output "impedances" (TFs) with covariances
      chead = ' '
      do i = 1,nsta_tot
         dr1(i) = dr(1,i)
      enddo
      call wrt_z(z,sig_s,sig_e,nbt,ngrp_tf,igrp_tf,cgrp,period,
     &  stc_tot,decl_tot,orient,chid_grp,sta,sta_grp,idl,ibandlim,
     &  arrayid,chead,ldf,dr1)
      endif

      stop
      end
c______________________________________________________________________
c
      subroutine printsym(s,n)
      complex s(*)
      i1 = 1
      do i = 1,n
         i2 = i1+i-1
         write(6,'(15f12.4)') (real(s(j)),j=i1,i2)
         write(6,'(15f12.4)') (aimag(s(j)),j=i1,i2)
         i1 = i1 + i
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine prtx(x,n)
      complex x(n,*)
      do i = 1,20
         write(0,*) (x(j,i),j=1,n,3)
      enddo
      return
      end
c____________________________________
c
      subroutine mvc(y,u,n)
      complex y(n),u(n)
      do i = 1,n
         u(i) = y(i)
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine wrt_v(iounit,irec,v,nt,icov,sigpwr,noisepwr)

      complex v(nt,2)
      real sigpwr(*),noisepwr(*)

      if(icov.eq.0) then
ccc         isotropic error variance; just output overall error scale
         write(iounit,rec=irec) v,sigpwr(1),sigpwr(2),noisepwr(1),
     &       noisepwr(2)
      else
         write(0,*) 'dont know how to deal with this option in wrt_v'
      end if
      return
      end
c______________________________________________________________________
c
      subroutine evout(ioev,period,u2,nt,nsig)
      complex u2(nt,nsig)
      real period
      write(ioev,'(e12.4,i6)') period,nsig
      do i = 1,nsig
         write(ioev,'(10f8.4)') (u2(j,i),j=1,nt)
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine wrt_uev(iouev,irec,u,ev,var,nt,nsig,nf,period)
      complex u(nt,nsig)
      real ev(nsig),var(nt),period
      integer nf
      write(iouev,rec=irec) period,nf,ev,var,u
      return
      end
c______________________________________________________________________
c
      subroutine wrt_s(iouev,irec,s,ns,var,nt,nf,period)
      complex s(ns)
      real var(nt),period
      integer nf
      write(iouev,rec=irec) period,nf,var,s
      return
      end
c______________________________________________________________________
c
      function n_eig_sig(nf,eval,nt,nsigmx)
      integer n_eig_sig,nt
      real eval(nt),sigmin
      parameter (sigmin = 3.,nsigmin=2)
ccc   should figure out how to choose sigmin based on sample size ....
      nsig = 0
      do i = 1,nt
         if(eval(i).ge.sigmin) nsig = nsig+1
      enddo
      nsig = max(nsig,nsigmin)
      n_eig_sig = min(nsig,nsigmx)
      return
      end
ccc_____________________________________________________________________
ccc
      subroutine init_lomit(lomit,ixs,nf,isuse,nuset,id,ndmax)
      integer isuse(2,ndmax,*),ixs(nf),nf,itemp(80),nuset,ndmax,id
      logical lomit(nf),lu

c      write(91,*) 'In init_lomit'
c      write(91,*) 'isuse = ',isuse

      do i = 1,nf
         lu = isuse(1,id,1).le.ixs(i) .and. isuse(2,id,1).ge.ixs(i)
         do k = 2,nuset
            lu = lu .or. 
     &       (isuse(1,id,k).le.ixs(i).and.isuse(2,id,k).ge.ixs(i))
         enddo
         lomit(i) = .not.lu
      enddo
      do i = 1,nf,80
         i2  = min(i+79,nf)
         do j = i,i2
            if(lomit(j)) then
               itemp(j-i+1)=0
            else
               itemp(j-i+1)=1
            endif
         enddo
c         write(91,'(80i1)') (itemp(j),j=1,i2-i+1)
      enddo
      return
      end
ccc_____________________________________________________________________
ccc
      subroutine usage(itmax_ln_d,itmax_cln_d,itmax_rrr_d)

      integer itmax_ln_d,itmax_cln_d,itmax_rrr_d

      write(0,*) 'Usage: multmrn    OPTIONS:'

      write(0,*) '-n : no robust features'

      write(0,*) '-m : only robust SDM'

      write(0,*) '-f<file> : change cf_array file to <file>'

      write(0,*) '-s<file> : change default grouping of channels'
      write(0,*) '           into stations'

      write(0,*) '-t : output local noise array TF files'

      write(0,*) '-u : assume signal and coh noise are uncorrelated'

      write(0,*) '-g : use spatial gradient model'

      write(0,*) '-N : dont transform eigenvectors in M_ file'

      write(0,*) '-c<?> : <?> = H for canonical covariance in units'
     &    ,' of mag fields'

      write(0,*) '      : <?> = N for canonical covariance in units'
     &    ,' of noise;'

      write(0,*) '      : <?> = R for canonical coherence squared'
      write(0,*) '               (default is nothing)'

      write(0,*) '-i : change default max # of iterations for robust'
      write(0,*) '     regression for each inner loop estimate of '
      write(0,*) '     local noise; Default is ',itmax_ln_d

      write(0,*) '-o : change default # of outer loops for local'
      write(0,*) '      noise estimation/individual channel outlier'
      write(0,*) '      cleanup; Default is ', itmax_cln_d

      write(0,*) '-I : change default max # of iterations of iterations'
      write(0,*) '      used for robust regression for final TF '
      write(0,*) '      estimate;  Default is ', itmax_rrr_d

      write(0,*) '-w : change default # of iterations for coherent'
      write(0,*) '         noise downweighting'

      write(0,*) '-p : change coherent noise cutoff'

      write(0,*) '-a : omit data set ranges specified in cf_array file'
      write(0,*) '     only on first coherent noise iteration'

      write(0,*) '-P : output principal component FCs in PC_ file'

      write(0,*) '-r : output raw data in PC_ file'

      write(0,*) '-z : rotate channels into common coordinate system'
      write(0,*) '     before outputing M_**** file'

      write(0,*) '-C : output correlation matrices (total+local noise'
      write(0,*) '       correlations)'

      write(0,*) '-L : dont output Pw_ file'
      write(0,*) '     (plane wave response space estimates)'

      write(0,*) '-T : turn ON automatic timing error correction'

      write(0,*) '-GX: change default channel grouping for incoherent'
      write(0,*) '   : noise variance estimation. Here X = T,S, or A  '
      write(0,*) '   : T => all components of a type at a site (def.)'
      write(0,*) '   : S => all components at a single site '
      write(0,*) '   : A => each component by itself'

      write(0,*) '-R#: Use projection of magnetic fields from '
      write(0,*) '      station # into coherent signal/noise space '
      write(0,*) '      to define plane wave reference'

      write(0,*) '-s<file_name> change default channel groupings for'
      write(0,*) '   single station TF output files (Z*****)'
      write(0,*) '   default is to use same grouping as FC files'

      write(0,*) '-S<decimation_level> make a list of set numbers'
      write(0,*) '   available for specified decimation level and quit'
      write(0,*) '-M<pwr_min> Only use sets with sig pwr (determined'
      write(0,*) '  by 2 dominant eigenvectors) exceeding pwr_min'

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
      if(c(1:3).eq.'NaN')inan=.true.
      return
      end
