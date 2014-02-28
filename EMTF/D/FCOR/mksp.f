      program mksp
ccc   program for making *.sp files from SEED volume header information
ccc   returned by running cshell script getSysRsp
ccc   
ccc   input to this program is via standard input, nch+2 records:
ccc
ccc   1) output *.sp file name
ccc   2) # of data channels
ccc   3-nch+2) *.rsp file names for each channel

      integer nchmax,nzmax,nfmax
      parameter (nchmax = 20,nzmax=10,nfmax=1)
      character*256 cfile,fparam(nfmax,nchmax),cfout
      character*4 sta
      character*2 ftype(nfmax,nchmax)
      character*3 qChid(nchmax)
      character*6 chid(nchmax)
      integer*4 nch,ich,nz,np,kz,kp,nf(nchmax)
      integer ktrim,k
      
      real*4 stcor(2),decl,dt,orient(nchmax),gain(nchmax),clock(2),
     &    sensitivity(nchmax),dipole_length(nchmax),a0,sens,clock
      complex*8 zeros(nzmax),poles(nzmax)


       read(5,'(a)') cfout
       read(5,'(i4)') nch
       do ich = 1,nch
          read(5,'(a)') cfile
           k=ktrim(cfile)
           write(0,*)cfile(1:k)
          call rdRsp(cfile,sta,chid(ich),qChid(ich),stcor,orient(ich),
     &      dt,nz,np,zeros,poles,a0,sens)
          if(chid(ich)(1:1).eq.'E' .or. chid(ich)(1:1).eq.'e') then
             sensitivity(ich) = 1.e+6/sens
          else
             sensitivity(ich) = 1.e+9/sens
          endif
          write(fparam(1,ich)(1:16),'(2i2,1x,e11.5)') nz,np,a0
          write(fparam(1,ich)(18:256),*) (zeros(kz),kz=1,nz),
     &        (poles(kp),kp=1,np)

          nf(ich) = 1
          ftype(1,ich) = 'PZ'
          gain(ich) = 1
          dipole_length(ich) = 1.
       enddo
     
       decl = 0
       clock(1) = 0
       clock(2) = 0
       open(unit = 2,file=cfout,status='unknown')
       call wrt_sp(2,sta,stcor,decl,nch,dt,clock,chid,orient,
     &  dipole_length,gain,sensitivity,nf,nfmax,ftype,fparam)
       close(2)
       end
c
       integer function ktrim(str)
       character*256 str
       integer k
       ktrim=1
       do k=2,256
        if(str(k:k).eq.' ')then
         ktrim=k-1
         return
        endif 
       enddo
       return
       end 
