
      subroutine wrt_sp(fid,stname,stcor,decl,nch,samprate,clock,
     &  chid,orient,dipole_length,gain,sensitivity,nf,nfmax,ftype,
     &   fparam)

ccc   writes out sp******* file in standard format ... 
ccc   dipole_length is in  km
ccc   
      integer fid,nch,nf(nch),nfmax,ich,k
      real stcor(2),decl,samprate,orient(nch),gain(nch),clock(2),
     &     sensitivity(nch),dipole_length(nch)
      character*4 stname
      character*2 ftype(nfmax,nch)
      character*256 fparam(nfmax,nch)
      character*6 chid(nch)
      write(fid,'(a20)') stname
      write(fid,'(2f10.4)') stcor
      write(fid,'(f10.4)') decl
      write(fid,'(i2)') nch
      write(fid,'(e12.4)') samprate
      write(fid,'(2f10.4)') clock
      do ich = 1,nch
         write(2,'(a6)') chid(ich)
         if(chid(ich)(1:1).eq.'E') then
            write(fid,'(4f10.4)') dipole_length(ich),orient(ich),0.,
     &            gain(ich)
         else
            write(fid,'(2f10.4)') orient(ich),0.
         endif
         write(fid,*) sensitivity(ich),nf(ich)
         do k = 1,nf(ich)
            write(fid,'(a2)') ftype(k,ich)
            write(fid,'(a)') fparam(k,ich)
         enddo 
      enddo
      return
      end
