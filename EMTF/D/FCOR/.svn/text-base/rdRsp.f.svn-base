      subroutine rdRsp(cfile,sta,chid,qChid,rlatlon,azimuth,dt,nzeros,
     &      npoles,zeros,poles,a0,sens)

      character*80 cline
      character*256 cfile
      character*3 qChid
      character*4 sta
      character*6 chid
      real rlatlon(2),azimuth,dt,a0,sens,temp(2)
      complex zeros(*),poles(*)
      integer nzeros,npoles,iunit,nCline,nword,i1,i2,kz,i0
      iunit = 1
      open(unit=iunit,file=cfile,status='old')
      read(1,'(a3,1x,a6)') qChid, chid
      nCline = 80
ccc   station code:
      read(1,'(a)') cline
      nWord = -1
      call findWord(cline,nCline,nword,i1,i2)
      read(cline(i1:i2),'(a)') sta
c      write(0,*) 'sta',sta

ccc   Lat
      read(1,'(a)') cline
      nWord = -1
      call findWord(cline,nCline,nword,i1,i2)
      read(cline(i1:i2),*) rlatlon(1)
c      write(0,*) 'rlatlon',rlatlon

ccc   Long
      read(1,'(a)') cline
      nWord = -1
      call findWord(cline,nCline,nword,i1,i2)
      read(cline(i1:i2),*) rlatlon(2)
      if(rlatlon(2) .lt.0) then
         rlatlon(2) = rlatlon(2)+360
      endif
c      write(0,*) 'rlatlon',rlatlon

ccc   azimuth
      read(1,'(a)') cline
      nWord = -1
      call findWord(cline,nCline,nword,i1,i2)
      read(cline(i1:i2),*) azimuth
c      write(0,*) 'azimuth',azimuth

ccc   dt
      read(1,'(a)') cline
      nWord = -1
      call findWord(cline,nCline,nword,i1,i2)
      read(cline(i1:i2),*) dt
      dt = 1./dt
c      write(0,*) 'dt',dt

ccc   TF type
      read(1,'(a)') cline
c      write(0,*) 'TF type',cline

ccc   Response in units lookup
      read(1,'(a)') cline

ccc   A0
      read(1,'(a)') cline
      nWord = -1
      call findWord(cline,nCline,nword,i1,i2)
      read(cline(i1:i2),*) a0
c      write(0,*) 'a0 ',a0

ccc   Normalization frequency
      read(1,'(a)') cline
c      write(0,*) 'normalization freq ',cline

ccc   nzeros
      read(1,'(a)') cline
      nWord = -1
      call findWord(cline,nCline,nword,i1,i2)
      read(cline(i1:i2),*) nzeros
c      write(0,*) 'nzeros ',nzeros

ccc   npoles  
      read(1,'(a)') cline
      nWord = -1
      call findWord(cline,nCline,nword,i1,i2)
      read(cline(i1:i2),*) npoles
c      write(0,*) 'npoles ',npoles
   
ccc   zeros
      do kz = 1,nzeros
         read(1,'(a)') cline
         nWord = -3
         call findWord(cline,nCline,nword,i1,i2)
         nWord = -4
         call findWord(cline,nCline,nword,i1,i0)
         read(cline(i1:i2),*) temp
         zeros(kz) = cmplx(temp(1),temp(2))
      enddo
c      write(0,*) 'zeros ',(zeros(kz),kz=1,nzeros)

ccc   poles
      do kz = 1,npoles
         read(1,'(a)') cline
         nWord = -3
         call findWord(cline,nCline,nword,i1,i2)
         nWord = -4
         call findWord(cline,nCline,nword,i1,i0)
         read(cline(i1:i2),*) temp
         poles(kz) = cmplx(temp(1),temp(2))
      enddo
c      write(0,*) 'poles ',(poles(kz),kz=1,npoles)

ccc   sensitivity
      read(1,'(a)') cline
      nWord = -1
      call findWord(cline,nCline,nword,i1,i2)
      read(cline(i1:i2),*) sens
c      write(0,*) 'sensitivity ',sens

      close(iunit)
      return
      end
ccc
ccc
ccc
      subroutine findWord(cline,nCline,nword,i1,i2)
ccc   finds range of characters in string containg
ccc   word number nword; words are assumed separated
ccc   by blanks
      integer nCline,nWord,i1,i2,n,nEnd,nBeg
      character*80 cline
      logical blanks

      nEnd = 0
      nBeg = 0
      blanks = .true.
      if(nword .lt. 0) then
ccc      count words from end of string
        n = nCline
        nWord = abs(nWord)
   
        do while (n .gt. 0)
           if(cline(n:n).ne.' ') then
               if(blanks) then
                  nEnd = nEnd + 1
                  blanks = .false.
                  if(nEnd .eq. nWord) then
                     i2 = n
                  endif
               endif
            else
               if(.not.blanks) then
                  nBeg = nBeg + 1
                  blanks = .true.
                  if(nBeg .eq. nWord) then
                     i1 = n
                  endif
               endif
            endif
            n = n-1
         enddo          
      else 
ccc      count words from begining of string
        n = 1
   
        do while (n .gt. 0)
           if(cline(n:n).ne.' ') then
               if(blanks) then
                  nBeg = nBeg + 1
                  blanks = .false.
                  if(nBeg .eq. nWord) then
                     i1 = n
                  endif
               endif
            else
               if(.not.blanks) then
                  nEnd = nEnd + 1
                  blanks = .true.
                  if(nEnd .eq. nWord) then
                     i2 = n
                  endif
               endif
            endif
            n = n+1
         enddo          
      endif 
      return
      end
