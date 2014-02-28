c__________________________________________________________
c
        subroutine bdrcsu(nd,nbdmx,ndp,idp,nds,ids,ioff,
     &       dt,time_int)
        
        include 'iounits.inc'

c       idp - array of record ranges to not process
c       ids - array of record ranges to not stack
c       ndp - array giving number of record ranges for each level
c       nds -   "     "      "            (for ids)
c       there is a seperate set of ranges for each level; level
c       zero corresponds to all levels  --   i.e. records specified in
c       level zero of idp ( ids) are not proceesed ( stacked)
c       at any level

c       now also gets record offset (i.e. the number which should be
c       added to the record numbers in the c. file to make them
c       correct (i.e. agree wtih other stations)); note: the offset
c       is added to the bad record intervals idp and ids; these should thus
c       be specified with the (erroneous) records which are on the plot

c       info is read from file brsta000x
      integer ndp(0:nd),nds(0:nd),idp(0:nd,2,nbdmx),ids(0:nd,2,nbdmx)
      integer time_int
cnew      character*80 cfbr

      parameter (nbmx = 200,ndmx=8)
c       need to change parameters in other subroutine in file !!!
      integer iflag(nbmx)
      integer nbr(2,0:ndmx),ibr(2,0:ndmx,2,nbmx) !,imap(2,4)
      real irecb(2,nbmx),dtBadRec,dt,multFac
      character*80 cstring

      open(unit = br_unit, file = cfbr,status='old',err=200)
        
c     new verision : also tries to read in sampling rate used
c     to define bad record times relative to clock zero; this
c     is required because plotting program might use decimated
c     data files for greater efficiency
c     ALSO: bad record usage was never adapted to work with 
c      the new approach to computing set numbers: record numbers
c      are relative to start of file; set offsets are computed
c      from time_int (number of seconds between "clock zero" and
c      start of data)
      read(br_unit,'(a)') cstring
      read(cstring,*,err=5) nseg, dtBadRec
      multFac = dtBadRec/dt
      go to 10
5     multFac = 1   
      read(cstring,*) nseg
10    continue
      do i = 1,nseg
         read(br_unit,*) irecb(1,i),irecb(2,i),iflag(i)
         irecb(1,i) = int(multfac*(irecb(1,i)-time_int/dtbadRec))
         irecb(2,i) = int(multfac*(irecb(2,i)-time_int/dtBadRec))
      enddo

      call mkbr(nseg,irecb,iflag,nd,nbr,ibr)
      
      read(br_unit,*,end = 112) ioff
      go to 115
112   print*,'no offset in br file; assuming zero offset'
      ioff = 0
115   continue
      close(br_unit)

c       add offset to bad record intervals
      do 125 id = 0,nd
         ndp(id) = nbr(1,id)
         do 120 j = 1,ndp(id)
            idp(id,1,j) = ibr(1,id,1,j) + ioff
            idp(id,2,j) = ibr(1,id,2,j) + ioff
120         continue
         nds(id) = nbr(2,id)
         do 122 j = 1,nds(id)
            ids(id,1,j) = ibr(2,id,1,j) + ioff
            ids(id,2,j) = ibr(2,id,2,j) + ioff
122         continue
125      continue
      return

c       this block is for no bad data (i.e., no bad record file found)
200   continue
      do 210 id = 0,nd
         ndp(id) = 0
210      continue
      do 310 id = 0,nd
         nds(id)= 0
310      continue
      ioff = 0
      close(69,status = 'delete')
      print*,'No bad record file found for this file;  ',
     &'Using all data '
      return
      end
c______________________________________________________________________
c
      subroutine mkbr(nseg,irecb,iflag,nd,nbr,ibr)
      parameter (nbmx=200,ndmx=8)
      integer iflag(nseg)
      integer nbr(2,0:ndmx),ibr(2,0:ndmx,2,nseg),imap(2,4)
      real irecb(2,nseg)
      data imap/1,0, 2,0, 2,4, 1,0/
c      imap(1,.) = 1 to not process, 2 to flag as bad
c      imap(2,.) = 0 for all levels, else = level to omit (only allows one
c             level or all for now)
c      imap(1,.) is for magnetics bad
c      imap(2,.) is for electrics bad
c      imap(3,.) is for long period bad
c      imap(4,.) is for all bad
c      My choices of imap for these four "buttons" are given in data statment
c        for imap

      do 5 i=0,nd
      do 5 k = 1,2
         nbr(k,i) = 0
5        continue

      do 100 j = 1,nseg
         nbr(imap(1,iflag(j)),imap(2,iflag(j))) = 
     &       nbr(imap(1,iflag(j)),imap(2,iflag(j))) + 1
         do 50 k = 1,2
            ibr(imap(1,iflag(j)),imap(2,iflag(j)),k,
     &             nbr(imap(1,iflag(j)),imap(2,iflag(j)))) = irecb(k,j)
50           continue
100      continue
 
      return
      end 
