      program sys_tbl_tf
c     stripped-down program for computing one
c     system parameter correction table, using system
c     parameters defined in an *.sp file, outputting into a *.tfs file

      include 'fcor.inc'
      integer nch,nfil(nchmx),iftype(nfilmax,nchmx)
      real sampr,sc(nchmx),decl,stcor(2),orient(2,nchmx),cda,cdb
      logical lclkd
      character*6 chid(nchmx)
      character*80 cfsp,cfsTF,sensdir
      character*120 afparam(nfilmax,nchmx)

      read(5,'(a80)') cfsp
      read(5,'(a80)') cfsTF
      read(5,'(a80)') sensdir

c  read system parameter file:  number of channels,
c  electrode line lengths, filter parameters, conversion
c  factors, etc. from file *.sp
      iounit=3
      open(unit=iounit,file=cfsp)
      call getsp1(nch,sampr,sc,nfil,iftype,afparam,
     &   iounit,decl,stcor,orient,chid,cda,cdb,lclkd)

      write(0,*) 'Done with getsp1'

ccc   make table with only system corrections
ccc   and output to file *.fts
      call systblsu(sampr,nch,nfil,afparam,iftype,sc,cda,cfsTF,sensdir)

      stop
      end
