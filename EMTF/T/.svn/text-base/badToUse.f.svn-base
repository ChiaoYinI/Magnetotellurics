c
      subroutine BadToUse(badRecFile,brUnit,ndMax,nd,dr,nuset,isuse,
     &   npts,fac)
        
ccc   reads badRecFile (format used for dnff) and converts to array
ccc   of set numbers to use (isuse) for tranmt and multmtrn
ccc   Ignores type of bad data flag (all segments listed are omitted)
ccc   Warning: to convert from record numbers to set, the # of points
ccc    per/set for the first site are used, and npts*fac is assumed
ccc    to be the offset between data sets. NOrmally fac is 0.75,
ccc    but this could be different if different overlaps were specified
ccc    in the cf_decset file used for FFT'ing the data
ccc   Also note that this routine just sets up set number limits for the
ccc    first decimation level.

      integer isuse(2,ndMax,*),nuset,ndMax,brUnit,nseg,npts(ndMax)
      integer nd,id,n
      character*80 badRecFile
      integer maxint

      parameter (nbmx = 500,ndmx=10)
      real irecb(2,nbmx),fac,dr(ndMax,*),dtBR

      maxint = 2**30 + (2**30 -1)
      write(0,*) badRecFile, brUnit
      open(unit = brUnit, file = badRecFile,status='old',err=200)
        
ccc   dtBR is sampling rate used when bad record segments
ccc   were marked (may differ from actual sampling rate dr(1,1) ...
ccc    due to possible use of decimated data files for marking bad
ccc    segments.
      read(brUnit,*) nseg,dtBR
      do 10 i = 1,nseg
         read(brUnit,*) irecb(1,i),irecb(2,i)
         irecb(1,i) = irecb(1,i)*(dtBR/dr(1,1))
         irecb(2,i) = irecb(2,i)*(dtBR/dr(1,1))
10       continue
      close(brUnit)

      nuset = nseg+1
      isuse(1,1,1) = -maxint
      isuse(2,1,nuset) = maxint
      write(0,*) 'Omitting level one sets :'
ccc   why k+1 below ???  (Because isuse gives the range of set numbers to use)
      do k = 1,nseg
         n = npts(1)*fac
         do id = 1,nd
            isuse(2,id,k) = int(irecb(1,k)/n-(1-fac)/fac)
            isuse(1,id,k+1) = int(irecb(2,k)/n+.999)+1
            n = n*dr(id+1,1)/dr(id,1)
         enddo
         write(0,*) isuse(2,1,k)+1, ' - ', isuse(1,1,k+1)-1
      enddo

      return

200   nuset = 0
      write(0,*) 'Error opening bad record file'
      return
      end
