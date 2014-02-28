program lpNIMS
!   modified from decNIMS to low-pass filter and decimate
!   1 hz data to a very low sampling rate to facilitate
!   exploration of electric field noise issues

implicit none

    integer (kind=4),parameter	:: fidIn = 11
    integer (kind=4),parameter	:: fidOut = 12
    integer (kind=4),parameter	:: nch = 5
    integer (kind=4),parameter	:: nblk = 256
    integer (kind=4),parameter	:: ndec = 500
    integer (kind=4),parameter	:: nfil2 = 501 
    integer (kind=4)		:: npts,msval,nHdrBytes,nHdrRec
    integer (kind=4)		:: nIn,nOut,jj,j,k,j1,j2,k
    integer (kind=4)		:: nExtra,nOutBlk,nScansOut
    integer (kind=4)		:: outRecLength,i,nHdrBytesOut,nDataBlk
    integer (kind=4)		:: iargc,narg
    integer (kind=4)		:: irlong
    real (kind=4)		:: dr,drOut,wsum,xsum
    real (kind=4)		:: fil(-nfil2:nfil2),sum
    character*80		:: cfileIn, cfileOut

    integer (kind=4), allocatable, dimension(:)	:: header
    integer (kind=4), allocatable, dimension(:,:)	:: ix, iy

    !  filter coefficients for low-pass filter (half cosine)
!      fil = (/.02,.04,.057,.071,.083,.090,.092,.092, &
!  		.092,.090,.083,.071,.057,.04,.02/)
    !  set up filter coefficients ...
    sum = 0.0
    do k = -nfil2,nfil2
       fil(k) = cos(3.14159*k/(2*nfil2))
       sum = sum+fil(k)
    enddo 
    do k = -nfil2,nfil2
       fil(k) = fil(k)/sum
    enddo 

    !  set input and output files
    if(iargc() .eq. 1) then
       !  strip off root
       call getarg(1,cFileIn)
       j = irlong(cFileIn)
       !  add .dbn to root for input
       cFileIn = cFileIn(1:j)//'.dbn'
       !  add .lpb to root for low-pass decimated binary output
       cFileOut = cFileIn(1:j)//'.lpb'
    elseif(iargc() .eq. 2) then
       call getarg(1,cfileIn)
       call getarg(2,cfileOut)
    else
       write(0,*) 'Usage: decNIMS input_file <output_file>'
       write(0,*) '       if only one argument output file'
       write(0,*) '       name is constructed from input'
       write(0,*) '       by appending .dbn to root'
       stop
    endif

    !  first get key info from NIMS binary file header
    open(file=cfileIn,unit=fidIn,access = 'direct',recl=4,&
		status = 'OLD')
    read(fidIn,rec=1) nHdrBytes 
    read(fidIn,rec=5) dr 
    read(fidIn,rec=19) nIn 
    read(fidIn,rec=21) msval 
    write(0,*) 'nHdrBytes,dr,nIn,msval', &
		nHdrBytes,dr,nIn,msval
    
    nHdrRec = nHdrBytes/4+3

    outRecLength = nblk*nch
    nHdrbytesOut = (outRecLength-3)*4
    !  open file for output
    open(unit=fidOut,file=cfileOut,access = 'direct',recl=outRecLength*4)
    !  copy header (including fortran sequential file markers,
    !    at start and end of header record + at start of single
    !    data record 
    !    NOTE: initially sampling rate in header of output file
    !      is set to that of input file
    allocate(header(outRecLength))
    header = 0
    do i = 1,nHdrRec
       read(fidIn,rec=i) header(i)
    enddo
    write(fidOut,rec=1) header
    
    !  close input file, reopen as sequential direct access
    close(fidIn)
    open(file=cfileIn,unit=fidIn,form = 'unformatted')
    read(fidIn)
    
    write(0,*) 'nch,nIn',nch,nIn

    ! allocate for and read input array
    allocate(ix(nch,nIn))
    read(fidIn) ix
    close(fidIn)

    ! allocate for output
    !   to make sure all data is output, allocate enough 
    !     to fill out last block with missing data
    nOut = int((nIn+1)/ndec)
    drOut = dr*ndec
    nOutBlk = nOut/nblk
    if(nOutBlk*nblk .lt. nOut) then
       nOutBlk = nOutBlk+1
    endif
    nExtra = nOutBlk*nblk - nOut 
    nOut = nOutBlk*nblk
    allocate(iy(nch,nOut))
   
    write(0,*) 'nOut,drOut,nOutBlk,nExtra', &
		nOut,drOut,nOutBlk,nExtra

    !  filter and decimate
    do i = 1,nch
       jj = 0
       do j = 1,nIn,ndec
          wsum = 0.0
          xsum = 0.0
          jj = jj+1
          j1 = max(j - nfil2,1)
          j2 = min(j + nfil2,nIn)
          do k = j1,j2
             if(ix(i,k).ne. msval) then
                wsum = wsum + fil(k-j)
                xsum = xsum + fil(k-j)*ix(i,k)
             endif
          enddo
          if(wsum.gt. 0.5) then
             iy(i,jj) = int(xsum/wsum)
          else
             iy(i,jj) = msval
          endif
       enddo
       !  pad end of last block with missing value code
       nScansOut = jj
       do j = nScansOut+1,nScansOut+nExtra
          iy(i,j) = msval
       enddo
    enddo

    do i = 1,nOutBlk
       j1 = (i-1)*nblk+1
       j2 = i*nblk
       write(fidOut,rec = i+1) ((iy(k,j),k=1,nch),j=j1,j2) 
    enddo

    nDataBlk = nOutBlk*outRecLength*4
    
    !   close output file, reopen with record length of 4, 
    !    output corrected header info
    close(fidOut)
    open(file=cfileOut,unit=fidOut,access = 'direct',recl=4)
    write(fidOut,rec=1)  nHdrBytesOut
    write(fidOut,rec=5) drOut 
    write(fidOut,rec=19) nScansOut 
    write(fidOut,rec=outRecLength-1)  nHdrBytesOut
    write(fidOut,rec=outRecLength)  nDataBlk
    write(fidOut,rec=outRecLength*(nOutBlk+1))  nDataBlk
    close(fidOut)

end program 
!___________________________________________________________
!
      function irlong(cstr)
!   find number of characters before first '.' or ' ' in a string
      integer(kind=4) irlong
      character(len=80),intent(in) :: cstr
      integer(kind=4) ::	n
      n = len_trim(cstr) 
      do 10 i = 1,n
      if((cstr(i:i).eq.' ').or.(cstr(i:i).eq.'.')) then
         irlong = i-1
         return
      end if
10    continue
      irlong = n
      return
      end function
