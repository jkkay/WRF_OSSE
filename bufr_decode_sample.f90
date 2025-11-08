program bufr_decode_sample
!
! example of reading observations from bufr
!
 implicit none

 integer, parameter :: mxmn=35, mxlv=1
 integer :: ireadmg,ireadsb
 character(len=40)            :: obstr,hdstr,qmstr,oestr, pcstr

 integer :: unit_in=10
 integer :: idate,iret,num_message,num_subset,num_msg_sb
 integer :: num_msg_sb_m3, num_msg_sb_00, num_msg_sb_p3

 real :: max_lat, min_lat, max_lon, min_lon 
!
!  namelist files
!
  character(180) :: infilename        ! file from GSI running directory
  character(180) :: outfilename_00       ! file name saving results
  character(180) :: outfilename_m3       ! file name saving results
  character(180) :: outfilename_p3       ! file name saving results
  character(180) :: obs_typ, obs_date
  namelist/iosetup/ infilename, outfilename_00, outfilename_m3, outfilename_p3, obs_typ, obs_date

 character(len=8)             :: subset
 character(len=14)            :: cdate, dmn, platform_name
 real*8             :: hdr(8)
 real*8             :: pmo(2,1)
 real*8             :: obs(8,255),qms(8,255),oes(8,255),pco(8,255)
 real               :: woe,toe,qoe,poe,pob,pwe
 real*8             :: obs_time
 integer              :: iyear, imonth, iday, ihour, imin

 real*8, allocatable :: data_obs_00(:,:,:,:), data_hdr_00(:,:)
 real*8, allocatable :: data_obs_m3(:,:,:,:), data_hdr_m3(:,:)
 real*8, allocatable :: data_obs_p3(:,:,:,:), data_hdr_p3(:,:)
 character(len=100) :: file_out

 real*8 :: a(2,4)
 integer :: i,j

 hdstr='SID XOB YOB DHR TYP ELV T29'
 obstr='POB QOB TOB ZOB UOB VOB PWO CAT' ! observation
 qmstr='PQM QQM TQM ZQM WQM NUL PWQ NUL' ! quality marker
 oestr='POE QOE TOE NUL WOE NUL PWE NUL' ! observation error
 pcstr='PPC QPC TPC ZPC WPC NUL PWP NUL' ! program code

! lat/lon
 min_lat=30 !27.7657
 max_lat=50 !51.0208

 min_lon=-123 !-127.598
 max_lon=-87 !-82.402

! decode
  open(11,file='namelist.conv')
   read(11,iosetup)
  close(11)
  print*, 'jk1'
! total number of each observation type
 open(unit_in,file=trim(infilename),action='read',form='unformatted')
 call openbf(unit_in,'IN',unit_in)
 call datelen(10)
   num_message=0
   num_msg_sb=0
   num_msg_sb_m3=0
   num_msg_sb_00=0
   num_msg_sb_p3=0
   msg_report0: do while (ireadmg(unit_in,subset,idate) == 0)
     num_message=num_message+1
     num_subset = 0
     !write(*,'(I10,2I8,a10)') idate,num_message,num_msg_sb,subset
     sb_report0: do while (ireadsb(unit_in) == 0)
       num_subset = num_subset+1
       call ufbint(unit_in,hdr,7,1 ,iret,hdstr)
       call ufbint(unit_in,obs,8,255,iret,obstr)
       call ufbint(unit_in,qms,8,255,iret,qmstr)
       call ufbint(unit_in,oes,8,255,iret,oestr)

       if (hdr(2) >= 180.0) hdr(2) = hdr(2) - 360.0

       if ((trim(subset) .eq. trim(obs_typ)) .and. &
           (hdr(2) .gt. min_lon) .and. (hdr(2) .lt. max_lon) .and. &
           (hdr(3) .gt. min_lat) .and. (hdr(3) .lt. max_lat) ) then
           !(hdr(3) .gt. min_lat) .and. (hdr(3) .lt. max_lat) .and. & 
              if ((hdr(4) .lt. -1.5)) then
                 num_msg_sb_m3 = num_msg_sb_m3+1
              else if ((hdr(4) .ge. -1.5) .and. (hdr(4) .lt. 1.5)) then
                 num_msg_sb_00 = num_msg_sb_00+1
              else if ((hdr(4) .ge. 1.5)) then
                 num_msg_sb_p3 = num_msg_sb_p3+1
              end if 
       end if

      ! write(*,'(2I10,f15.1)') num_subset, num_msg_sb, hdr(8) 
     enddo sb_report0
   enddo msg_report0
 call closbf(unit_in)
 write(*,"(A,5I10)") 'num_msg_sb=',num_message,num_msg_sb,num_msg_sb_m3,num_msg_sb_00,num_msg_sb_p3 

 open(10,file=trim(obs_date)//"_"//trim(obs_typ)//"_obs_num.txt",status='new')
 write(10,'(3I10)') num_msg_sb_m3, num_msg_sb_00, num_msg_sb_p3
 close(10)
!----------------------------------------------------------------------
! READ OBS1 : -1.5 hr <= OBS < 1.5 hr
 allocate(data_hdr_00(num_msg_sb_00,8))
 allocate(data_obs_00(num_msg_sb_00,4,8,255))
 data_hdr_00 = 100000000000.0
 data_obs_00 = 100000000000.0

 open(unit_in,file=trim(infilename),action='read',form='unformatted')
 call openbf(unit_in,'IN',unit_in)
 call datelen(10)
   num_message=0
   num_msg_sb_00=0
   msg_report_00: do while (ireadmg(unit_in,subset,idate) == 0)
     num_message=num_message+1
     num_subset = 0
     write(*,'(I10,I8,a10)') idate,num_message,subset
     sb_report_00: do while (ireadsb(unit_in) == 0)
       num_subset = num_subset+1
       call ufbint(unit_in,hdr,7,1 ,iret,hdstr)
       call ufbint(unit_in,obs,8,255,iret,obstr)
       call ufbint(unit_in,qms,8,255,iret,qmstr)
       call ufbint(unit_in,pco,8,255,iret,pcstr)
       call ufbint(unit_in,oes,8,255,iret,oestr)

       if (hdr(2) >= 180.0) hdr(2) = hdr(2) - 360.0

       if ((trim(subset) .eq. trim(obs_typ)) .and. &
           (hdr(2) .gt. min_lon) .and. (hdr(2) .lt. max_lon) .and. &
           (hdr(3) .gt. min_lat) .and. (hdr(3) .lt. max_lat) .and. &
           (hdr(4) .ge. -1.5) .and. (hdr(4) .lt. 1.5)) then 
           !(hdr(3) .gt. min_lat) .and. (hdr(3) .lt. max_lat)) then
           hdr(8) = iret
           num_msg_sb_00 = num_msg_sb_00+1
           data_hdr_00(num_msg_sb_00,1:8) = hdr(1:8)
           data_obs_00(num_msg_sb_00,1,:,:) = obs(:,:)
           data_obs_00(num_msg_sb_00,2,:,:) = qms(:,:)
           data_obs_00(num_msg_sb_00,3,:,:) = pco(:,:)
           data_obs_00(num_msg_sb_00,4,:,:) = oes(:,:)

       end if

       write(*,'(A,3I7,13f15.1)')'sb_report', num_subset,num_msg_sb_00,iret,hdr(2:7),obs(1:6,1)
     enddo sb_report_00
           
   enddo msg_report_00
 call closbf(unit_in)

!----------------------------------------------------------------------
! READ OBS2 :  OBS < -1.5 hr
 allocate(data_hdr_m3(num_msg_sb_m3,8))
 allocate(data_obs_m3(num_msg_sb_m3,4,8,255))
 data_hdr_m3 = 100000000000.0
 data_obs_m3 = 100000000000.0

 open(unit_in,file=trim(infilename),action='read',form='unformatted')
 call openbf(unit_in,'IN',unit_in)
 call datelen(10)
   num_message=0
   num_msg_sb_m3=0
   msg_report_m3: do while (ireadmg(unit_in,subset,idate) == 0)
     num_message=num_message+1
     num_subset = 0
     write(*,'(I10,I8,a10)') idate,num_message,subset
     sb_report_m3: do while (ireadsb(unit_in) == 0)
       num_subset = num_subset+1
       call ufbint(unit_in,hdr,7,1 ,iret,hdstr)
       call ufbint(unit_in,obs,8,255,iret,obstr)
       call ufbint(unit_in,qms,8,255,iret,qmstr)
       call ufbint(unit_in,pco,8,255,iret,pcstr)
       call ufbint(unit_in,oes,8,255,iret,oestr)

       if (hdr(2) >= 180.0) hdr(2) = hdr(2) - 360.0

       if ((trim(subset) .eq. trim(obs_typ)) .and. &
           (hdr(2) .gt. min_lon) .and. (hdr(2) .lt. max_lon) .and. &
           (hdr(3) .gt. min_lat) .and. (hdr(3) .lt. max_lat) .and. &
           (hdr(4) .lt. -1.5)) then
           !(hdr(3) .gt. min_lat) .and. (hdr(3) .lt. max_lat)) then
           hdr(8) = iret
           num_msg_sb_m3 = num_msg_sb_m3+1
           data_hdr_m3(num_msg_sb_m3,1:8) = hdr(1:8)
           data_obs_m3(num_msg_sb_m3,1,:,:) = obs(:,:)
           data_obs_m3(num_msg_sb_m3,2,:,:) = qms(:,:)
           data_obs_m3(num_msg_sb_m3,3,:,:) = pco(:,:)
           data_obs_m3(num_msg_sb_m3,4,:,:) = oes(:,:)
       end if

       write(*,'(A,3I7,13f15.1)')'sb_report', num_subset,num_msg_sb_m3,iret,hdr(2:7),obs(1:6,1)
     enddo sb_report_m3

   enddo msg_report_m3
 call closbf(unit_in)

!----------------------------------------------------------------------
! READ OBS3 : 1.5 hr <= OBS 
 allocate(data_hdr_p3(num_msg_sb_p3,8))
 allocate(data_obs_p3(num_msg_sb_p3,4,8,255))
 data_hdr_p3 = 100000000000.0
 data_obs_p3 = 100000000000.0

 open(unit_in,file=trim(infilename),action='read',form='unformatted')
 call openbf(unit_in,'IN',unit_in)
 call datelen(10)
   num_message=0
   num_msg_sb_p3=0
   msg_report_p3: do while (ireadmg(unit_in,subset,idate) == 0)
     num_message=num_message+1
     num_subset = 0
     write(*,'(I10,I8,a10)') idate,num_message,subset
     sb_report_p3: do while (ireadsb(unit_in) == 0)
       num_subset = num_subset+1
       call ufbint(unit_in,hdr,7,1 ,iret,hdstr)
       call ufbint(unit_in,obs,8,255,iret,obstr)
       call ufbint(unit_in,qms,8,255,iret,qmstr)
       call ufbint(unit_in,pco,8,255,iret,pcstr)
       call ufbint(unit_in,oes,8,255,iret,oestr)

       if (hdr(2) >= 180.0) hdr(2) = hdr(2) - 360.0

       if ((trim(subset) .eq. trim(obs_typ)) .and. &
           (hdr(2) .gt. min_lon) .and. (hdr(2) .lt. max_lon) .and. &
           (hdr(3) .gt. min_lat) .and. (hdr(3) .lt. max_lat) .and. &
           (hdr(4) .ge. 1.5)) then
           !(hdr(3) .gt. min_lat) .and. (hdr(3) .lt. max_lat)) then
           hdr(8) = iret
           num_msg_sb_p3 = num_msg_sb_p3+1
           data_hdr_p3(num_msg_sb_p3,1:8) = hdr(1:8)
           data_obs_p3(num_msg_sb_p3,1,:,:) = obs(:,:)
           data_obs_p3(num_msg_sb_p3,2,:,:) = qms(:,:)
           data_obs_p3(num_msg_sb_p3,3,:,:) = pco(:,:)
           data_obs_p3(num_msg_sb_p3,4,:,:) = oes(:,:)
       end if

       write(*,'(A,3I7,13f15.1)')'sb_report_p3', num_subset,num_msg_sb_p3,iret,hdr(2:7),obs(1:6,1)
     enddo sb_report_p3

   enddo msg_report_p3
 call closbf(unit_in)
 write(*,"(2I5,8f15.1)") num_message,num_msg_sb_00,data_hdr_00(num_msg_sb,:)
 write(*,"(2I5,8f15.1)") num_message,num_msg_sb_00,data_obs_00(2,1,:,1)

!---------------------------------------------------------------------
! write output
!---------------------------------------------------------------------
!  00
if (num_msg_sb_00 .ne. 0) then
   file_out="/glade/work/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_00)//"_obs.dat"
   open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_00*4*8*255)
   write(53,rec=1) data_obs_00
   close(53)

   file_out="/glade/work/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_00)//"_hdr.dat"
   open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_00*8)
   write(53,rec=1) data_hdr_00
   close(53)
end if

!  m03
if (num_msg_sb_m3 .ne. 0) then
   file_out="/glade/work/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_m3)//"_obs.dat"
   open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_m3*4*8*255)
   write(53,rec=1) data_obs_m3
   close(53)

   file_out="/glade/work/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_m3)//"_hdr.dat"
   open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_m3*8)
   write(53,rec=1) data_hdr_m3
   close(53)
end if

!  p03
if (num_msg_sb_p3 .ne. 0) then
   file_out="/glade/work/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_p3)//"_obs.dat"
   open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_p3*4*8*255)
   write(53,rec=1) data_obs_p3
   close(53)

   file_out="/glade/work/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_p3)//"_hdr.dat"
   open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_p3*8)
   write(53,rec=1) data_hdr_p3
   close(53)
end if
end program
