program prepbufr_encode_upperair
!
!  example of writing one value into a bufr file
!
 implicit none
 
 character(len=40)            :: obstr,hdstr,qmstr,oestr, pcstr
 
 character(8) subset
 integer :: unit_out=10,unit_table=20
 integer :: idate,iret

 integer, parameter :: mxmn=35, mxlv=1

!  namelist files
!
 character(180) :: infilename        ! file from GSI running directory
 character(180) :: obs_typ, obs_date
 character(10)  :: obs_date2
 character(200) :: input_dir, input_file, output_dir, output_file
 namelist/iosetup/ obs_typ, obs_date, input_dir, input_file, output_dir, output_file

 real*8             :: hdr(7)
 real*8             :: pmo(2,1)
 real*8             :: obs(8,255),qms(8,255),oes(8,255),pco(8,255)
 real               :: woe,toe,qoe,poe,pob,pwe
 real*8             :: obs_time
 integer              :: iyear, imonth, iday, ihour, imin

 real*8, allocatable :: data_obs(:,:,:,:), data_hdr(:,:)
 character(len=200) :: file_out
 
 integer :: num_msg_sb, nlvl, i4obs,nstn 

 hdstr='SID XOB YOB DHR TYP ELV T29'
 obstr='POB QOB TOB ZOB UOB VOB PWO CAT' ! observation
 qmstr='PQM QQM TQM ZQM WQM NUL PWQ NUL' ! quality marker
 oestr='POE QOE TOE NUL WOE NUL PWE NUL' ! observation error
 pcstr='PPC QPC TPC ZPC WPC NUL PWP NUL' ! program code

! Read namelist 

  open(11,file='namelist.mpd')
   read(11,iosetup)
  close(11)

! Read observation number saved by decoder
 open(10,file=trim(input_dir)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_obs_num.txt",status='old')
 read(10,'(I10)') num_msg_sb ! num_msg_sb
 close(10)
 print*, "obs_number = ", num_msg_sb

! Read OBS
 allocate(data_hdr(num_msg_sb,8))
 allocate(data_obs(num_msg_sb,3,8,255))

 data_hdr = 100000000000.0
 data_obs = 100000000000.0

!-----------------------------------------------------------
! V00 Read obs & head data
!-----------------------------------------------------------
if (num_msg_sb .ne. 0) then
 file_out=trim(input_dir)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_obs_"//trim(input_file)
 write(*,"(A,A)") "simulated MPD obs file=",trim(file_out)
 open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb*3*8*255)
 read(53,rec=1) data_obs
 close(53)

 file_out=trim(input_dir)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_hdr_"//trim(input_file)
 write(*,"(A,A)") "obs hdr file=",trim(file_out)
 open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb*8)
 read(53,rec=1) data_hdr
 close(53)

 !print*, data_hdr(10,:)
 print*, "data_obs"
 print*, data_obs(num_msg_sb,1,1,:)

 ! Open for writing
 open(unit_table,file='prepbufr.table',action='read')
 open(unit_out,file=trim(output_dir)//'/bufr_'//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(output_file),action='write',form='unformatted')

 subset = "SATEMP" !trim(obs_typ)
 obs_date2=trim(obs_date)
 Read(obs_date2, '(I10)' )  idate 
 print*, "idate= ", trim(obs_date), idate

 call datelen(10)
 call openbf(unit_out,'OUT',unit_table)
 call openmb(unit_out,subset,idate)

! observation loop
 DO i4obs=1, num_msg_sb
    hdr=10.0e10
    hdr(1:7) = data_hdr(i4obs,1:7)
    print*, hdr

    nlvl=26
    !obs=10.0e10;qcf=10.0e10;oer=10.0e10
    obs=10.0e10;qms=10.0e10;pco=10.0e10;oes=10.0e10
    obs(:,:) = data_obs(i4obs,1,:,:)
    qms(:,:) = data_obs(i4obs,2,:,:)
    oes(:,:) = data_obs(i4obs,3,:,:)

    !if ((idate .eq. 2015071412) .and. (hdr(4) .lt. 0.0)) then
    !   hdr(4) = 0.0
    !   write(*,"(A,i10)") "very first time, earier forecasts than 2015071412 are changed to the nature run valid at 2015071412 at ", i4obs
    !end if

    call ufbint(unit_out,hdr,7,1 ,iret,hdstr)
    call ufbint(unit_out,obs,8,nlvl,iret,obstr)
    call ufbint(unit_out,qms,8,nlvl,iret,qmstr)
    call ufbint(unit_out,pco,8,nlvl,iret,pcstr)
    call ufbint(unit_out,oes,8,nlvl,iret,oestr)
    call writsb(unit_out)
    print*, i4obs, " iret=", iret, obs(:,1)
 End Do

   call closmg(unit_out)
 call closbf(unit_out)
end if

!
end program
