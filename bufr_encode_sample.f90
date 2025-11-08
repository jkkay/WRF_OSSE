program bufr_encode_sample
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
 character(180) :: inputdir       ! file from GSI running directory
 character(180) :: outputdir       ! file from GSI running directory
 character(180) :: outfilename_00       ! file name saving results
 character(180) :: outfilename_m3       ! file name saving results
 character(180) :: outfilename_p3       ! file name saving results
 character(180) :: obs_typ, obs_date
 character(10)  :: obs_date2
 namelist/iosetup/ inputdir, outputdir, outfilename_00, outfilename_m3, outfilename_p3, obs_typ, obs_date

 real*8             :: hdr(8)
 real*8             :: pmo(2,1)
 real*8             :: obs(8,255),qms(8,255),oer(8,255),pco(8,255)
 real               :: woe,toe,qoe,poe,pob,pwe
 real*8             :: obs_time
 integer              :: iyear, imonth, iday, ihour, imin

 real*8, allocatable :: data_obs(:,:,:,:), data_hdr(:,:)
 real*8, allocatable :: data_obs_00(:,:,:,:), data_hdr_00(:,:)
 real*8, allocatable :: data_obs_m3(:,:,:,:), data_hdr_m3(:,:)
 real*8, allocatable :: data_obs_p3(:,:,:,:), data_hdr_p3(:,:)
 character(len=200) :: file_out
 
 integer :: num_msg_sb, nlvl, i4obs 
 integer :: num_msg_sb_m3, num_msg_sb_00, num_msg_sb_p3

 hdstr='SID XOB YOB DHR TYP ELV T29'
 obstr='POB QOB TOB ZOB UOB VOB PWO CAT' ! observation
 qmstr='PQM QQM TQM ZQM WQM NUL PWQ NUL' ! quality marker
 oestr='POE QOE TOE NUL WOE NUL PWE NUL' ! observation error
 pcstr='PPC QPC TPC ZPC WPC NUL PWP NUL' ! program code

! Read namelist 
  open(11,file='namelist.conv_out')
   read(11,iosetup)
  close(11)

! Read observation number saved by decoder
 !open(10,file=trim(obs_date)//"_"//trim(obs_typ)//"_obs_num.txt",status='old')
 open(10,file=trim(inputdir)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_obs_num.txt",status='old')
 read(10,'(3I10)') num_msg_sb_m3, num_msg_sb_00, num_msg_sb_p3 ! num_msg_sb
 close(10)
 print*, "obs_number = ", num_msg_sb_00

! Read OBS
 allocate(data_hdr_00(num_msg_sb_00,8))
 allocate(data_obs_00(num_msg_sb_00,4,8,255))
 allocate(data_hdr_m3(num_msg_sb_m3,8))
 allocate(data_obs_m3(num_msg_sb_m3,4,8,255))
 allocate(data_hdr_p3(num_msg_sb_p3,8))
 allocate(data_obs_p3(num_msg_sb_p3,4,8,255))

 data_hdr_00 = 100000000000.0
 data_obs_00 = 100000000000.0
 data_hdr_p3 = 100000000000.0
 data_obs_p3 = 100000000000.0
 data_hdr_m3 = 100000000000.0
 data_obs_m3 = 100000000000.0

!-----------------------------------------------------------
! V00 Read obs & head data
!-----------------------------------------------------------
if (num_msg_sb_00 .ne. 0) then
 file_out=trim(inputdir)//"/"//trim(obs_date)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_00)//"_sim_obs.dat"
 !file_out="/glade/u/home/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/"//trim(obs_date)//"_"//trim(obs_typ)//"_obs.dat"
 write(*,"(A,A)") "simulated obs file=",trim(file_out)
 open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_00*4*8*255)
 read(53,rec=1) data_obs_00
 close(53)

 file_out=trim(inputdir)//"/"//trim(obs_date)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_00)//"_sim_hdr.dat"
 !file_out=trim(inputdir)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_00)//"_hdr.dat"
 write(*,"(A,A)") "obs hdr file=",trim(file_out)
 open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_00*8)
 read(53,rec=1) data_hdr_00
 close(53)

 !print*, data_hdr_00(1,:)
 !print*, data_hdr_00(2,:)
 !print*, data_obs_00(num_msg_sb_00-1,1,3,:)
 ! Open for writing
 open(unit_table,file='prepbufr.table',action='read')
 open(unit_out,file=trim(outputdir)//"/bufr_"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_00),action='write',form='unformatted')

 subset = trim(obs_typ)
 !idate = ICHAR(trim(obs_date))
 obs_date2=trim(obs_date)
 Read(obs_date2, '(I10)' )  idate 
 print*, "idate= ", trim(obs_date), idate

 call datelen(10)
 call openbf(unit_out,'OUT',unit_table)
 call openmb(unit_out,subset,idate)

! observation loop
 DO i4obs=1, num_msg_sb_00
    hdr=10.0e10
    hdr(1:8) = data_hdr_00(i4obs,1:8)
    print*, hdr

    nlvl=IDINT(hdr(8)) 
    print*, "nlvl", hdr(8), nlvl
    !obs=10.0e10;qcf=10.0e10;oer=10.0e10
    obs=10.0e10;qms=10.0e10;pco=10.0e10;oer=10.0e10
    obs(:,:) = data_obs_00(i4obs,1,:,:)
    qms(:,:) = data_obs_00(i4obs,2,:,:)
    pco(:,:) = data_obs_00(i4obs,3,:,:)
    oer(:,:) = data_obs_00(i4obs,4,:,:)

    if ((idate .eq. 2015071412) .and. (hdr(4) .lt. 0.0)) then
       hdr(4) = 0.0
       write(*,"(A,i10)") "very first time, earier forecasts than 2015071412 are changed to the nature run valid at 2015071412 at ", i4obs
    end if

    !if ((hdr(5) .eq. 281) .and. (qms(5,1) .eq. 9) ) then
    if ((hdr(5) .eq. 281) .and. ((qms(5,1) .eq. 9) .or. (qms(5,1) .eq. 14))) then
       print*,"jk_281 ", i4obs,obs(1,1), oer(5,1), qms(5,1)
       oer(5,1) = 1.8854*2.5
       qms(5,1) = 2
    !else if ((hdr(5) .eq. 284) .and. (qms(5,1) .eq. 9)) then
    else if ((hdr(5) .eq. 284) .and. ((qms(5,1) .eq. 9) .or. (qms(5,1) .eq. 14))) then
       print*, "jk_284 ",i4obs, obs(1,1),oer(5,1), qms(5,1)
       oer(5,1) = 1.761*2.5
       qms(5,1) = 2
    else if ((hdr(5) .eq. 287) .and. ((qms(5,1) .eq. 9) .or. (qms(5,1) .eq. 14))) then
       print*, "jk_287", i4obs, obs(1,1),oer(5,1), qms(5,1)
       oer(5,1) = 1.9162*2.5
       qms(5,1) = 2
    end if

    !jkmod4
    if ((hdr(5) .eq. 281) .or. (hdr(5) .eq. 284) .or. (hdr(5) .eq. 287)) then
       if ((i4obs .lt. num_msg_sb_00) .and. (data_hdr_00(i4obs+1,2) .eq. hdr(2)) .and. (data_hdr_00(i4obs+1,3) .eq. hdr(3))) then
           print*, "jk_copyP_before", obs(1,1)
           obs(1,1) = data_obs_00(i4obs+1,1,1,1)
           print*, "jk_copyP_after", obs(1,1)
       end if
    end if
    !jkmod4

    if ((hdr(5) .eq. 181) .and. ((qms(2,1) .eq. 9) .or. (qms(2,1) .eq. 14))) then
       print*,"jk_181 ", i4obs, obs(2,1), oer(2,1), qms(2,1)
       oer(1,1) = oer(1,1)*2.5
       oer(2,1) = 0.39136*2.5 
       qms(2,1) = 2
    else if ((hdr(5) .eq. 184) .and. ((qms(2,1) .eq. 9) .or. (qms(2,1) .eq. 14))) then
       print*, "jk_184 ",i4obs, obs(2,1), oer(2,1), qms(2,1)
       oer(1,1) = oer(1,1)*2.5
       oer(2,1) = 0.59120*2.5 
       qms(2,1) = 2
    else if ((hdr(5) .eq. 187) .and. ((qms(2,1) .eq. 9) .or. (qms(2,1) .eq. 14))) then
       print*, "jk_187", i4obs, obs(2,1), oer(2,1), qms(2,1)
       oer(1,1) = oer(1,1)*2.5
       oer(2,1) = 0.27059*2.5 
       qms(2,1) = 2
    end if

    if ((hdr(5) .eq. 181) .and. ((qms(3,1) .eq. 9) .or. (qms(3,1) .eq. 14))) then
       print*,"jk_t181 ", i4obs, obs(3,1), oer(3,1), qms(3,1)
       oer(3,1) = oer(3,1)*2.5 
       qms(3,1) = 2
    else if ((hdr(5) .eq. 184) .and. ((qms(3,1) .eq. 9) .or. (qms(3,1) .eq. 14))) then
       print*, "jk_t184 ",i4obs, obs(3,1), oer(3,1), qms(3,1)
       oer(3,1) = oer(3,1)*2.5 
       qms(3,1) = 2
    else if ((hdr(5) .eq. 187) .and. ((qms(3,1) .eq. 9) .or. (qms(3,1) .eq. 14))) then
       print*, "jk_t187", i4obs, obs(3,1), oer(3,1), qms(3,1)
       oer(3,1) = oer(3,1)*2.5 
       qms(3,1) = 2
    end if

! RawinUVTP, no Q ...jkmod1
    if (hdr(5) .eq. 120) then
       qms(2,:) = 14
    end if
! RawinUVTP, no Q ...jkmod1

    call ufbint(unit_out,hdr,7,1 ,iret,hdstr)
    call ufbint(unit_out,obs,8,nlvl,iret,obstr)
    call ufbint(unit_out,qms,8,nlvl,iret,qmstr)
    call ufbint(unit_out,pco,8,nlvl,iret,pcstr)
    call ufbint(unit_out,oer,8,nlvl,iret,oestr)
    call writsb(unit_out)
    print*, i4obs, " iret=", iret
 End Do

   call closmg(unit_out)
 call closbf(unit_out)
end if
!-----------------------------------------------------------
! m3 Read obs & head data
!-----------------------------------------------------------
if (num_msg_sb_m3 .ne. 0) then
 file_out=trim(inputdir)//"/"//trim(obs_date)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_m3)//"_sim_obs.dat"
 !file_out="/glade/u/home/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/"//trim(obs_date)//"_"//trim(obs_typ)//"_obs.dat"
 write(*,"(A,A)") "simulated obs file=",trim(file_out)
 open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_m3*4*8*255)
 read(53,rec=1) data_obs_m3
 close(53)

 file_out=trim(inputdir)//"/"//trim(obs_date)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_m3)//"_sim_hdr.dat"
 !file_out=trim(inputdir)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_m3)//"_hdr.dat"
 write(*,"(A,A)") "obs hdr file=",trim(file_out)
 open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_m3*8)
 read(53,rec=1) data_hdr_m3
 close(53)

 print*, data_hdr_m3(1,:)
 print*, data_hdr_m3(2,:)
 print*, data_obs_m3(num_msg_sb_m3-1,1,3,:)
 ! Open for writing
 open(unit_out,file=trim(outputdir)//"/bufr_"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_m3),action='write',form='unformatted')

 call datelen(10)
 call openbf(unit_out,'OUT',unit_table)
 call openmb(unit_out,subset,idate)

! observation loop
 DO i4obs=1, num_msg_sb_m3
    hdr=10.0e10
    hdr(1:8) = data_hdr_m3(i4obs,1:8)
    print*, hdr

    nlvl=IDINT(hdr(8))
    print*, "nlvl", hdr(8), nlvl
    !obs=10.0e10;qcf=10.0e10;oer=10.0e10
    obs=10.0e10;qms=10.0e10;pco=10.0e10;oer=10.0e10
    obs(:,:) = data_obs_m3(i4obs,1,:,:)
    qms(:,:) = data_obs_m3(i4obs,2,:,:)
    pco(:,:) = data_obs_m3(i4obs,3,:,:)
    oer(:,:) = data_obs_m3(i4obs,4,:,:)

    if ((idate .eq. 2015071412) .and. (hdr(4) .lt. 0.0)) then
       hdr(4) = 0.0
       write(*,"(A,i10)") "very first time, earier forecasts than 2015071412 are changed to the nature run valid at 2015071412 at ", i4obs
    end if

    !if ((hdr(5) .eq. 281) .and. (qms(5,1) .eq. 9)) then
    if ((hdr(5) .eq. 281) .and. ((qms(5,1) .eq. 9) .or. (qms(5,1) .eq. 14))) then
       print*,"jk_281 ", i4obs, oer(5,1), qms(5,1)
       oer(5,1) = 1.8854*2.5
       qms(5,1) = 2
    !else if ((hdr(5) .eq. 284) .and. (qms(5,1) .eq. 9)) then
    else if ((hdr(5) .eq. 284) .and. ((qms(5,1) .eq. 9) .or. (qms(5,1) .eq. 14))) then
       print*, "jk_284 ",i4obs, oer(5,1), qms(5,1)
       oer(5,1) = 1.761*2.5
       qms(5,1) = 2
    !else if ((hdr(5) .eq. 287) .and. (qms(5,1) .eq. 9)) then
    else if ((hdr(5) .eq. 287) .and. ((qms(5,1) .eq. 9) .or. (qms(5,1) .eq. 14))) then
       print*, "jk_287", i4obs, oer(5,1), qms(5,1)
       oer(5,1) = 1.9162*2.5
       qms(5,1) = 2
    end if
    !jkmod4
    if ((hdr(5) .eq. 281) .or. (hdr(5) .eq. 284) .or. (hdr(5) .eq. 287)) then
       if ((i4obs .lt. num_msg_sb_m3) .and. (data_hdr_m3(i4obs+1,2) .eq. hdr(2)) .and. (data_hdr_m3(i4obs+1,3) .eq. hdr(3))) then
           print*, "jk_copyP_before", obs(1,1)
           obs(1,1) = data_obs_m3(i4obs+1,1,1,1)
           print*, "jk_copyP_after", obs(1,1)
       end if
    end if
    !jkmod4


   if ((hdr(5) .eq. 181) .and. ((qms(2,1) .eq. 9) .or. (qms(2,1) .eq. 14))) then
       print*,"jk_181 ", i4obs, obs(2,1), oer(2,1), qms(2,1)
       oer(1,1) = oer(1,1)*2.5
       oer(2,1) = 0.39136*2.5
       qms(2,1) = 2
    else if ((hdr(5) .eq. 184) .and. ((qms(2,1) .eq. 9) .or. (qms(2,1) .eq. 14))) then
       print*, "jk_184 ",i4obs, obs(2,1), oer(2,1), qms(2,1)
       oer(1,1) = oer(1,1)*2.5
       oer(2,1) = 0.59120*2.5
       qms(2,1) = 2
    else if ((hdr(5) .eq. 187) .and. ((qms(2,1) .eq. 9) .or. (qms(2,1) .eq. 14))) then
       print*, "jk_187", i4obs, obs(2,1), oer(2,1), qms(2,1)
       oer(1,1) = oer(1,1)*2.5
       oer(2,1) = 0.27059*2.5
       qms(2,1) = 2
    end if

    if ((hdr(5) .eq. 181) .and. ((qms(3,1) .eq. 9) .or. (qms(3,1) .eq. 14))) then
       print*,"jk_t181 ", i4obs, obs(3,1), oer(3,1), qms(3,1)
       oer(3,1) = oer(3,1)*2.5
       qms(3,1) = 2
    else if ((hdr(5) .eq. 184) .and. ((qms(3,1) .eq. 9) .or. (qms(3,1) .eq. 14))) then
       print*, "jk_t184 ",i4obs, obs(3,1), oer(3,1), qms(3,1)
       oer(3,1) = oer(3,1)*2.5
       qms(3,1) = 2
    else if ((hdr(5) .eq. 187) .and. ((qms(3,1) .eq. 9) .or. (qms(3,1) .eq. 14))) then
       print*, "jk_t187", i4obs, obs(3,1), oer(3,1), qms(3,1)
       oer(3,1) = oer(3,1)*2.5
       qms(3,1) = 2
    end if

! RawinUVTP, no Q ...jkmod1
    if (hdr(5) .eq. 120) then
       qms(2,:) = 14
    end if
! RawinUVTP, no Q ...jkmod1

    call ufbint(unit_out,hdr,7,1 ,iret,hdstr)
    call ufbint(unit_out,obs,8,nlvl,iret,obstr)
    call ufbint(unit_out,qms,8,nlvl,iret,qmstr)
    call ufbint(unit_out,pco,8,nlvl,iret,pcstr)
    call ufbint(unit_out,oer,8,nlvl,iret,oestr)
    call writsb(unit_out)
    print*, i4obs, " iret=", iret
 End Do

   call closmg(unit_out)
 call closbf(unit_out)
end if

!-----------------------------------------------------------
! p3 Read obs & head data
!-----------------------------------------------------------
if (num_msg_sb_p3 .ne. 0) then
 file_out=trim(inputdir)//"/"//trim(obs_date)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_p3)//"_sim_obs.dat"
 !file_out="/glade/u/home/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/"//trim(obs_date)//"_"//trim(obs_typ)//"_obs.dat"
 write(*,"(A,A)") "simulated obs file=",trim(file_out)
 open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_p3*4*8*255)
 read(53,rec=1) data_obs_p3
 close(53)

 !jkmod3 uses changed hdr file that include modified elevation from the nature run
 file_out=trim(inputdir)//"/"//trim(obs_date)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_p3)//"_sim_hdr.dat"  
 !file_out=trim(inputdir)//"/"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_p3)//"_hdr.dat"
 write(*,"(A,A)") "obs hdr file=",trim(file_out)
 open(53,file=trim(file_out),form='unformatted',access='direct',recl=8*num_msg_sb_p3*8)
 read(53,rec=1) data_hdr_p3
 close(53)

 print*, data_hdr_p3(1,:)
 print*, data_hdr_p3(2,:)
 print*, data_obs_p3(num_msg_sb_p3-1,1,3,:)
 ! Open for writing
 open(unit_out,file=trim(outputdir)//"/bufr_"//trim(obs_date)//"_"//trim(obs_typ)//"_"//trim(outfilename_p3),action='write',form='unformatted')

 call datelen(10)
 call openbf(unit_out,'OUT',unit_table)
 call openmb(unit_out,subset,idate)

! observation loop
 DO i4obs=1, num_msg_sb_p3
    hdr=10.0e10
    hdr(1:8) = data_hdr_p3(i4obs,1:8)
    print*, hdr

    nlvl=IDINT(hdr(8))
    print*, "nlvl", hdr(8), nlvl
    !obs=10.0e10;qcf=10.0e10;oer=10.0e10
    obs=10.0e10;qms=10.0e10;pco=10.0e10;oer=10.0e10
    obs(:,:) = data_obs_p3(i4obs,1,:,:)
    qms(:,:) = data_obs_p3(i4obs,2,:,:)
    pco(:,:) = data_obs_p3(i4obs,3,:,:)
    oer(:,:) = data_obs_p3(i4obs,4,:,:)

    if ((hdr(5) .eq. 281) .and. ((qms(5,1) .eq. 9) .or. (qms(5,1) .eq. 14))) then
    !if ((hdr(5) .eq. 281) .and. (qms(5,1) .eq. 9)) then
       print*,"jk_281 ", i4obs, oer(5,1), qms(5,1)
       oer(5,1) = 1.8854*2.5
       qms(5,1) = 2
    !else if ((hdr(5) .eq. 284) .and. (qms(5,1) .eq. 9)) then
    else if ((hdr(5) .eq. 284) .and. ((qms(5,1) .eq. 9) .or. (qms(5,1) .eq. 14))) then
       print*, "jk_284 ",i4obs, oer(5,1), qms(5,1)
       oer(5,1) = 1.761*2.5
       qms(5,1) = 2
    !else if ((hdr(5) .eq. 287) .and. (qms(5,1) .eq. 9)) then
    else if ((hdr(5) .eq. 287) .and. ((qms(5,1) .eq. 9) .or. (qms(5,1) .eq. 14))) then
       print*, "jk_287", i4obs, oer(5,1), qms(5,1)
       oer(5,1) = 1.9162*2.5
       qms(5,1) = 2
    end if
    !jkmod4
    if ((hdr(5) .eq. 281) .or. (hdr(5) .eq. 284) .or. (hdr(5) .eq. 287)) then
       if ((i4obs .lt. num_msg_sb_p3) .and. (data_hdr_p3(i4obs+1,2) .eq. hdr(2)) .and. (data_hdr_p3(i4obs+1,3) .eq. hdr(3))) then
           print*, "jk_copyP_before", obs(1,1)
           obs(1,1) = data_obs_p3(i4obs+1,1,1,1)
           print*, "jk_copyP_after", obs(1,1)
       end if
    end if
    !jkmod4


   if ((hdr(5) .eq. 181) .and. ((qms(2,1) .eq. 9) .or. (qms(2,1) .eq. 14))) then
       print*,"jk_181 ", i4obs, obs(2,1), oer(2,1), qms(2,1)
       oer(1,1) = oer(1,1)*2.5
       oer(2,1) = 0.39136*2.5
       qms(2,1) = 2
    else if ((hdr(5) .eq. 184) .and. ((qms(2,1) .eq. 9) .or. (qms(2,1) .eq. 14))) then
       print*, "jk_184 ",i4obs, obs(2,1), oer(2,1), qms(2,1)
       oer(1,1) = oer(1,1)*2.5
       oer(2,1) = 0.59120*2.5
       qms(2,1) = 2
    else if ((hdr(5) .eq. 187) .and. ((qms(2,1) .eq. 9) .or. (qms(2,1) .eq. 14))) then
       print*, "jk_187", i4obs, obs(2,1), oer(2,1), qms(2,1)
       oer(1,1) = oer(1,1)*2.5
       oer(2,1) = 0.27059*2.5
       qms(2,1) = 2
    end if

    if ((hdr(5) .eq. 181) .and. ((qms(3,1) .eq. 9) .or. (qms(3,1) .eq. 14))) then
       print*,"jk_t181 ", i4obs, obs(3,1), oer(3,1), qms(3,1)
       oer(3,1) = oer(3,1)*2.5
       qms(3,1) = 2
    else if ((hdr(5) .eq. 184) .and. ((qms(3,1) .eq. 9) .or. (qms(3,1) .eq. 14))) then
       print*, "jk_t184 ",i4obs, obs(3,1), oer(3,1), qms(3,1)
       oer(3,1) = oer(3,1)*2.5
       qms(3,1) = 2
    else if ((hdr(5) .eq. 187) .and. ((qms(3,1) .eq. 9) .or. (qms(3,1) .eq. 14))) then
       print*, "jk_t187", i4obs, obs(3,1), oer(3,1), qms(3,1)
       oer(3,1) = oer(3,1)*2.5
       qms(3,1) = 2
    end if

! RawinUVTP, no Q ...jkmod1
    if (hdr(5) .eq. 120) then
       qms(2,:) = 14
    end if
! RawinUVTP, no Q ...jkmod1

    call ufbint(unit_out,hdr,7,1 ,iret,hdstr)
    call ufbint(unit_out,obs,8,nlvl,iret,obstr)
    call ufbint(unit_out,qms,8,nlvl,iret,qmstr)
    call ufbint(unit_out,pco,8,nlvl,iret,pcstr)
    call ufbint(unit_out,oer,8,nlvl,iret,oestr)
    call writsb(unit_out)
    print*, i4obs, " iret=", iret
 End Do

   call closmg(unit_out)
 call closbf(unit_out)
end if
!
end program
