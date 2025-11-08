program prepbufr_decode_all
!
! read all observations out from prepbufr. 
! read bufr table from prepbufr file
!
 implicit none

 integer, parameter :: mxmn=35, mxlv=250
 character(80):: hdstr='SID XOB YOB DHR TYP ELV SAID T29'
 character(80):: obstr='POB QOB TOB ZOB UOB VOB PWO CAT PRSS'
 character(80):: qcstr='PQM QQM TQM ZQM WQM NUL PWQ     '
 character(80):: oestr='POE QOE TOE NUL WOE NUL PWE     '
 real(8) :: hdr(mxmn),obs(mxmn,mxlv),qcf(mxmn,mxlv),oer(mxmn,mxlv)

 INTEGER        :: ireadmg,ireadsb

 character(8)   :: subset
 integer        :: unit_in=10,idate,nmsg,ntb

 character(8)   :: c_sid
 real(8)        :: rstation_id
 equivalence(rstation_id,c_sid)

 integer        :: i,k,iret
!
!
 open(24,file='prepbufr.table')
 !open(unit_in,file='prepbufr',form='unformatted',status='old')
 !open(unit_in,file='/glade/u/home/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/nr_rawin/bufr_2015071500_MPD_R2.0',form='unformatted',status='old')
 open(unit_in,file='/glade/work/junkyung/data/obs/LIUZ/20150714/gdas1.t12z.prepbufr.nr',form='unformatted',status='old')
 !open(unit_in,file='/glade/work/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/using_GSI_obsR_mod3/sim_bufr/bufr_2015071412_2015071506_conv',form='unformatted',status='old')
 call openbf(unit_in,'IN',unit_in)
 call dxdump(unit_in,24)
 call datelen(10)
   nmsg=0
   msg_report: do while (ireadmg(unit_in,subset,idate) == 0)
     nmsg=nmsg+1
     ntb = 0
     write(*,*)
     !write(*,'(3a,i10)') 'subset=',subset,' cycle time =',idate
     sb_report: do while (ireadsb(unit_in) == 0)
       ntb = ntb+1
       call ufbint(unit_in,hdr,mxmn,1   ,iret,hdstr)
       call ufbint(unit_in,obs,mxmn,mxlv,iret,obstr)
       call ufbint(unit_in,oer,mxmn,mxlv,iret,oestr)
       call ufbint(unit_in,qcf,mxmn,mxlv,iret,qcstr)
       rstation_id=hdr(1)
       write(*,*)

       if (hdr(2) >= 180.0) hdr(2) = hdr(2) - 360.0 
       !if( (hdr(8) .eq. 511) .or. (hdr(8) .eq. 514)) then 
       !if( (hdr(5) .eq. 126) ) then 
       !   write(*,'(2I10,a14,8f14.1)') ntb,iret,c_sid,(hdr(i),i=2,8)
       !   DO k=1,iret
       !  !write(*,'(i3,a10,3f14.1)') k,'obs=',(obs(i,k),i=1,3)
       !       write(13,'(2f14.4,f14.1,a5,f14.1,a5,f14.1,a5,2f14.1)') hdr(2), hdr(3), obs(1,k),' hPa ', obs(2,k),' w/o ',oer(2,k),' err ', qcf(2,k), obs(4,k)
       !  !write(*,'(i3,a10,3f14.1)') k,'qcf=',(qcf(i,k),i=1,3)
       !   ENDDO
       !else if( (hdr(5) .eq. 223) .or. (hdr(5) .eq. 227) .or. (hdr(5) .eq. 229)) then 
          !DO k=1,iret
          !    write(14,'(i10,3f14.4,f14.1,a5,f14.1,a5,f14.1,a5,f14.1,a5,f14.1,a5,f14.1)') k, hdr(2), hdr(3), hdr(5), obs(1,k),' hPa ', obs(4,k),'Z m ',obs(5,k), 'U m/s',oer(5,k),' err ', qcf(5,k)
          !ENDDO       
       !end if
       if (subset .eq. 'ADPSFC') then
          write(*,'(8f14.4)') hdr(1:8)
          write(*,'(9f14.4)') obs(1:9,1) 
          write(*,'(9f14.4)') qcf(1:9,1) 
          write(*,'(9f14.4)') oer(1:9,1) 
       end if
     enddo sb_report
   enddo msg_report
 call closbf(unit_in)

end program
