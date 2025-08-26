#!/bin/ksDh

export WORK_DIR=/glade/work/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/knu
export OBS_DATES="2022011100" #"2015071500" #  2015071506"
export OBS_TYPES="ADPUPA" # PROFLR"
#export OBS_TYPES="ADPUPA AIRCFT PROFLR SATWND VADWND"
export OBS_DHRS="00" # m3 p3"
#export OBS_DHRS="00 m3 p3"

for OBS_DATE in $OBS_DATES; do
mkdir -p $WORK_DIR/$OBS_DATE
cd $WORK_DIR/$OBS_DATE

for OBS_TYPE in $OBS_TYPES; do
for OBS_DHR in $OBS_DHRS; do

rm -f $WORK_DIR/$OBS_DATE/mk_${OBS_TYPE}_${OBS_DATE}_nrobs_${OBS_DHR}.ncl
cat > mk_${OBS_TYPE}_${OBS_DATE}_nrobs_${OBS_DHR}.ncl << EOF
load "\$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "\$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

undef("da_get_q_error")
function da_get_q_error(RH, T_K, RH0, T_C0, P_MB0)
begin

    ;P_MB0 = P_PA0 / 100.

    ES  = 6.112 * 17.67 * 243.5 * T_K / ((T_C0+243.5)*(T_C0+243.5)) * exp (17.67*(T_C0)/(T_C0+243.5))
    ES0 = 6.112 * exp (17.67*(T_C0)/(T_C0+243.5))

    QS  = 0.622 * (P_MB0 * ES ) /  \
                 ((P_MB0 - ES0) * (P_MB0 - ES0))
    QS0 = 0.622 * ES0 /(P_MB0-ES0)

    QV_KG  = 0.01 * (RH0 * QS + RH * QS0)
    QV_KG0 = 0.01 * RH0 * QS0

    if (QV_KG0 .lt. 0.0) then
         QV_KG = 0.
    end if

return(QV_KG)

end

begin

  date="${OBS_DATE}"
  obs_typ="${OBS_TYPE}"; 
  path_hdr = "${WORK_DIR}/"
  file_hdr = path_hdr+date+"_"+obs_typ+"_V${OBS_DHR}_gdas1_hdr.dat"
  file_obs = path_hdr+date+"_"+obs_typ+"_V${OBS_DHR}_gdas1_obs.dat"

  R_factor=1.0

  wrf_dir="/glade/scratch/junkyung/KNU/ICBC/rundir/test/era5/usphy/wrfrun/fc/2022011012.e042/"
  print("file_hdr= "+file_hdr)

  FILESIZE=systemfunc("wc -c "+file_hdr+" | awk '{print \$1}'")
  num_msg_sb = stringtointeger(FILESIZE)/64
  print("number of observation = "+num_msg_sb)

  ; Read observation binary files
  setfileoption("bin","ReadByteOrder","BigEndian")
  nrec = 0
  dims_hdr = (/8,num_msg_sb/)
  dims_obs = (/255,8,4,num_msg_sb/) ; (lev,var,obs-qms-pco,number)
  hdr = fbindirread(file_hdr, nrec, dims_hdr, "double")
  obs = fbindirread(file_obs, nrec, dims_obs, "double")

  hdr_model = hdr
  obs_model = obs

  print(hdr(:,0))
  print(hdr(:,num_msg_sb-1))
;----------------------------------------------------------
; Non-regular Sonde observation location
;----------------------------------------------------------
lon_asonde=(/ -97.40, -99.00, -99.30, -99.60, -101.40/)
lat_asonde=(/  38.10,  40.50,  37.60,  38.90,   39.40/)
nasonde=dimsizes(lon_asonde)

;----------------------------------------------------------
; Check Time of the observation within window
;----------------------------------------------------------
  time_00 = new(1,integer)
  time_m3 = new(1,integer)
  time_p3 = new(1,integer)
  time_00 = 0
  time_m3 = 0
  time_p3 = 0

  do i4o=0, num_msg_sb-1
     if (hdr(3,i4o) .lt. -1.5) then    ; -3.0 hr ~ -1.5 hr
        time_m3 = time_m3 + 1
     else if ((hdr(3,i4o) .ge. -1.5) .and. (hdr(3,i4o) .lt. 1.5)) then  ; -1.5 hr ~ 1.5 hr
        time_00 = time_00 + 1
     else if (hdr(3,i4o) .ge. 1.5) then    ; 1.5 hr ~ 3.0 hr
        time_p3 = time_p3 + 1
     end if
     end if
     end if
  end do

  print("T+00hr obs number = "+time_00)
  print("T-03hr obs number = "+time_m3)
  print("T+03hr obs number = "+time_p3)

;----------------------------------------------------------
; Read Observation error
;----------------------------------------------------------
; 1) 120 ; T, Q
  Rname = "${WORK_DIR}/120_observation_R.txt"  ;jkmod
  Rdata  = asciiread(Rname,-1,"string")
  Rlev  = stringtofloat(str_get_cols(Rdata,0,12))
  np = dimsizes(Rlev)

  delete(Rname)
  delete(Rdata)

  Rname = "${WORK_DIR}/nam_errtable.r3dv"
  Rdata = asciiread(Rname,-1,"string")
  print(Rdata)

  typ0 = str_get_cols(Rdata, 0,3)
  typ = stringtofloat(typ0)
  print(typ)
;----------------------------------------------------------
  print("      LAT       LON       NLEV       POB         ZOB        TYP        T29")
  do i4o=0, num_msg_sb-1
    print (sprintf("%9.2f", hdr(1,i4o))    +" " \
          +sprintf("%9.2f", hdr(2,i4o)) +"  " \
          +sprintf("%9.2f", hdr(7,i4o)) +"  " \
          +sprintf("%9.2f", obs(0,0,0,i4o)) +"  " \
          +sprintf("%9.2f", obs(0,3,0,i4o)) +"  " \
          +sprintf("%9.2f", hdr(4,i4o))+"  " \
          +sprintf("%9.2f", hdr(6,i4o))    )
  end do

  ;print(obs(0,:,0,0))
  ;print(obs(0,:,1,0))
;----------------------------------------------------------
; Read nature run according to the observation time
;----------------------------------------------------------

  do i4o=0 , num_msg_sb-1
     ;dhr_min=hdr(3,i4o)*60.
     dhr_min=round(hdr(3,i4o)*60.,0)
     obs_time=systemfunc("/glade/u/home/junkyung/WRF/V4.1.2/WRFDA/var/build/da_advance_time.exe " + date + " "+dhr_min+"m")
     print("input obs time : "+hdr(3,i4o)+" "+dhr_min+" "+obs_time)

     yyyy = str_get_cols(obs_time, 0, 3)
     mm = str_get_cols(obs_time, 4, 5)
     dd = str_get_cols(obs_time, 6, 7)
     hh = str_get_cols(obs_time, 8, 9)

     mmm = new(1,string)

     ;print(strlen(obs_time))

     if (strlen(obs_time) .eq. 10) then
         mmm="00"
     else
         mmm=str_get_cols(obs_time,10,11)
     end if

     print("final obs time is "+yyyy+" "+mm+" "+dd+" "+hh+" "+mmm)

     ;-------------------------------------------------------------
     ; READ WRF output
     ;-------------------------------------------------------------
     in_file=wrf_dir+"/wrfout_d01_"+yyyy+"-"+mm+"-"+dd+"_"+hh+":"+mmm+":00"
     print("input wrf file = "+in_file)

     ; For the very first cycling valid at 1200 UTC 14 July
     if (date .eq. "2022011012") then
        if (dhr_min .lt. 0.0) then
           in_file=wrf_dir+"/wrfout_d01_2022-01-10_12:00:00"
        end if
     end if

     f = addfile(in_file+".nc","r")

     z = wrf_user_getvar(f,"z",0)        ; full model height [m]
     p = wrf_user_getvar(f,"pressure",0) ; full model pressure [hPa]

     if (i4o .eq. 0) then
        lat_wrf = wrf_user_getvar(f,"lat",0)
        lon_wrf = wrf_user_getvar(f,"lon",0)

        nj = dimsizes(z(0,:,0)) ; 320
        ni = dimsizes(z(0,0,:)) ; 410

        print(ni+" "+nj)
     end if

     flon = hdr(1,i4o)
     flat = hdr(2,i4o)
     fnz = doubletoint(hdr(7,i4o))
     fpob = obs(0:fnz-1,0,0,i4o)
     fzob = obs(0:fnz-1,3,0,i4o)

     print(i4o+" th nz is "+fnz)

     obsij = wrf_latlon_to_ij(lat_wrf,lon_wrf,flat,flon)
     obsi = obsij(0)      ; lat
     obsj = obsij(1)      ; lon

     zin = z(:,obsi,obsj)
     pin = p(:,obsi,obsj)

  ;-------------------------------------------------------------
  ; MASS observation
  ;-------------------------------------------------------------
  if (hdr(4,i4o) .lt. 200.) then
     t = wrf_user_getvar(f,"tc",0)        ; Temp [C]
     tv= wrf_user_getvar(f,"tv",0)        ; Virtual Temp [K]
     tv=tv-273.15                         ; Virtual Temp [c]

     qv= wrf_user_getvar(f,"QVAPOR",0)    ; QV [kg/kg]
     q = qv
     q = qv/(qv+1)                        ; specific humidity [kg/kg]
     q = q*1000000.                       ; specific humidity [mg/kg]

     rh= wrf_user_getvar(f,"rh",0)        ; relative humidity [%]

     qin  = q(:,obsi,obsj) ; level, lat, lon
     tvin = tv(:,obsi,obsj)
     ;zin  = z(:,obsi,obsj)
     rhin = rh(:,obsi,obsj)

     q_d  = wrf_interp_1d(qin,pin,fpob)
     tv_d = wrf_interp_1d(tvin,pin,fpob)
     z_d  = wrf_interp_1d(zin,pin,fpob)

     p_sfc = flt2dble(pin(0))  ; jkmod3 
     q_sfc = flt2dble(qin(0))
     tv_sfc = flt2dble(tvin(0))
     z_sfc = flt2dble(zin(0))

     print("jk_ij="+obsi+" "+obsj)
     print("jk_lon_lat="+flon+" "+flat)
     print("jk_zin="+zin(0))
     print("jk_pin="+pin(0))
     print("jkobs_p="+fpob(0))
     print("jkobs_z="+fzob(0))

     do i4z=0, fnz-1
        print(sprintf("%8.2f",obs(i4z,0,0,i4o))+" hPa, MODEL Q: "+sprintf("%15.2f", q_d(i4z))+" / OBS : "+sprintf("%15.2f",obs(i4z,1,0,i4o)))
     end do

    ; READ OBSERVATION ERROR
     ;Rname = hdr(4,i4o)+"_observation_R.txt"  ;jkmod
     ;Rdata  = asciiread(Rname,-1,"string")
     ;R_P  = stringtofloat(str_get_cols(Rdata, 0,12))
     ;R_T  = stringtofloat(str_get_cols(Rdata,13,24))
     ;R_RH = stringtofloat(str_get_cols(Rdata,25,36))    ; RH error, not SH error at this time
     ;print("RH_Rerr="+R_RH)

     i = 0
     do while(typ(i) .ne. hdr(4,i4o))
       i=i+1
     end do
     print(hdr(4,i4o)+" is "+i)
     P0  = stringtofloat(str_get_cols(Rdata(i+1:i+33), 0,12))    ; [mb]
     R_P  = stringtofloat(str_get_cols(Rdata(i+1:i+33),49,60))   ; [mb]
     R_T  = stringtofloat(str_get_cols(Rdata(i+1:i+33),13,24))
     R_RH = stringtofloat(str_get_cols(Rdata(i+1:i+33),25,36))   ; [%]/10
     R_RH = R_RH * 10.                                           ; [%]

    ; CONVERT rh observation error (R_RH) to specific humidity observation error (R_Q) (da_get_q_error.inc)
     tv_r = wrf_interp_1d(tvin,pin,P0)
     rh_r = wrf_interp_1d(rhin,pin,P0)
     q_r  = wrf_interp_1d(qin,pin,P0)

     R_Q = new((/np/),float,0.10000E+10)
     print("P_level        R_RH        R_T        rh         tv         q        R_Q ")

     do i4p=0, np-1
        if ((ismissing(tv_r(i4p)) .eq. False) .and. (ismissing(rh_r(i4p)) .eq. False)) then
            ; da_get_q_error : R_RH [%], R_T[K], rh_r[%], tv_r[C], P0[mb] in f_qv_from_rh.f
            R_Q(i4p) = da_get_q_error(R_RH(i4p), R_T(i4p), rh_r(i4p), tv_r(i4p), P0(i4p))  ; (f_qv_from_rh)
            R_Q(i4p) = R_Q(i4p)*1000000.  ; [kg/kg] --> [mg/kg]
        end if

        print(sprintf("%10.6f",P0(i4p))+" "+sprintf("%10.6f",R_RH(i4p))+" "+sprintf("%10.6f",R_T(i4p))+" "+ \
              sprintf("%10.6f",rh_r(i4p))+" "+sprintf("%10.6f",tv_r(i4p))+" "+sprintf("%10.f",q_r(i4p))+" "+sprintf("%10.6f",R_Q(i4p)))
     end do

     ; normal distribution random noise (av = average, sd = standard deviation)
     ; random_setallseed(36484749, 9494848)               ; Set seeds (suggested, NOT required)
     rand1 = toint(systemfunc(" date +%s"))       ; jkmod7
     rand2 = toint(systemfunc(" date +%s"))       ; jkmod7
     random_setallseed(rand1,rand2)               ; Set seeds (suggested, NOT required)
     mi=generate_sample_indices( np, 0 )         ;jkmod7

     obserr_T = new((/np,np/),float)
     obserr_Q = new((/np,np/),float)
     obserr_P = new((/np,np/),float)    ; jkmod3
     obserr_T = 0.0
     obserr_Q = 0.0
     obserr_P = 0.0                      ; jkmod3
     do i4p=0,np-1
        obserr_P(:,i4p) = random_normal(0, (R_P(i4p)) , (/np/))  ;jkmod3   ; jkmod7
        obserr_T(:,i4p) = random_normal(0, (R_T(i4p)) , (/np/))  ;               ; jkmod7
        ;obserr_Q(:,i4p) = random_normal(0, 0.5*(R_Q(i4p)) , (/np/))  ; [mg/kg]       ; jkmod7
        obserr_Q(:,i4p) = random_normal(0, (R_Q(i4p)) , (/np/))  ; [mg/kg]       ; jkmod7
     end do

;jkmod11 prevent the radom_error obserr_P being too larger than the standard devation of the normal dist.
     do i4p=0, np-1
        if ((.not. ismissing(obserr_T(mi(i4p),i4p))) .and.  (obserr_T(mi(i4p),i4p) .gt. R_T(i4p))) then
           print(i4p+" level before obserr_T="+obserr_T(mi(i4p),i4p))
           obserr_T(mi(i4p),i4p) = abs(obserr_T(mi(i4p),i4p))/obserr_T(mi(i4p),i4p)*R_T(i4p)
           print(i4p+" level after  obserr_T="+obserr_T(mi(i4p),i4p))
        end if
     end do

     do i4p=0, np-1
        if ((.not. ismissing(obserr_P(mi(i4p),i4p))) .and.  (obserr_P(mi(i4p),i4p) .gt. R_P(i4p))) then
        ;if (obserr_P(mi(i4p),i4p) .gt. R_P(i4p)) then
           print(i4p+" level before obserr_P="+obserr_P(mi(i4p),i4p))
           obserr_P(mi(i4p),i4p) = abs(obserr_P(mi(i4p),i4p))/obserr_P(mi(i4p),i4p)*R_P(i4p)
           print(i4p+" level after  obserr_P="+obserr_P(mi(i4p),i4p))
        end if
     end do

     do i4p=0, np-1
        if (.not. ismissing(obserr_Q(mi(i4p),i4p)) .and.(.not. ismissing(R_Q(i4p))))  then
           if ( (obserr_Q(mi(i4p),i4p) .gt. R_Q(i4p))) then
               print(i4p+" level before obserr_Q="+obserr_Q(mi(i4p),i4p))
               obserr_Q(mi(i4p),i4p) = abs(obserr_Q(mi(i4p),i4p))/obserr_Q(mi(i4p),i4p)*R_Q(i4p)
               print(i4p+" level after  obserr_Q="+obserr_Q(mi(i4p),i4p))
          end if
        end if
     end do


;-------------------------------------------------------------------------------------------------------

     R_Tz = new((/fnz/),float)
     R_RHz = new((/fnz/),float)
     R_Pz = new((/fnz/),float)
     R_Tz = 0.0
     R_RHz = 0.0
     R_Pz = 0.0

    ; calculate observation error corresponding to pressure level
     do i=0, fnz-1
        if ((ismissing(tv_d(i)) .eq. False) .and. (ismissing(q_d(i)) .eq. False) .and. \
            (q_d(i) .gt. 0.0)) then
          
             do i4p=0, np-2
                if ((fpob(i) .le. Rlev(i4p)) .and. (fpob(i) .gt. Rlev(i4p+1))) then
                    print("level "+i+" : "+fpob(i)+" plevel Q: "+q_d(i))

                    tv_d(i) = tv_d(i)+obserr_T(mi(i4p),i4p)
                    if (q_d(i) .lt. obserr_Q(mi(i4p),i4p)) then
                        q_d(i) = q_d(i)*1.5
                    else if (q_d(i) + obserr_Q(mi(i4p),i4p) .lt. 0.0) then
                        q_d(i) = q_d(i)*0.5
                    else 
                        q_d(i) = q_d(i)+obserr_Q(mi(i4p),i4p)
                    end if
                    end if

                    if (q_d(i) .lt. 0.0) then
                        q_d(i) = 0.0
                    end if

                    R_Tz(i) = R_factor*R_T(i4p)
                    R_RHz(i) = R_factor*R_RH(i4p)/10.
                    R_Pz(i) = R_factor*R_P(i4p)

                    print("level "+i+" : "+fpob(i)+" plevel Q: "+q_d(i))
                end if
            end do
        else
            tv_d(i) = 100000000000.0 
            q_d(i)  = 100000000000.0

            R_Tz(i) = 100000000000.0 
            R_RHz(i)= 100000000000.0
            R_Pz(i) = 100000000000.0   ; Psfc observation will not be assimilated
        end if
     end do

     ; Set SONDE_SFC observation
     p_sfc = p_sfc+obserr_P(mi(0),0)    ; jkmod7
     tv_sfc = tv_sfc+obserr_T(mi(0),0)  ; jkmod7
     q_sfc = q_sfc+obserr_Q(mi(0),0)    ; jkmod7

     obs_model(0:fnz-1,1,0,i4o) = q_d(0:fnz-1)
     obs_model(0:fnz-1,2,0,i4o) = tv_d(0:fnz-1)
     obs_model(0:fnz-1,3,0,i4o) = z_d(0:fnz-1)
     ;obs_model(0,0,0,i4o) = p_sfc
     ;obs_model(0,1,0,i4o) = q_sfc
     ;obs_model(0,2,0,i4o) = tv_sfc
     ;obs_model(0,3,0,i4o) = z_sfc

     obs_model(0:fnz-1,0,3,i4o) = (/ R_Pz(0:fnz-1) /)     ; [mb]
     obs_model(0:fnz-1,1,3,i4o) = (/ R_RHz(0:fnz-1) /)    ; [%/10]
     obs_model(0:fnz-1,2,3,i4o) = (/ R_Tz(0:fnz-1) /)

     print("before_Z : "+hdr_model(5,i4o))
     hdr_model(5,i4o) = z_sfc
     print("after_Z  : "+hdr_model(5,i4o))

; remove two extra sonde observation
     print("jk_lonlat="+hdr(1,i4o)+" "+hdr(2,i4o))
     do i4s=0, nasonde-1
        if (abs(hdr(1,i4o)-lon_asonde(i4s)) .le. 0.01 .and. abs(hdr(2,i4o) - lat_asonde(i4s)) .le. 0.01) then
           print("jkjkjk")
           obs_model(:,0,1,i4o) = 14
           obs_model(:,1,1,i4o) = 14
           obs_model(:,2,1,i4o) = 14
           obs_model(:,3,1,i4o) = 14
        end if
     end do
 
     delete(p_sfc)
     delete(q_sfc)
     delete(tv_sfc)
     delete(z_sfc)

     delete(mi)
     delete(i)

     delete(R_Tz)
     delete(R_RHz)
     delete(R_Pz)

     delete(obserr_P)
     delete(obserr_T)
     delete(obserr_Q)
     ;delete(Rname)
     ;delete(Rdata)
     delete(R_T)
     delete(R_Q)
     delete(R_RH)

     delete(t)
     delete(tv)
     delete(qv)
     delete(q)
     delete(qin)
     delete(tvin)
     ;delete(zin)
     delete(q_d)
     delete(tv_d)
     delete(z_d)

  ;-------------------------------------------------------------
  ; WIND observation
  ;-------------------------------------------------------------
  else if (hdr(4,i4o) .ge. 200.) then      ; WIND

     uvmet = wrf_user_getvar(f,"uvmet",0)       ; Wind rotated to earth coord.
     u = uvmet(0,:,:,:)
     v = uvmet(1,:,:,:)

     uin = u(:,obsi,obsj) ; level, lat, lon
     vin = v(:,obsi,obsj)

     u_d = wrf_interp_1d(uin,pin,fpob)
     v_d = wrf_interp_1d(vin,pin,fpob)

     p_sfc = flt2dble(pin(0))  ; jkmod3
     z_sfc = flt2dble(zin(0))

     ;u_d(ind(ismissing(u_d)))=100000000000.0
     ;v_d(ind(ismissing(v_d)))=100000000000.0

     do i4z=0, fnz-1
        print(sprintf("%8.2f",obs(i4z,0,0,i4o))+" hPa, MODEL U: "+sprintf("%15.2f", u_d(i4z))+" / OBS : "+sprintf("%15.2f",obs(i4z,4,0,i4o)))
     end do

    ; READ OBSERVATION ERROR
     ;Rname = hdr(4,i4o)+"_observation_R.txt"  ;jkmod
     ;Rdata  = asciiread(Rname,-1,"string")
     ;R_W  = stringtofloat(str_get_cols(Rdata,37,48))
     ;print("RH_Rerr="+R_W)

     i = 0
     do while(typ(i) .ne. hdr(4,i4o))
       i=i+1
     end do
     print(hdr(4,i4o)+" is "+i)
     R_W  = stringtofloat(str_get_cols(Rdata(i+1:i+33),37,48))
     ;print("RH_Rerr="+R_W)
     

     ; normal distribution random noise (av = average, sd = standard deviation)
     ; random_setallseed(36484749, 9494848)               ; Set seeds (suggested, NOT required)
     rand1 = toint(systemfunc(" date +%s"))       ; jkmod7
     rand2 = toint(systemfunc(" date +%s"))       ; jkmod7
     random_setallseed(rand1,rand2)               ; Set seeds (suggested, NOT required)
     mi=generate_sample_indices(np , 0 )         ;jkmod7


     obserr_U = new((/np,np/),float)
     obserr_V = new((/np,np/),float)
     obserr_U = 0.0
     obserr_V = 0.0
     do i4p=0,np-1
        obserr_U(:,i4p) = random_normal(0, (R_W(i4p)) , (/np/))  ; jkmod7
        obserr_V(:,i4p) = random_normal(0, (R_W(i4p)) , (/np/))  ; use R_W sd, not varinace
     end do

;jkmod11 prevent the radom_error obserr_P being too larger than the standard devation of the normal dist.
     do i4p=0, np-1
        if ((.not. ismissing(R_W(i4p))) .and.  (obserr_U(mi(i4p),i4p) .gt. R_W(i4p))) then
           print(i4p+" level before obserr_U="+obserr_U(mi(i4p),i4p))
           obserr_U(mi(i4p),i4p) = abs(obserr_U(mi(i4p),i4p))/obserr_U(mi(i4p),i4p)*R_W(i4p)
           print(i4p+" level after  obserr_U="+obserr_U(mi(i4p),i4p))
        end if
     end do

     do i4p=0, np-1
        if ((.not. ismissing(R_W(i4p))) .and.  (obserr_V(mi(i4p),i4p) .gt. R_W(i4p))) then
           print(i4p+" level before obserr_V="+obserr_V(mi(i4p),i4p))
           obserr_V(mi(i4p),i4p) = abs(obserr_V(mi(i4p),i4p))/obserr_V(mi(i4p),i4p)*R_W(i4p)
           print(i4p+" level after  obserr_V="+obserr_V(mi(i4p),i4p))
        end if
     end do

    ; calculate observation error corresponding to pressure level
     R_Wz = new((/fnz/),float)
     R_Wz = 0.0
     do i=0, fnz-1
        if ((ismissing(u_d(i)) .eq. False) .and. (ismissing(v_d(i)) .eq. False)) then 
            do i4p=0, np-2
                if  ((fpob(i) .le. Rlev(i4p)) .and. (fpob(i) .gt. Rlev(i4p+1))) then
                    print("level "+i+" : "+fpob(i)+" plevel U: "+u_d(i))

                    u_d(i) = u_d(i)+obserr_U(mi(i4p),i4p)
                    v_d(i) = v_d(i)+obserr_V(mi(i4p),i4p)
                    R_Wz(i) = R_factor*R_W(i4p)

                    print("level "+i+" : "+fpob(i)+" plevel U: "+u_d(i))
                end if
            end do
        else 
            u_d(i) = 100000000000.0 
            v_d(i) = 100000000000.0
            R_Wz(i)= 100000000000.0
        end if
     end do

     obs_model(0:fnz-1,4,0,i4o) = u_d(0:fnz-1)
     obs_model(0:fnz-1,5,0,i4o) = v_d(0:fnz-1)
     obs_model(0,0,0,i4o) = p_sfc   ; jkmod3

     obs_model(0:fnz-1,4,3,i4o) = (/ R_Wz(0:fnz-1) /)
     obs_model(0:fnz-1,5,3,i4o) = (/ R_Wz(0:fnz-1) /)

     print("before : "+hdr_model(5,i4o))
     hdr_model(5,i4o) = z_sfc
     print("after  : "+hdr_model(5,i4o))

; remove two extra sonde observation
     print("jk_lonlat="+hdr(1,i4o)+" "+hdr(2,i4o))
     do i4s=0, nasonde-1
        if (abs(hdr(1,i4o)-lon_asonde(i4s)) .le. 0.01 .and. abs(hdr(2,i4o) - lat_asonde(i4s)) .le. 0.01) then
           print("jkjk")
           obs_model(:,4,1,i4o) = 14
           obs_model(:,5,1,i4o) = 14
        end if
     end do

     delete(p_sfc)
     delete(z_sfc)

     delete(mi)
     delete(i)

     delete(R_Wz)
     delete(R_W)
     ;delete(Rname)
     ;delete(Rdata)

     delete(obserr_U)
     delete(obserr_V)
     delete(u_d)
     delete(v_d)
     delete(uin)
     delete(vin)

     delete(uvmet)
     delete(u)
     delete(v)

   end if
   end if

     delete(obsij)
     delete(obsi)
     delete(obsj)

     delete(pin)
     delete(zin)

     delete(flon)
     delete(flat)
     delete(fnz)
     delete(fpob)
     delete(fzob)

     delete(in_file)
     delete(mmm)
     delete(f)
     delete(z)
     delete(p)

  end do ; for i4o

  obs_model_1D      = ndtooned (obs_model)
  printVarSummary (obs_model)
  printVarSummary (obs_model_1D)

  if (any(ismissing(obs_model_1D))) then
     obs_model_1D(ind(ismissing(obs_model_1D))) = 100000000000.0
  end if

  obs_model = onedtond(obs_model_1D, dimsizes(obs_model))
  delete(obs_model_1D)
;  print(obs_model(:,2,0,num_msg_sb-2))
;-------------------------------------------------------
; Write out interpolated data
;-------------------------------------------------------
  file_sim_obs = path_hdr+"/"+date+"/"+date+"_"+obs_typ+"_V${OBS_DHR}_gdas1_sim_obs.dat"

  system("/bin/rm -f "+file_sim_obs)

  setfileoption("bin","WriteByteOrder","BigEndian")
  fbindirwrite(file_sim_obs,obs_model)

;--header file
  file_sim_hdr = path_hdr+"/"+date+"/"+date+"_"+obs_typ+"_V${OBS_DHR}_gdas1_sim_hdr.dat"
  system("/bin/rm -f "+file_sim_hdr)

  setfileoption("bin","WriteByteOrder","BigEndian")
  fbindirwrite(file_sim_hdr,hdr_model)
end

EOF

echo "Calculating simulated observation of " $OBS_TYPE "at " $OBS_DATE "V"$OBS_DHR 
/glade/u/apps/ch/opt/ncl/6.6.2/intel/18.0.5/bin/ncl $WORK_DIR/$OBS_DATE/mk_${OBS_TYPE}_${OBS_DATE}_nrobs_${OBS_DHR}.ncl 
echo "Finish" $OBS_TYPE "at " $OBS_DATE "V"$OBS_DHR

done
done
done
