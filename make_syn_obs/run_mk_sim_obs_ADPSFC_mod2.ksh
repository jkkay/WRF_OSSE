#!/bin/ksh

export WORK_DIR=/glade/work/junkyung/GSI/comGSIv3.7_EnKFv1.3/util/bufr_tools/knu
export OBS_DATES="2022011100" #  2015071418 2015071500 2015071506" 
export OBS_TYPES="ADPSFC" #"PROFLR" #"SFCSHP" # "ADPSFC SFCSHP"
export OBS_DHRS="00 m3 p3"

for OBS_DATE in $OBS_DATES; do
mkdir -p $WORK_DIR/$OBS_DATE
cd $WORK_DIR/$OBS_DATE

for OBS_TYPE in $OBS_TYPES; do
for OBS_DHR in $OBS_DHRS; do
cat > mk_${OBS_TYPE}_${OBS_DATE}_nrobs_${OBS_DHR}.ncl << EOF
load "\$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "\$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "\$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
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
  obs_typ="${OBS_TYPE}"
  path_hdr = "${WORK_DIR}/"
  file_hdr = path_hdr+date+"_"+obs_typ+"_V${OBS_DHR}_gdas1_hdr.dat"
  file_obs = path_hdr+date+"_"+obs_typ+"_V${OBS_DHR}_gdas1_obs.dat"

  ;wrf_dir="/glade/scratch/junkyung/KNU/ICBC/rundir/test/era5/usphy/wrfrun/fc/2022011012.e042/"
  wrf_dir="/glade/scratch/junkyung/KNU/ICBC/rundir/test/era5/usphy/run/2022011012/wrf.e042_1min/"
  print("file_hdr= "+file_hdr)

  ;FILESIZE=systemfunc("wc -c "+file_hdr+" | awk '{print $1}'")
  FILESIZE=systemfunc("wc -c "+file_hdr)
  delim = " "
  num_msg_sb = stringtointeger(str_get_field(FILESIZE, 1, delim))/64
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
printVarSummary(obs_model)

  ;print(hdr(:,0))
  ;print(hdr(:,num_msg_sb-1))

  ;print(obs(:,2,0,num_msg_sb-2))
;----------------------------------------------------------
; Check Time of the observation within window
;----------------------------------------------------------
  time_00 = new(1,integer)
  time_m3 = new(1,integer)
  time_p3 = new(1,integer)
  time_tot = new(1,integer)
  time_00 = 0
  time_m3 = 0
  time_p3 = 0
  time_tot = 0

  do i4o=0, num_msg_sb-1
     time_tot = time_tot + 1
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
  print("Total time number = "+time_tot)
  print("T+00hr +/- 3hr obs number = "+time_00)
  print("< T-03hr obs number = "+time_m3)
  print("> T+03hr obs number = "+time_p3)

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

;----------------------------------------------------------
; Read nature run according to the observation time
;----------------------------------------------------------
  do i4o=0 , num_msg_sb-1
     ;dhr_min=hdr(3,i4o)*60.
     dhr_min=round(hdr(3,i4o)*60.,0)
     obs_time=systemfunc("/glade/u/home/junkyung/WRF/V4.1.2/WRFDA/var/build/da_advance_time.exe " + date + " "+dhr_min+"m")
     ;print("input obs time : "+hdr(3,i4o)+" "+dhr_min+" "+obs_time)

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

     ;print("final obs time is "+yyyy+" "+mm+" "+dd+" "+hh+" "+mmm)

     ;-------------------------------------------------------------
     ; READ WRF output
     ;-------------------------------------------------------------
     in_file=wrf_dir+"/wrfout_d01_"+yyyy+"-"+mm+"-"+dd+"_"+hh+":"+mmm+":00"
     print(i4o+"th input wrf file = "+in_file)

     ; For the very first cycling valid at 1200 UTC 14 July (071412)
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

     ;   print(ni+" "+nj)
     end if

     flon = hdr(1,i4o)
     flat = hdr(2,i4o)
     fnz = doubletoint(hdr(7,i4o))
     fpob = obs(0:fnz-1,0,0,i4o)
     fzob = obs(0:fnz-1,3,0,i4o)

     ;print(i4o+" th nz is "+fnz)

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

     rh= wrf_user_getvar(f,"rh",0)

     qin  = q(:,obsi,obsj) ; level, lat, lon
     tvin = tv(:,obsi,obsj)
     ;zin  = z(:,obsi,obsj)
     rhin = rh(:,obsi,obsj)

     ;q_d  = wrf_interp_1d(qin,pin,fpob)
     ;tv_d = wrf_interp_1d(tvin,pin,fpob)
     ;z_d  = wrf_interp_1d(zin,pin,fpob)
     p_d = flt2dble(pin(0))  ; jkmod3
     q_d = flt2dble(qin(0))
     tv_d = flt2dble(tvin(0))
     z_d = flt2dble(zin(0))


     ;print("jk_ij="+obsi+" "+obsj)
     ;print("jk_lon_lat="+flon+" "+flat)
     ;print("jk_zin="+zin(0))
     ;print("jk_pin="+pin(0))
     ;print("jkobs_p="+fpob(0))
     ;print("jkobs_z="+fzob(0))

     ;q_d(ind(ismissing(q_d)))=100000000000.0
     ;z_d(ind(ismissing(z_d)))=100000000000.0
     ;tv_d(ind(ismissing(tv_d)))=100000000000.0

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
     ;print(hdr(4,i4o)+" is "+i)
     P0  = stringtofloat(str_get_cols(Rdata(i+1:i+33), 0,12))  ; [mb] in prepbufr & Rtable
     R_P  = stringtofloat(str_get_cols(Rdata(i+1:i+33),49,60)) ; [mb]
     R_T  = stringtofloat(str_get_cols(Rdata(i+1:i+33),13,24)) ; [K]
     R_RH = stringtofloat(str_get_cols(Rdata(i+1:i+33),25,36)) ; [%/10] in prepbufr & Rtable
     R_RH=R_RH*10. ; [%/10]-->[%}

    ; CONVERT rh observation error (R_RH) to specific humidity observation error (R_Q) (da_get_q_error.inc)
     ;tv_r = wrf_interp_1d(tvin,pin,R_P)
     ;rh_r = wrf_interp_1d(rhin,pin,R_P)
     ;q_r  = wrf_interp_1d(qin,pin,R_P)
     tv_r = tvin(0)
     rh_r = rhin(0)
     q_r = qin(0)

     R_Q = new((/1/),float,0.10000E+10)
     print("P_level        R_RH        R_T        rh         tv         q        R_Q ")


     if ((ismissing(tv_r) .eq. False) .and. (ismissing(rh_r) .eq. False)) then
           ;R_Q(0) = da_get_q_error(R_RH(0), R_T(0), rh_r, tv_r, R_P(0))  ; (f_qv_from_rh)
           
           ; da_get_q_error : R_RH [%], R_T[K], rh_r[%], tv_r[C], P0[mb] in f_qv_from_rh.f 
           R_Q(0) = da_get_q_error(R_RH(0), R_T(0), rh_r, tv_r, P0(0))  ; (f_qv_from_rh)
           R_Q(0) = R_Q(0)*1000000.  ; [kg/kg] --> [mg/kg]

           print(sprintf("%10.6f",P0(0))+" "+sprintf("%10.6f",R_RH(0))+" "+sprintf("%10.6f",R_T(0))+" "+ \
                 sprintf("%10.6f",rh_r)+" "+sprintf("%10.6f",tv_r)+" "+sprintf("%10.f",q_r)+" "+sprintf("%10.6f",R_Q(0)))
     end if

     ; normal distribution random noise (av = average, sd = standard deviation)
     ; random_setallseed(36484749, 9494848)               ; Set seeds (suggested, NOT required)
     rand1 = toint(systemfunc(" date +%s"))       ; jkmod7
     rand2 = toint(systemfunc(" date +%s"))       ; jkmod7
     random_setallseed(rand1,rand2)               ; Set seeds (suggested, NOT required)
     fns = 100 ; dimension number of seed ;       jkmod10
     mi=generate_sample_indices( fns, 0 )         ;jkmod7

     obserr_P = new((/fns,1/),float)                          ;jkmod3
     obserr_T = new((/fns,1/),float)
     obserr_Q = new((/fns,1/),float)
     obserr_T = 0.0
     obserr_Q = 0.0
     do i4p=0,0
        obserr_P(:,i4p) = random_normal(0, 1.0*R_P(i4p) , (/fns/))  ;jkmod3
        obserr_T(:,i4p) = random_normal(0, 1.0*R_T(i4p) , (/fns/))
        obserr_Q(:,i4p) = random_normal(0, 1.0*R_Q(i4p) , (/fns/))
     end do

;jkmod11 prevent the radom_error obserr_P being too larger than the standard devation of the normal dist.
     if (obserr_T(mi(0),0) .gt. R_T(0)) then
         print("before obserr_T="+obserr_T(mi(0),0))
         obserr_T(mi(0),0) = abs(obserr_T(mi(0),0))/obserr_T(mi(0),0)*R_T(0)
         print("after  obserr_T="+obserr_T(mi(0),0))
     end if

     if (obserr_P(mi(0),0) .gt. R_P(0)) then
         print("before obserr_P="+obserr_P(mi(0),0))
         obserr_P(mi(0),0) = abs(obserr_P(mi(0),0))/obserr_P(mi(0),0)*R_P(0)
         print("after  obserr_P="+obserr_P(mi(0),0))
     end if

     if (obserr_Q(mi(0),0) .gt. R_Q(0)) then
         print("before obserr_Q="+obserr_Q(mi(0),0))
         obserr_Q(mi(0),0) = abs(obserr_Q(mi(0),0))/obserr_Q(mi(0),0)*R_Q(0)
         print("after obserr_Q="+obserr_Q(mi(0),0))
     end if
;--------------------------------------------------------------------------------------------------------

     tv_d = tv_d+obserr_T(mi(0),0)
     q_d = q_d+obserr_Q(mi(0),0)
     p_d = p_d+obserr_P(mi(0),0)
    ; print("level 0 : "+fpob(0)+" plevel Q: "+q_d(0))

     obs_model(0:fnz-1,0,0,i4o) = p_d ;p_d(0:fnz-1)   ; [mb]  ;jkmod3
     obs_model(0:fnz-1,1,0,i4o) = q_d ;q_d(0:fnz-1)   ; [mg/kg]
     obs_model(0:fnz-1,2,0,i4o) = tv_d ;tv_d(0:fnz-1) ; [c]
     obs_model(0:fnz-1,3,0,i4o) = z_d ;z_d(0:fnz-1)   ; [m]

     obs_model(0:fnz-1,0,3,i4o) = R_P(0)  ; [mb]
     obs_model(0:fnz-1,1,3,i4o) = R_RH(0)/10. ; [%/10]
     obs_model(0:fnz-1,2,3,i4o) = R_T(0)

     print("before : "+hdr_model(5,i4o))
     hdr_model(5,i4o) = z_d
     print("after  : "+hdr_model(5,i4o))
     delete(i)
     delete(mi)

     delete(tv_r)
     delete(q_r)
     delete(rh_r)

     delete(obserr_P)
     delete(obserr_T)
     delete(obserr_Q)
     ;delete(Rname)
     ;delete(Rdata)
     delete(R_P)
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
     delete(p_d)
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

     p_d = flt2dble(pin(0))  ; jkmod3
     u_d = flt2dble(uin(0))
     v_d = flt2dble(vin(0))
     z_d = flt2dble(zin(0))

     ;u_d = wrf_interp_1d(uin,pin,fpob)
     ;v_d = wrf_interp_1d(vin,pin,fpob)
;printVarSummary(u_d)

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
     ;print(hdr(4,i4o)+" is "+i)
     R_W  = stringtofloat(str_get_cols(Rdata(i+1:i+33),37,48))
     R_P  = stringtofloat(str_get_cols(Rdata(i+1:i+33),49,60)) ; [mb]       ;jkmod3 ; this value is 0.10000E+10 
     ;print("W_Rerr="+R_W(0))

     ; normal distribution random noise (av = average, sd = standard deviation)
     ; random_setallseed(36484749, 9494848)               ; Set seeds (suggested, NOT required)
     rand1 = toint(systemfunc(" date +%s"))       ; jkmod7
     rand2 = toint(systemfunc(" date +%s"))       ; jkmod7
     random_setallseed(rand1,rand2)               ; Set seeds (suggested, NOT required)
     fns = 100 ; dimension number of seed ;       jkmod10
     mi=generate_sample_indices( fns, 0 )         ;jkmod7

     obserr_P = new((/fns,1/),float)     ;jkmod3
     obserr_U = new((/fns,1/),float)
     obserr_V = new((/fns,1/),float)
     obserr_P = 0.0
     obserr_U = 0.0
     obserr_V = 0.0
     do i4p=0,0
        obserr_P(:,i4p) = random_normal(0, R_P(i4p) , (/fns/))   ;jkmod3
        obserr_U(:,i4p) = random_normal(0, R_W(i4p) , (/fns/))
        obserr_V(:,i4p) = random_normal(0, R_W(i4p) , (/fns/))
     end do

;jkmod11 prevent the radom_error obserr_P being too larger than the standard devation of the normal dist.
     if (obserr_P(mi(0),0) .gt. R_P(0)) then
         print("before obserr_P="+obserr_P(mi(0),0))
         obserr_P(mi(0),0) = abs(obserr_P(mi(0),0))/obserr_P(mi(0),0)*R_P(0)
         print("after  obserr_P="+obserr_P(mi(0),0))
     end if

     if (obserr_U(mi(0),0) .gt. R_W(0)) then
         print("before obserr_U="+obserr_U(mi(0),0))
         obserr_U(mi(0),0) = abs(obserr_U(mi(0),0))/obserr_U(mi(0),0)*R_W(0)
         print("after obserr_U="+obserr_U(mi(0),0))
     end if

     if (obserr_V(mi(0),0) .gt. R_W(0)) then
         print("before obserr_V="+obserr_V(mi(0),0))
         obserr_V(mi(0),0) = abs(obserr_V(mi(0),0))/obserr_V(mi(0),0)*R_W(0)
         print("after obserr_V="+obserr_V(mi(0),0))
     end if
;--------------------------------------------------------------------------------------------------------

    ; print("obserr_U :"+obserr_U(0,0))
    ; calculate observation error corresponding to pressure level
     do i=0, fnz-1
        do i4p=0, 0 ;np-2
            ;if ((ismissing(u_d(i)) .eq. False) .and. (ismissing(v_d(i)) .eq. False) .and. \
            ;    (fpob(i) .le. Rlev(i4p)) .and. (fpob(i) .gt. Rlev(i4p+1))) then
            if ((ismissing(u_d(i)) .eq. False) .and. (ismissing(v_d(i)) .eq. False)) then
                print("before level "+i+" : "+fpob(i)+" plevel U: "+u_d)

                p_d = p_d+obserr_P(mi(i),i4p)     ;jkmod3
                u_d = u_d+obserr_U(mi(i),i4p)
                v_d = v_d+obserr_V(mi(i),i4p)

                print("after  level "+i+" : "+fpob(i)+" plevel U: "+u_d)
            end if
        end do
     end do


     obs_model(0:fnz-1,0,0,i4o) = p_d(0:fnz-1)
     obs_model(0:fnz-1,4,0,i4o) = u_d(0:fnz-1)
     obs_model(0:fnz-1,5,0,i4o) = v_d(0:fnz-1)

     obs_model(0:fnz-1,0,3,i4o) = R_P(0)
     obs_model(0:fnz-1,4,3,i4o) = R_W(0)
     obs_model(0:fnz-1,5,3,i4o) = R_W(0) 

     print("before : "+hdr_model(5,i4o))
     hdr_model(5,i4o) = z_d
     print("after  : "+hdr_model(5,i4o))

;printVarSummary(obs_model)

     delete(i)
     delete(mi)

     delete(R_P)
     delete(R_W)
     ;delete(Rname)
     ;delete(Rdata)

     delete(obserr_P)
     delete(obserr_U)
     delete(obserr_V)
     delete(p_d)
     delete(z_d)
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


  print("jkcheck1")
  obs_model_1D      = ndtooned (obs_model(0,:,:,:))
  printVarSummary (obs_model)
  printVarSummary (obs_model_1D)
  print(obs_model_1D)
  print(ismissing(obs_model_1D))

  if(any(ismissing(obs_model_1D))) then
    print("x contains one or more missing values, cannot continue.")
    ;return
  end if

  print(ind(ismissing(obs_model_1D)))

;-------------------------------------------------------
; Write out interpolated data
;-------------------------------------------------------
;--observation file
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
#/glade/u/apps/ch/opt/ncl/6.6.2/intel/18.0.5/bin/ncl $WORK_DIR/$OBS_DATE/mk_${OBS_TYPE}_${OBS_DATE}_nrobs_${OBS_DHR}.ncl #> log_${OBS_TYPE}_${OBS_DATE}_${OBS_DHR}
echo "Finish" $OBS_TYPE "at " $OBS_DATE "V"$OBS_DHR

done
done
done
