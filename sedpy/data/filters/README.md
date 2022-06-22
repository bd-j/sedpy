# Filters in `sedpy`

## Space-Based

|Filter|Telescope|Instrument|`sedpy` Name|
|:----:|:-------:|:--------:|:----------:|
|$G$|Gaia||[`gaia_g`](gaia_g.par)|
|$G_\textrm{BP}$|Gaia|BP|[`gaia_bp`](gaia_bp.par)|
|$G_\textrm{RP}$|Gaia|RP|[`gaia_rp`](gaia_rp.par)|
|$\textrm{FUV}$|GALEX|FUV|[`galex_FUV`](galex_FUV.par)|
|$\textrm{NUV}$|GALEX|NUV|[`galex_NUV`](galex_NUV.par)|
|$70\ \textrm{μm}$|Herschel|PACS|[`herschel_pacs_70`](herschel_pacs_70.par)|
|$100\ \textrm{μm}$|Herschel|PACS|[`herschel_pacs_100`](herschel_pacs_100.par)|
|$160\ \textrm{μm}$|Herschel|PACS|[`herschel_pacs_160`](herschel_pacs_160.par)|
|$250\ \textrm{μm}$|Herschel|SPIRE|[`herschel_spire_250`](herschel_spire_250.par)|
|$350\ \textrm{μm}$|Herschel|SPIRE|[`herschel_spire_350`](herschel_spire_350.par)|
|$500\ \textrm{μm}$|Herschel|SPIRE|[`herschel_spire_500`](herschel_spire_500.par)|
|$\textrm{F435W}$|Hubble|ACS (WFC)|[`acs_wfc_f435w`](acs_wfc_f435w.par)|
|$\textrm{F475W}$|Hubble|ACS (WFC)|[`acs_wfc_f475w`](acs_wfc_f475w.par)|
|$\textrm{F555W}$|Hubble|ACS (WFC)|[`acs_wfc_f555w`](acs_wfc_f555w.par)|
|$\textrm{F606W}$|Hubble|ACS (WFC)|[`acs_wfc_f606w`](acs_wfc_f606w.par)|
|$\textrm{F625W}$|Hubble|ACS (WFC)|[`acs_wfc_f625w`](acs_wfc_f625w.par)|
|$\textrm{F775W}$|Hubble|ACS (WFC)|[`acs_wfc_f775w`](acs_wfc_f775w.par)|
|$\textrm{F814W}$|Hubble|ACS (WFC)|[`acs_wfc_f814w`](acs_wfc_f814w.par)|
|$\textrm{F435LP}$|Hubble|ACS (WFC)|[`acs_wfc_f850lp`](acs_wfc_f850lp.par)|
|$\textrm{F275W}$|Hubble|WFC3 (UVIS)|[`wfc3_uvis_f275w`](wfc3_uvis_f275w.par)|
|$\textrm{F336W}$|Hubble|WFC3 (UVIS)|[`wfc3_uvis_f336w`](wfc3_uvis_f336w.par)|
|$\textrm{F475W}$|Hubble|WFC3 (UVIS)|[`wfc3_uvis_f475w`](wfc3_uvis_f475w.par)|
|$\textrm{F555W}$|Hubble|WFC3 (UVIS)|[`wfc3_uvis_f555w`](wfc3_uvis_f555w.par)|
|$\textrm{F606W}$|Hubble|WFC3 (UVIS)|[`wfc3_uvis_f606w`](wfc3_uvis_f606w.par)|
|$\textrm{F814W}$|Hubble|WFC3 (UVIS)|[`wfc3_uvis_f814w`](wfc3_uvis_f814w.par)|
|$\textrm{F105W}$|Hubble|WFC3 (IR)|[`wfc3_ir_f105w`](wfc3_ir_f105w.par)|
|$\textrm{F110W}$|Hubble|WFC3 (IR)|[`wfc3_ir_f110w`](wfc3_ir_f110w.par)|
|$\textrm{F125W}$|Hubble|WFC3 (IR)|[`wfc3_ir_f125w`](wfc3_ir_f125w.par)|
|$\textrm{F140W}$|Hubble|WFC3 (IR)|[`wfc3_ir_f140w`](wfc3_ir_f140w.par)|
|$\textrm{F160W}$|Hubble|WFC3 (IR)|[`wfc3_ir_f160w`](wfc3_ir_f160w.par)|
|$\textrm{Band}\ 1$|Spitzer|IRAC|[`spitzer_irac_ch1`](spitzer_irac_ch1.par)|
|$\textrm{Band}\ 2$|Spitzer|IRAC|[`spitzer_irac_ch2`](spitzer_irac_ch2.par)|
|$\textrm{Band}\ 3$|Spitzer|IRAC|[`spitzer_irac_ch3`](spitzer_irac_ch3.par)|
|$\textrm{Band}\ 4$|Spitzer|IRAC|[`spitzer_irac_ch4`](spitzer_irac_ch4.par)|
|$24\ \textrm{μm}$|Spitzer|MIPS|[`spitzer_mips_24`](spitzer_mips_24.par)|
|$70\ \textrm{μm}$|Spitzer|MIPS|[`spitzer_mips_70`](spitzer_mips_70.par)|
|$160\ \textrm{μm}$|Spitzer|MIPS|[`spitzer_mips_160`](spitzer_mips_160.par)|
|$\textrm{F070W}$|Webb|NIRCam|[`jwst_f070w`](jwst_f070w.par)|
|$\textrm{F090W}$|Webb|NIRCam|[`jwst_f090w`](jwst_f090w.par)|
|$\textrm{F115W}$|Webb|NIRCam|[`jwst_f115w`](jwst_f115w.par)|
|$\textrm{F150W}$|Webb|NIRCam|[`jwst_f150w`](jwst_f150w.par)|
|$\textrm{F200W}$|Webb|NIRCam|[`jwst_f200w`](jwst_f200w.par)|
|$\textrm{F277W}$|Webb|NIRCam|[`jwst_f277w`](jwst_f277w.par)|
|$\textrm{F356W}$|Webb|NIRCam|[`jwst_f356w`](jwst_f356w.par)|
|$\textrm{F444W}$|Webb|NIRCam|[`jwst_f444w`](jwst_f444w.par)|
|$\textrm{F140M}$|Webb|NIRCam|[`jwst_f140m`](jwst_f140m.par)|
|$\textrm{F162M}$|Webb|NIRCam|[`jwst_f162m`](jwst_f162m.par)|
|$\textrm{F182M}$|Webb|NIRCam|[`jwst_f182m`](jwst_f182m.par)|
|$\textrm{F210M}$|Webb|NIRCam|[`jwst_f210m`](jwst_f210m.par)|
|$\textrm{F250M}$|Webb|NIRCam|[`jwst_f250m`](jwst_f250m.par)|
|$\textrm{F300M}$|Webb|NIRCam|[`jwst_f300m`](jwst_f300m.par)|
|$\textrm{F335M}$|Webb|NIRCam|[`jwst_f335m`](jwst_f335m.par)|
|$\textrm{F360M}$|Webb|NIRCam|[`jwst_f360m`](jwst_f360m.par)|
|$\textrm{F410M}$|Webb|NIRCam|[`jwst_f410m`](jwst_f410m.par)|
|$\textrm{F430M}$|Webb|NIRCam|[`jwst_f430m`](jwst_f430m.par)|
|$\textrm{F460M}$|Webb|NIRCam|[`jwst_f460m`](jwst_f460m.par)|
|$\textrm{F480M}$|Webb|NIRCam|[`jwst_f480m`](jwst_f480m.par)|
|$\textrm{F560W}$|Webb|MIRI|[`jwst_f560w`](jwst_f560w.par)|
|$\textrm{F770W}$|Webb|MIRI|[`jwst_f770w`](jwst_f770w.par)|
|$\textrm{F1000W}$|Webb|MIRI|[`jwst_f1000w`](jwst_f1000w.par)|
|$\textrm{F1130W}$|Webb|MIRI|[`jwst_f1130w`](jwst_f1130w.par)|
|$\textrm{F1280W}$|Webb|MIRI|[`jwst_f1280w`](jwst_f1280w.par)|
|$\textrm{F1500W}$|Webb|MIRI|[`jwst_f1500w`](jwst_f1500w.par)|
|$\textrm{F1800W}$|Webb|MIRI|[`jwst_f1800w`](jwst_f1800w.par)|
|$\textrm{F2100W}$|Webb|MIRI|[`jwst_f2100w`](jwst_f2100w.par)|
|$\textrm{F2550W}$|Webb|MIRI|[`jwst_f2550w`](jwst_f2550w.par)|
|$\textrm{W1}$|WISE|WISE|[`wise_w1`](wise_w1.par)|
|$\textrm{W2}$|WISE|WISE|[`wise_w2`](wise_w2.par)|
|$\textrm{W3}$|WISE|WISE|[`wise_w3`](wise_w3.par)|
|$\textrm{W4}$|WISE|WISE|[`wise_w4`](wise_w4.par)|

## Ground-Based

|Filter|Telescope|Instrument|`sedpy` Name|
|:----:|:-------:|:--------:|:----------:|
|$J$|2MASS|2MASS|[`twomass_J`](twomass_J.par)|
|$H$|2MASS|2MASS|[`twomass_H`](twomass_H.par)|
|$Ks$|2MASS|2MASS|[`twomass_Ks`](twomass_Ks.par)|
|$u$|Blanco|DECam|[`decam_u`](decam_u.par)|
|$g$|Blanco|DECam|[`decam_g`](decam_g.par)|
|$r$|Blanco|DECam|[`decam_r`](decam_r.par)|
|$i$|Blanco|DECam|[`decam_i`](decam_i.par)|
|$z$|Blanco|DECam|[`decam_z`](decam_z.par)|
|$Y$|Blanco|DECam|[`decam_Y`](decam_Y.par)|
|$u_s\ \textrm{(9301)}$|CFHT|MegaCam|[`cfht_megacam_us_9301`](cfht_megacam_us_9301.par)|
|$g_s\ \textrm{(9401)}$|CFHT|MegaCam|[`cfht_megacam_gs_9401`](cfht_megacam_gs_9401.par)|
|$r_s\ \textrm{(9601)}$|CFHT|MegaCam|[`cfht_megacam_rs_9601`](cfht_megacam_rs_9601.par)|
|$i_s\ \textrm{(9701)}$|CFHT|MegaCam|[`cfht_megacam_is_9701`](cfht_megacam_is_9701.par)|
|$z_s\ \textrm{(9801)}$|CFHT|MegaCam|[`cfht_megacam_zs_9801`](cfht_megacam_zs_9801.par)|
|$J\ \textrm{(8101)}$|CFHT|WIRCam|[`cfht_wircam_J_8101`](cfht_wircam_J_8101.par)|
|$H\ \textrm{(8201)}$|CFHT|WIRCam|[`cfht_wircam_H_8201`](cfht_wircam_H_8201.par)|
|$K_s\ \textrm{(8302)}$|CFHT|WIRCam|[`cfht_wircam_Ks_8302`](cfht_wircam_Ks_8302.par)|
|$g$|Keck|LRIS|[`keck_lris_g`](keck_lris_g.par)|
|$R_s$|Keck|LRIS|[`keck_lris_Rs`](keck_lris_Rs.par)|
|$U\ \textrm{(k1001)}$|Mayall|Mosaic|[`mayall_mosaic_U_k1001`](mayall_mosaic_U_k1001.par)|
|$J_1$|Mayall|NEWFIRM|[`mayall_newfirm_J1`](mayall_newfirm_J1.par)|
|$J_2$|Mayall|NEWFIRM|[`mayall_newfirm_J2`](mayall_newfirm_J2.par)|
|$J_3$|Mayall|NEWFIRM|[`mayall_newfirm_J3`](mayall_newfirm_J3.par)|
|$H_1$|Mayall|NEWFIRM|[`mayall_newfirm_H1`](mayall_newfirm_H1.par)|
|$H_2$|Mayall|NEWFIRM|[`mayall_newfirm_H2`](mayall_newfirm_H2.par)|
|$K$|Mayall|NEWFIRM|[`mayall_newfirm_K`](mayall_newfirm_K.par)|
|$U_{38}\ \textrm{(ESO 841)}$|MPG/ESO|WFI|[`mpgeso_wfi_U38_eso841`](mpgeso_wfi_U38_eso841.par)|
|$B\ \textrm{(ESO 842)}$|MPG/ESO|WFI|[`mpgeso_wfi_B_eso842`](mpgeso_wfi_B_eso842.par)|
|$V\ \textrm{(ESO 843)}$|MPG/ESO|WFI|[`mpgeso_wfi_V_eso843`](mpgeso_wfi_V_eso843.par)|
|$R_c\ \textrm{(ESO 844)}$|MPG/ESO|WFI|[`mpgeso_wfi_Rc_eso844`](mpgeso_wfi_Rc_eso844.par)|
|$I_c\ \textrm{(ESO 845)}$|MPG/ESO|WFI|[`mpgeso_wfi_Ic_eso845`](mpgeso_wfi_Ic_eso845.par)|
|$u$|SDSS|SDSS|[`sdss_u`](sdss_u.par)|
|$g$|SDSS|SDSS|[`sdss_g`](sdss_g.par)|
|$r$|SDSS|SDSS|[`sdss_r`](sdss_r.par)|
|$i$|SDSS|SDSS|[`sdss_i`](sdss_i.par)|
|$z$|SDSS|SDSS|[`sdss_z`](sdss_z.par)|
|$J$|Subaru|MOIRCS|[`subaru_moircs_J`](subaru_moircs_J.par)|
|$H$|Subaru|MOIRCS|[`subaru_moircs_H`](subaru_moircs_H.par)|
|$K_s$|Subaru|MOIRCS|[`subaru_moircs_Ks`](subaru_moircs_Ks.par)|
|$B$|Subaru|Suprime-Cam|[`subaru_suprimecam_B`](subaru_suprimecam_B.par)|
|$V$|Subaru|Suprime-Cam|[`subaru_suprimecam_V`](subaru_suprimecam_V.par)|
|$R_c$|Subaru|Suprime-Cam|[`subaru_suprimecam_Rc`](subaru_suprimecam_Rc.par)|
|$r^\prime$|Subaru|Suprime-Cam|[`subaru_suprimecam_rp`](subaru_suprimecam_rp.par)|
|$i^\prime$|Subaru|Suprime-Cam|[`subaru_suprimecam_ip`](subaru_suprimecam_ip.par)|
|$z^\prime$|Subaru|Suprime-Cam|[`subaru_suprimecam_zp`](subaru_suprimecam_zp.par)|
|$\textrm{IA427}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia427`](subaru_suprimecam_ia427.par)|
|$\textrm{IA445}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia445`](subaru_suprimecam_ia445.par)|
|$\textrm{IA464}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia464`](subaru_suprimecam_ia464.par)|
|$\textrm{IA484}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia484`](subaru_suprimecam_ia484.par)|
|$\textrm{IA505}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia505`](subaru_suprimecam_ia505.par)|
|$\textrm{IA527}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia527`](subaru_suprimecam_ia527.par)|
|$\textrm{IA550}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia550`](subaru_suprimecam_ia550.par)|
|$\textrm{IA574}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia574`](subaru_suprimecam_ia574.par)|
|$\textrm{IA598}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia598`](subaru_suprimecam_ia598.par)|
|$\textrm{IA624}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia624`](subaru_suprimecam_ia624.par)|
|$\textrm{IA651}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia651`](subaru_suprimecam_ia651.par)|
|$\textrm{IA679}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia679`](subaru_suprimecam_ia679.par)|
|$\textrm{IA709}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia709`](subaru_suprimecam_ia709.par)|
|$\textrm{IA738}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia738`](subaru_suprimecam_ia738.par)|
|$\textrm{IA767}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia767`](subaru_suprimecam_ia767.par)|
|$\textrm{IA797}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia797`](subaru_suprimecam_ia797.par)|
|$\textrm{IA827}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia827`](subaru_suprimecam_ia827.par)|
|$\textrm{IA856}$|Subaru|Suprime-Cam|[`subaru_suprimecam_ia856`](subaru_suprimecam_ia856.par)|
|$J$|UKIRT|WFCAM|[`ukirt_wfcam_J`](ukirt_wfcam_J.par)|
|$H$|UKIRT|WFCAM|[`ukirt_wfcam_H`](ukirt_wfcam_H.par)|
|$K$|UKIRT|WFCAM|[`ukirt_wfcam_K`](ukirt_wfcam_K.par)|
|$Z$|VISTA|VIRCAM|[`vista_vircam_Z`](vista_vircam_Z.par)|
|$Y$|VISTA|VIRCAM|[`vista_vircam_Y`](vista_vircam_Y.par)|
|$J$|VISTA|VIRCAM|[`vista_vircam_J`](vista_vircam_J.par)|
|$H$|VISTA|VIRCAM|[`vista_vircam_H`](vista_vircam_H.par)|
|$K_s$|VISTA|VIRCAM|[`vista_vircam_Ks`](vista_vircam_Ks.par)|
|$J$|VLT|ISAAC|[`vlt_isaac_J`](vlt_isaac_J.par)|
|$H$|VLT|ISAAC|[`vlt_isaac_H`](vlt_isaac_H.par)|
|$K_s$|VLT|ISAAC|[`vlt_isaac_Ks`](vlt_isaac_Ks.par)|
|$U$|VLT|VIMOS|[`vlt_vimos_U`](vlt_vimos_U.par)|
|$R$|VLT|VIMOS|[`vlt_vimos_R`](vlt_vimos_R.par)|
