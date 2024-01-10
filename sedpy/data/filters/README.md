# Filters in `sedpy`

## Idealized Systems

### Bessell

|Filter|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:----------:|---------------------:|---------------------------:|:-------|
|$U$|[`bessell_U`](filters/bessell_U.par)|357.1 nm|52.4 nm||
|$B$|[`bessell_B`](filters/bessell_B.par)|434.4 nm|79.7 nm||
|$V$|[`bessell_V`](filters/bessell_V.par)|545.6 nm|80.7 nm||
|$R$|[`bessell_R`](filters/bessell_R.par)|644.2 nm|137.7 nm||
|$I$|[`bessell_I`](filters/bessell_I.par)|799.4 nm|107.9 nm||

### Stromgren

|Filter|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:----------:|---------------------:|---------------------------:|:-------|
|$u$|[`stromgren_u`](filters/stromgren_u.par)|344.4 nm|29.7 nm||
|$v$|[`stromgren_v`](filters/stromgren_v.par)|410.4 nm|23.1 nm||
|$b$|[`stromgren_b`](filters/stromgren_b.par)|466.7 nm|21.2 nm||
|$y$|[`stromgren_y`](filters/stromgren_y.par)|547.5 nm|23.9 nm||

## Space-Based Telescopes

### Gaia

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$G$||[`gaia_g`](filters/gaia_g.par)|590.2 nm|333.4 nm|DR2 revised|
|$G_\textrm{BP}$|BP|[`gaia_bp`](filters/gaia_bp.par)|488.8 nm|211.7 nm|DR2 revised|
|$G_\textrm{RP}$|RP|[`gaia_rp`](filters/gaia_rp.par)|762.7 nm|216.9 nm|DR2 revised|

### Galaxy Evolution Explorer (GALEX)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{FUV}$|FUV|[`galex_FUV`](filters/galex_FUV.par)|152.8 nm|242.2 nm||
|$\textrm{NUV}$|NUV|[`galex_NUV`](filters/galex_NUV.par)|227.1 nm|612.5 nm||

### Herschel Space Observatory

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$70\ \textrm{μm}$|PACS|[`herschel_pacs_70`](filters/herschel_pacs_70.par)|70.02 μm|16.94 μm||
|$100\ \textrm{μm}$|PACS|[`herschel_pacs_100`](filters/herschel_pacs_100.par)|99.60 μm|25.41 μm||
|$160\ \textrm{μm}$|PACS|[`herschel_pacs_160`](filters/herschel_pacs_160.par)|158.6 μm|53.13 μm||
|$250\ \textrm{μm}$|SPIRE|[`herschel_spire_250`](filters/herschel_spire_250.par)|246.1 μm|54.32 μm||
|$350\ \textrm{μm}$|SPIRE|[`herschel_spire_350`](filters/herschel_spire_350.par)|345.5 μm|74.97 μm||
|$500\ \textrm{μm}$|SPIRE|[`herschel_spire_500`](filters/herschel_spire_500.par)|493.1 μm|139.6 μm||

### Hubble Space Telescope (HST)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{F435W}$|ACS (WFC)|[`acs_wfc_f435w`](filters/acs_wfc_f435w.par)|430.7 nm|70.22 nm||
|$\textrm{F475W}$|ACS (WFC)|[`acs_wfc_f475w`](filters/acs_wfc_f475w.par)|470.8 nm|98.93 nm||
|$\textrm{F555W}$|ACS (WFC)|[`acs_wfc_f555w`](filters/acs_wfc_f555w.par)|533.6 nm|84.78 nm||
|$\textrm{F606W}$|ACS (WFC)|[`acs_wfc_f606w`](filters/acs_wfc_f606w.par)|584.2 nm|158.3 nm||
|$\textrm{F625W}$|ACS (WFC)|[`acs_wfc_f625w`](filters/acs_wfc_f625w.par)|628.4 nm|97.82 nm||
|$\textrm{F775W}$|ACS (WFC)|[`acs_wfc_f775w`](filters/acs_wfc_f775w.par)|766.9 nm|102.3 nm||
|$\textrm{F814W}$|ACS (WFC)|[`acs_wfc_f814w`](filters/acs_wfc_f814w.par)|800.6 nm|153.8 nm||
|$\textrm{F850LP}$|ACS (WFC)|[`acs_wfc_f850lp`](filters/acs_wfc_f850lp.par)|900.5 nm|124.1 nm||
|$\textrm{F225W}$|WFC3 (UVIS)|[`wfc3_uvis_f225w`](filters/wfc3_uvis_f225w.par)|235.9 nm|41.8 nm||
|$\textrm{F275W}$|WFC3 (UVIS)|[`wfc3_uvis_f275w`](filters/wfc3_uvis_f275w.par)|270.0 nm|38.71 nm||
|$\textrm{F336W}$|WFC3 (UVIS)|[`wfc3_uvis_f336w`](filters/wfc3_uvis_f336w.par)|334.7 nm|37.31 nm||
|$\textrm{F390W}$|WFC3 (UVIS)|[`wfc3_uvis_f390w`](filters/wfc3_uvis_f390w.par)|390.1 nm|68.57 nm|Average of UVIS1 and UVIS2|
|$\textrm{F475W}$|WFC3 (UVIS)|[`wfc3_uvis_f475w`](filters/wfc3_uvis_f475w.par)|473.6 nm|99.19 nm||
|$\textrm{F555W}$|WFC3 (UVIS)|[`wfc3_uvis_f555w`](filters/wfc3_uvis_f555w.par)|525.6 nm|121.8 nm||
|$\textrm{F606W}$|WFC3 (UVIS)|[`wfc3_uvis_f606w`](filters/wfc3_uvis_f606w.par)|581.3 nm|154.6 nm||
|$\textrm{F814W}$|WFC3 (UVIS)|[`wfc3_uvis_f814w`](filters/wfc3_uvis_f814w.par)|797.3 nm|156.2 nm||
|$\textrm{F098M}$|WFC3 (IR)|[`wfc3_ir_f098m`](filters/wfc3_ir_f098m.par)|0.984 μm|117.8 nm||
|$\textrm{F105W}$|WFC3 (IR)|[`wfc3_ir_f105w`](filters/wfc3_ir_f105w.par)|1.048 μm|199.1 nm||
|$\textrm{F110W}$|WFC3 (IR)|[`wfc3_ir_f110w`](filters/wfc3_ir_f110w.par)|1.136 μm|336.4 nm||
|$\textrm{F125W}$|WFC3 (IR)|[`wfc3_ir_f125w`](filters/wfc3_ir_f125w.par)|1.241 μm|204.2 nm||
|$\textrm{F140W}$|WFC3 (IR)|[`wfc3_ir_f140w`](filters/wfc3_ir_f140w.par)|1.383 μm|267.6 nm||
|$\textrm{F160W}$|WFC3 (IR)|[`wfc3_ir_f160w`](filters/wfc3_ir_f160w.par)|1.532 μm|194.5 nm||

### James Webb Space Telescope (JWST)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{F070W}$|NIRCam|[`jwst_f070w`](filters/jwst_f070w.par)|700.7 nm|116.4 nm||
|$\textrm{F090W}$|NIRCam|[`jwst_f090w`](filters/jwst_f090w.par)|898.1 nm|144.8 nm||
|$\textrm{F115W}$|NIRCam|[`jwst_f115w`](filters/jwst_f115w.par)|1.149 μm|188.0 nm||
|$\textrm{F150W}$|NIRCam|[`jwst_f150w`](filters/jwst_f150w.par)|1.494 μm|232.1 nm||
|$\textrm{F200W}$|NIRCam|[`jwst_f200w`](filters/jwst_f200w.par)|1.979 μm|322.6 nm||
|$\textrm{F277W}$|NIRCam|[`jwst_f277w`](filters/jwst_f277w.par)|2.747 μm|480.3 nm||
|$\textrm{F356W}$|NIRCam|[`jwst_f356w`](filters/jwst_f356w.par)|3.552 μm|577.4 nm||
|$\textrm{F444W}$|NIRCam|[`jwst_f444w`](filters/jwst_f444w.par)|4.381 μm|750.8 nm||
|$\textrm{F140M}$|NIRCam|[`jwst_f140m`](filters/jwst_f140m.par)|1.404 μm|105.4 nm||
|$\textrm{F162M}$|NIRCam|[`jwst_f162m`](filters/jwst_f162m.par)|1.626 μm|120.5 nm||
|$\textrm{F182M}$|NIRCam|[`jwst_f182m`](filters/jwst_f182m.par)|1.842 μm|170.5 nm||
|$\textrm{F210M}$|NIRCam|[`jwst_f210m`](filters/jwst_f210m.par)|2.094 μm|146.4 nm||
|$\textrm{F250M}$|NIRCam|[`jwst_f250m`](filters/jwst_f250m.par)|2.502 μm|129.3 nm||
|$\textrm{F300M}$|NIRCam|[`jwst_f300m`](filters/jwst_f300m.par)|2.986 μm|233.7 nm||
|$\textrm{F335M}$|NIRCam|[`jwst_f335m`](filters/jwst_f335m.par)|3.359 μm|261.1 nm||
|$\textrm{F360M}$|NIRCam|[`jwst_f360m`](filters/jwst_f360m.par)|3.620 μm|275.9 nm||
|$\textrm{F410M}$|NIRCam|[`jwst_f410m`](filters/jwst_f410m.par)|4.078 μm|311.5 nm||
|$\textrm{F430M}$|NIRCam|[`jwst_f430m`](filters/jwst_f430m.par)|4.280 μm|164.8 nm||
|$\textrm{F460M}$|NIRCam|[`jwst_f460m`](filters/jwst_f460m.par)|4.628 μm|182.9 nm||
|$\textrm{F480M}$|NIRCam|[`jwst_f480m`](filters/jwst_f480m.par)|4.815 μm|256.1 nm||
|$\textrm{F070W}$|NIRCam|[`jwst_f150w2`](filters/jwst_f150w2.par)| 1.578 μm | 908.6 nm ||
|$\textrm{F090W}$|NIRCam|[`jwst_f322w2`](filters/jwst_f322w2.par)| 3.179 μm | 1099.2 nm ||
|$\textrm{F560W}$|MIRI|[`jwst_f560w`](filters/jwst_f560w.par)|5.615 μm|791.9 nm||
|$\textrm{F770W}$|MIRI|[`jwst_f770w`](filters/jwst_f770w.par)|7.591 μm|1.426 μm||
|$\textrm{F1000W}$|MIRI|[`jwst_f1000w`](filters/jwst_f1000w.par)|9.923 μm|1.287 μm||
|$\textrm{F1130W}$|MIRI|[`jwst_f1130w`](filters/jwst_f1130w.par)|11.30 μm|558.2 nm||
|$\textrm{F1280W}$|MIRI|[`jwst_f1280w`](filters/jwst_f1280w.par)|12.77 μm|1.732 μm||
|$\textrm{F1500W}$|MIRI|[`jwst_f1500w`](filters/jwst_f1500w.par)|15.01 μm|2.151 μm||
|$\textrm{F1800W}$|MIRI|[`jwst_f1800w`](filters/jwst_f1800w.par)|17.94 μm|2.105 μm||
|$\textrm{F2100W}$|MIRI|[`jwst_f2100w`](filters/jwst_f2100w.par)|20.70 μm|3.296 μm||
|$\textrm{F2550W}$|MIRI|[`jwst_f2550w`](filters/jwst_f2550w.par)|25.28 μm|3.461 μm||

The default NIRCam filters are averages over the A and B modules. For every
JWST/NIRCam filter there are also module specific curves available as
`jwst_mod[a,b]_f*w.par`.  For the SW filters these correspond to detector '1' of
each module.

### Spitzer Space Telescope

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{Band}\ 1$|IRAC|[`spitzer_irac_ch1`](filters/spitzer_irac_ch1.par)|3.538 μm|504.5 nm||
|$\textrm{Band}\ 2$|IRAC|[`spitzer_irac_ch2`](filters/spitzer_irac_ch2.par)|4.478 μm|667.1 nm||
|$\textrm{Band}\ 3$|IRAC|[`spitzer_irac_ch3`](filters/spitzer_irac_ch3.par)|5.696 μm|946.9 nm||
|$\textrm{Band}\ 4$|IRAC|[`spitzer_irac_ch4`](filters/spitzer_irac_ch4.par)|7.798 μm|1.937 μm||
|$24\ \textrm{μm}$|MIPS|[`spitzer_mips_24`](filters/spitzer_mips_24.par)|23.43 μm|4.497 μm||
|$70\ \textrm{μm}$|MIPS|[`spitzer_mips_70`](filters/spitzer_mips_70.par)|69.86 μm|19.50 μm||
|$160\ \textrm{μm}$|MIPS|[`spitzer_mips_160`](filters/spitzer_mips_160.par)|154.3 μm|30.53 μm||

### Wide-field Infrared Survey Explorer (WISE)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{W1}$|WISE|[`wise_w1`](filters/wise_w1.par)|3.346 μm|647.9 nm||
|$\textrm{W2}$|WISE|[`wise_w2`](filters/wise_w2.par)|4.595 μm|759.7 nm||
|$\textrm{W3}$|WISE|[`wise_w3`](filters/wise_w3.par)|11.55 μm|5.851 μm||
|$\textrm{W4}$|WISE|[`wise_w4`](filters/wise_w4.par)|22.08 μm|3.720 μm||

### Hipparcos

|Filter|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:----------:|---------------------:|---------------------------:|:-------|
|$H$|[`hipparcos_H`](filters/hipparcos_H.par)|519.0 nm|198.5 nm||
|$B$|[`hipparcos_B`](filters/hipparcos_B.par)|413.1 nm|71.4 nm||
|$V$|[`hipparcos_V`](filters/hipparcos_V.par)|527.8 nm|91.2 nm||

### SWIFT/UVOT

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$w1$|UVOT|[`uvot_w1`](filters/uvot_w1.par)|251.7 nm|93.8 nm||
|$w2$|UVOT|[`uvot_w2`](filters/uvot_w2.par)|201.0 nm|67.2 nm||
|$m2$|UVOT|[`uvot_m2`](filters/uvot_m2.par)|223.0 nm|45.9 nm||

## Ground-Based Telescopes

### Two Micron All-Sky Survey (2MASS) Telescope

|Filter|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:----------:|---------------------:|---------------------------:|:-------|
|$J$|[`twomass_J`](filters/twomass_J.par)|1.232 μm|156.4 nm||
|$H$|[`twomass_H`](filters/twomass_H.par)|1.642 μm|184.7 nm||
|$Ks$|[`twomass_Ks`](filters/twomass_Ks.par)|2.157 μm|206.9 nm||

### Víctor M. Blanco Telescope

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$u$|DECam|[`decam_u`](filters/decam_u.par)|380.9 nm|33.07 nm||
|$g$|DECam|[`decam_g`](filters/decam_g.par)|477.1 nm|99.45 nm||
|$r$|DECam|[`decam_r`](filters/decam_r.par)|639.2 nm|101.5 nm||
|$i$|DECam|[`decam_i`](filters/decam_i.par)|778.9 nm|102.0 nm||
|$z$|DECam|[`decam_z`](filters/decam_z.par)|914.9 nm|101.4 nm||
|$Y$|DECam|[`decam_Y`](filters/decam_Y.par)|988.5 nm|62.56 nm||

### Canada–France–Hawaii Telescope (CFHT)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$u_s\ \textrm{(9301)}$|MegaCam|[`cfht_megacam_us_9301`](filters/cfht_megacam_us_9301.par)|381.3 nm|55.39 nm||
|$g_s\ \textrm{(9401)}$|MegaCam|[`cfht_megacam_gs_9401`](filters/cfht_megacam_gs_9401.par)|483.3 nm|99.01 nm||
|$r_s\ \textrm{(9601)}$|MegaCam|[`cfht_megacam_rs_9601`](filters/cfht_megacam_rs_9601.par)|622.4 nm|86.20 nm||
|$i_s\ \textrm{(9701)}$|MegaCam|[`cfht_megacam_is_9701`](filters/cfht_megacam_is_9701.par)|765.1 nm|102.0 nm||
|$z_s\ \textrm{(9801)}$|MegaCam|[`cfht_megacam_zs_9801`](filters/cfht_megacam_zs_9801.par)|884.7 nm|108.7 nm||
|$J\ \textrm{(8101)}$|WIRCam|[`cfht_wircam_J_8101`](filters/cfht_wircam_J_8101.par)|1.251 μm|108.0 nm||
|$H\ \textrm{(8201)}$|WIRCam|[`cfht_wircam_H_8201`](filters/cfht_wircam_H_8201.par)|1.625 μm|196.5 nm||
|$K_s\ \textrm{(8302)}$|WIRCam|[`cfht_wircam_Ks_8302`](filters/cfht_wircam_Ks_8302.par)|2.154 μm|213.8 nm||

### W. M. Keck Telescope I (Keck I)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$g$|LRIS|[`keck_lris_g`](filters/keck_lris_g.par)|473.3 nm|66.87 nm||
|$R_s$|LRIS|[`keck_lris_Rs`](filters/keck_lris_Rs.par)|678.7 nm|103.6 nm||

### Nicholas U. Mayall Telescope

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$U\ \textrm{(k1001)}$|Mosaic|[`mayall_mosaic_U_k1001`](filters/mayall_mosaic_U_k1001.par)|358.2 nm|46.70 nm||
|$J_1$|NEWFIRM|[`mayall_newfirm_J1`](filters/mayall_newfirm_J1.par)|1.044 μm|102.8 nm||
|$J_2$|NEWFIRM|[`mayall_newfirm_J2`](filters/mayall_newfirm_J2.par)|1.193 μm|103.6 nm||
|$J_3$|NEWFIRM|[`mayall_newfirm_J3`](filters/mayall_newfirm_J3.par)|1.276 μm|98.84 nm||
|$H_1$|NEWFIRM|[`mayall_newfirm_H1`](filters/mayall_newfirm_H1.par)|1.559 μm|115.9 nm||
|$H_2$|NEWFIRM|[`mayall_newfirm_H2`](filters/mayall_newfirm_H2.par)|1.705 μm|120.9 nm||
|$K$|NEWFIRM|[`mayall_newfirm_K`](filters/mayall_newfirm_K.par)|2.164 μm|219.6 nm||

### MPG/ESO Telescope

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$U_{38}\ \textrm{(ESO 841)}$|WFI|[`mpgeso_wfi_U38_eso841`](filters/mpgeso_wfi_U38_eso841.par)|368.2 nm|29.64 nm||
|$B\ \textrm{(ESO 842)}$|WFI|[`mpgeso_wfi_B_eso842`](filters/mpgeso_wfi_B_eso842.par)|457.3 nm|68.16 nm||
|$V\ \textrm{(ESO 843)}$|WFI|[`mpgeso_wfi_V_eso843`](filters/mpgeso_wfi_V_eso843.par)|536.0 nm|61.65 nm||
|$R_c\ \textrm{(ESO 844)}$|WFI|[`mpgeso_wfi_Rc_eso844`](filters/mpgeso_wfi_Rc_eso844.par)|646.5 nm|110.6 nm||
|$I_c\ \textrm{(ESO 845)}$|WFI|[`mpgeso_wfi_Ic_eso845`](filters/mpgeso_wfi_Ic_eso845.par)|859.5 nm|139.9 nm||

### Sloan Digital Sky Survey (SDSS) Telescope

|Filter|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:----------:|---------------------:|---------------------------:|:-------|
|$u$|[`sdss_u0`](filters/sdss_u0.par)|354.6 nm|45.72 nm||
|$g$|[`sdss_g0`](filters/sdss_g0.par)|467.0 nm|92.79 nm||
|$r$|[`sdss_r0`](filters/sdss_r0.par)|615.6 nm|81.28 nm||
|$i$|[`sdss_i0`](filters/sdss_i0.par)|747.2 nm|89.08 nm||
|$z$|[`sdss_z0`](filters/sdss_z0.par)|891.7 nm|118.3 nm||

### Subaru Telescope

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$J$|MOIRCS|[`subaru_moircs_J`](filters/subaru_moircs_J.par)|1.250 μm|110.9 nm||
|$H$|MOIRCS|[`subaru_moircs_H`](filters/subaru_moircs_H.par)|1.631 μm|191.3 nm||
|$K_s$|MOIRCS|[`subaru_moircs_Ks`](filters/subaru_moircs_Ks.par)|2.154 μm|207.1 nm||
|$B$|Suprime-Cam|[`subaru_suprimecam_B`](filters/subaru_suprimecam_B.par)|442.7 nm|71.91 nm||
|$V$|Suprime-Cam|[`subaru_suprimecam_V`](filters/subaru_suprimecam_V.par)|545.5 nm|68.03 nm||
|$R_c$|Suprime-Cam|[`subaru_suprimecam_Rc`](filters/subaru_suprimecam_Rc.par)|649.0 nm|81.83 nm||
|$r^\prime$|Suprime-Cam|[`subaru_suprimecam_rp`](filters/subaru_suprimecam_rp.par)|624.9 nm|95.78 nm||
|$i^\prime$|Suprime-Cam|[`subaru_suprimecam_ip`](filters/subaru_suprimecam_ip.par)|764.6 nm|103.3 nm||
|$z^\prime$|Suprime-Cam|[`subaru_suprimecam_zp`](filters/subaru_suprimecam_zp.par)|901.1 nm|92.26 nm||
|$\textrm{IA427}$|Suprime-Cam|[`subaru_suprimecam_ia427`](filters/subaru_suprimecam_ia427.par)|425.9 nm|15.13 nm||
|$\textrm{IA445}$|Suprime-Cam|[`subaru_suprimecam_ia445`](filters/subaru_suprimecam_ia445.par)|444.2 nm|15.00 nm||
|$\textrm{IA464}$|Suprime-Cam|[`subaru_suprimecam_ia464`](filters/subaru_suprimecam_ia464.par)|463.2 nm|15.81 nm||
|$\textrm{IA484}$|Suprime-Cam|[`subaru_suprimecam_ia484`](filters/subaru_suprimecam_ia484.par)|484.6 nm|16.83 nm||
|$\textrm{IA505}$|Suprime-Cam|[`subaru_suprimecam_ia505`](filters/subaru_suprimecam_ia505.par)|506.0 nm|17.22 nm||
|$\textrm{IA527}$|Suprime-Cam|[`subaru_suprimecam_ia527`](filters/subaru_suprimecam_ia527.par)|525.8 nm|18.49 nm||
|$\textrm{IA550}$|Suprime-Cam|[`subaru_suprimecam_ia550`](filters/subaru_suprimecam_ia550.par)|549.4 nm|20.35 nm||
|$\textrm{IA574}$|Suprime-Cam|[`subaru_suprimecam_ia574`](filters/subaru_suprimecam_ia574.par)|576.2 nm|20.49 nm||
|$\textrm{IA598}$|Suprime-Cam|[`subaru_suprimecam_ia598`](filters/subaru_suprimecam_ia598.par)|600.6 nm|22.20 nm||
|$\textrm{IA624}$|Suprime-Cam|[`subaru_suprimecam_ia624`](filters/subaru_suprimecam_ia624.par)|622.9 nm|22.55 nm||
|$\textrm{IA651}$|Suprime-Cam|[`subaru_suprimecam_ia651`](filters/subaru_suprimecam_ia651.par)|649.7 nm|24.03 nm||
|$\textrm{IA679}$|Suprime-Cam|[`subaru_suprimecam_ia679`](filters/subaru_suprimecam_ia679.par)|678.0 nm|25.04 nm||
|$\textrm{IA709}$|Suprime-Cam|[`subaru_suprimecam_ia709`](filters/subaru_suprimecam_ia709.par)|707.2 nm|23.69 nm||
|$\textrm{IA738}$|Suprime-Cam|[`subaru_suprimecam_ia738`](filters/subaru_suprimecam_ia738.par)|735.8 nm|23.74 nm||
|$\textrm{IA767}$|Suprime-Cam|[`subaru_suprimecam_ia767`](filters/subaru_suprimecam_ia767.par)|767.9 nm|26.27 nm||
|$\textrm{IA797}$|Suprime-Cam|[`subaru_suprimecam_ia797`](filters/subaru_suprimecam_ia797.par)|796.4 nm|26.61 nm||
|$\textrm{IA827}$|Suprime-Cam|[`subaru_suprimecam_ia827`](filters/subaru_suprimecam_ia827.par)|824.5 nm|24.62 nm||
|$\textrm{IA856}$|Suprime-Cam|[`subaru_suprimecam_ia856`](filters/subaru_suprimecam_ia856.par)|856.3 nm|26.78 nm||
|$g$|Hyper Suprime-Cam|[`hsc_g`](filters/hsc_g.par)|475.5 nm|97.7 nm||
|$r$|Hyper Suprime-Cam|[`hsc_r`](filters/hsc_r.par)|618.4 nm|101.4 nm||
|$i$|Hyper Suprime-Cam|[`hsc_i`](filters/hsc_i.par)|766.1 nm|107.9 nm||
|$z$|Hyper Suprime-Cam|[`hsc_z`](filters/hsc_z.par)|889.7 nm|55.2 nm||
|$y$|Hyper Suprime-Cam|[`hsc_y`](filters/hsc_y.par)|976.2 nm|73.1 nm||

### United Kingdom Infrared Telescope (UKIRT)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$J$|WFCAM|[`ukirt_wfcam_J`](filters/ukirt_wfcam_J.par)|1.248 μm|112.2 nm||
|$H$|WFCAM|[`ukirt_wfcam_H`](filters/ukirt_wfcam_H.par)|1.631 μm|204.6 nm||
|$K$|WFCAM|[`ukirt_wfcam_K`](filters/ukirt_wfcam_K.par)|2.201 μm|246.5 nm||

### Visible and Infrared Survey Telescope for Astronomy (VISTA)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$Z$|VIRCAM|[`vista_vircam_Z`](filters/vista_vircam_Z.par)|876.6 nm|68.04 nm||
|$Y$|VIRCAM|[`vista_vircam_Y`](filters/vista_vircam_Y.par)|1.020 μm|75.26 nm||
|$J$|VIRCAM|[`vista_vircam_J`](filters/vista_vircam_J.par)|1.250 μm|120.9 nm||
|$H$|VIRCAM|[`vista_vircam_H`](filters/vista_vircam_H.par)|1.639 μm|200.6 nm||
|$K_s$|VIRCAM|[`vista_vircam_Ks`](filters/vista_vircam_Ks.par)|2.145 μm|258.7 nm||

### Very Large Telescope (VLT) Array

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$J$|ISAAC|[`vlt_isaac_J`](filters/vlt_isaac_J.par)|1.231 μm|178.0 nm||
|$H$|ISAAC|[`vlt_isaac_H`](filters/vlt_isaac_H.par)|1.645 μm|199.2 nm||
|$K_s$|ISAAC|[`vlt_isaac_Ks`](filters/vlt_isaac_Ks.par)|2.164 μm|185.5 nm||
|$U$|VIMOS|[`vlt_vimos_U`](filters/vlt_vimos_U.par)|374.2 nm|39.07 nm||
|$R$|VIMOS|[`vlt_vimos_R`](filters/vlt_vimos_R.par)|641.7 nm|94.24 nm||

### BASS and MzLS (https://legacysurvey.org)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$g$|90-Prime|[`BASS-g`](filters/BASS-g.par)|472.1 nm|100.2 nm||
|$r$|90-Prime|[`BASS-r`](filters/BASS-r.par)|637.0 nm|97.76 nm||
|$z$|Mosaic-3|[`MzLS-z`](filters/MzLS-z.par)|917.9 nm|100.7 nm||
