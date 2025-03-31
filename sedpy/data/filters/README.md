# Filters in `sedpy`

## Idealized Systems

### Bessell

|Filter|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:----------:|---------------------:|---------------------------:|:-------|
|$U$|[`bessell_U`](./bessell_U.par)|357.1 nm|52.4 nm||
|$B$|[`bessell_B`](./bessell_B.par)|434.4 nm|79.7 nm||
|$V$|[`bessell_V`](./bessell_V.par)|545.6 nm|80.7 nm||
|$R$|[`bessell_R`](./bessell_R.par)|644.2 nm|137.7 nm||
|$I$|[`bessell_I`](./bessell_I.par)|799.4 nm|107.9 nm||

### Stromgren

|Filter|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:----------:|---------------------:|---------------------------:|:-------|
|$u$|[`stromgren_u`](./stromgren_u.par)|344.4 nm|29.7 nm||
|$v$|[`stromgren_v`](./stromgren_v.par)|410.4 nm|23.1 nm||
|$b$|[`stromgren_b`](./stromgren_b.par)|466.7 nm|21.2 nm||
|$y$|[`stromgren_y`](./stromgren_y.par)|547.5 nm|23.9 nm||

## Space-Based Telescopes

### Gaia

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$G$||[`gaia_g`](./gaia_g.par)|588.2 nm|332 nm|EDR3 (Riello21)|
|$G_\textrm{BP}$|BP|[`gaia_bp`](./gaia_bp.par)|496.4 nm|202.4 nm|EDR3 (Riello21)|
|$G_\textrm{RP}$|RP|[`gaia_rp`](./gaia_rp.par)|765.9 nm|214.3 nm|EDR3 (Riello21)|

### Galaxy Evolution Explorer (GALEX)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{FUV}$|FUV|[`galex_FUV`](./galex_FUV.par)|152.8 nm|242.2 nm||
|$\textrm{NUV}$|NUV|[`galex_NUV`](./galex_NUV.par)|227.1 nm|612.5 nm||

### Herschel Space Observatory

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$70\ \textrm{μm}$|PACS|[`herschel_pacs_70`](./herschel_pacs_70.par)|70.02 μm|16.94 μm||
|$100\ \textrm{μm}$|PACS|[`herschel_pacs_100`](./herschel_pacs_100.par)|99.60 μm|25.41 μm||
|$160\ \textrm{μm}$|PACS|[`herschel_pacs_160`](./herschel_pacs_160.par)|158.6 μm|53.13 μm||
|$250\ \textrm{μm}$|SPIRE|[`herschel_spire_250`](./herschel_spire_250.par)|246.1 μm|54.32 μm||
|$350\ \textrm{μm}$|SPIRE|[`herschel_spire_350`](./herschel_spire_350.par)|345.5 μm|74.97 μm||
|$500\ \textrm{μm}$|SPIRE|[`herschel_spire_500`](./herschel_spire_500.par)|493.1 μm|139.6 μm||

### Hipparcos

|Filter|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:----------:|---------------------:|---------------------------:|:-------|
|$H$|[`hipparcos_H`](./hipparcos_H.par)|519.0 nm|198.5 nm||
|$B$|[`hipparcos_B`](./hipparcos_B.par)|413.1 nm|71.4 nm||
|$V$|[`hipparcos_V`](./hipparcos_V.par)|527.8 nm|91.2 nm||

### Hubble Space Telescope (HST)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{F435W}$|ACS (WFC)|[`acs_wfc_f435w`](./acs_wfc_f435w.par)|430.7 nm|70.22 nm||
|$\textrm{F475W}$|ACS (WFC)|[`acs_wfc_f475w`](./acs_wfc_f475w.par)|470.8 nm|98.93 nm||
|$\textrm{F555W}$|ACS (WFC)|[`acs_wfc_f555w`](./acs_wfc_f555w.par)|533.6 nm|84.78 nm||
|$\textrm{F606W}$|ACS (WFC)|[`acs_wfc_f606w`](./acs_wfc_f606w.par)|584.2 nm|158.3 nm||
|$\textrm{F625W}$|ACS (WFC)|[`acs_wfc_f625w`](./acs_wfc_f625w.par)|628.4 nm|97.82 nm||
|$\textrm{F775W}$|ACS (WFC)|[`acs_wfc_f775w`](./acs_wfc_f775w.par)|766.9 nm|102.3 nm||
|$\textrm{F814W}$|ACS (WFC)|[`acs_wfc_f814w`](./acs_wfc_f814w.par)|800.6 nm|153.8 nm||
|$\textrm{F850LP}$|ACS (WFC)|[`acs_wfc_f850lp`](./acs_wfc_f850lp.par)|900.5 nm|124.1 nm||
|$\textrm{F225W}$|WFC3 (UVIS)|[`wfc3_uvis_f225w`](./wfc3_uvis_f225w.par)|235.9 nm|41.8 nm||
|$\textrm{F275W}$|WFC3 (UVIS)|[`wfc3_uvis_f275w`](./wfc3_uvis_f275w.par)|270.0 nm|38.71 nm||
|$\textrm{F336W}$|WFC3 (UVIS)|[`wfc3_uvis_f336w`](./wfc3_uvis_f336w.par)|334.7 nm|37.31 nm||
|$\textrm{F390W}$|WFC3 (UVIS)|[`wfc3_uvis_f390w`](./wfc3_uvis_f390w.par)|390.1 nm|68.57 nm|Average of UVIS1 and UVIS2|
|$\textrm{F475W}$|WFC3 (UVIS)|[`wfc3_uvis_f475w`](./wfc3_uvis_f475w.par)|473.6 nm|99.19 nm||
|$\textrm{F555W}$|WFC3 (UVIS)|[`wfc3_uvis_f555w`](./wfc3_uvis_f555w.par)|525.6 nm|121.8 nm||
|$\textrm{F606W}$|WFC3 (UVIS)|[`wfc3_uvis_f606w`](./wfc3_uvis_f606w.par)|581.3 nm|154.6 nm||
|$\textrm{F814W}$|WFC3 (UVIS)|[`wfc3_uvis_f814w`](./wfc3_uvis_f814w.par)|797.3 nm|156.2 nm||
|$\textrm{F098M}$|WFC3 (IR)|[`wfc3_ir_f098m`](./wfc3_ir_f098m.par)|0.984 μm|117.8 nm||
|$\textrm{F105W}$|WFC3 (IR)|[`wfc3_ir_f105w`](./wfc3_ir_f105w.par)|1.048 μm|199.1 nm||
|$\textrm{F110W}$|WFC3 (IR)|[`wfc3_ir_f110w`](./wfc3_ir_f110w.par)|1.136 μm|336.4 nm||
|$\textrm{F125W}$|WFC3 (IR)|[`wfc3_ir_f125w`](./wfc3_ir_f125w.par)|1.241 μm|204.2 nm||
|$\textrm{F140W}$|WFC3 (IR)|[`wfc3_ir_f140w`](./wfc3_ir_f140w.par)|1.383 μm|267.6 nm||
|$\textrm{F160W}$|WFC3 (IR)|[`wfc3_ir_f160w`](./wfc3_ir_f160w.par)|1.532 μm|194.5 nm||

### James Webb Space Telescope (JWST)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{F070W}$||[`jwst_f070w`](./jwst_f070w.par)|701.4 nm|114.7 nm|May 2024|
|$\textrm{F090W}$||[`jwst_f090w`](./jwst_f090w.par)|898.3 nm|143.9 nm|May 2024|
|$\textrm{F115W}$||[`jwst_f115w`](./jwst_f115w.par)|1148.6 nm|188.0 nm|May 2024|
|$\textrm{F150W}$||[`jwst_f150w`](./jwst_f150w.par)|1494.4 nm|232.1 nm|May 2024|
|$\textrm{F200W}$||[`jwst_f200w`](./jwst_f200w.par)|1978.2 nm|321.9 nm|May 2024|
|$\textrm{F277W}$||[`jwst_f277w`](./jwst_f277w.par)|2760.5 nm|490.3 nm|May 2024|
|$\textrm{F356W}$||[`jwst_f356w`](./jwst_f356w.par)|3548.4 nm|581.0 nm|May 2024|
|$\textrm{F444W}$||[`jwst_f444w`](./jwst_f444w.par)|4378.7 nm|745.5 nm|May 2024|
|$\textrm{F140M}$||[`jwst_f140m`](./jwst_f140m.par)|1403.9 nm|105.1 nm|May 2024|
|$\textrm{F162M}$||[`jwst_f162m`](./jwst_f162m.par)|1625.6 nm|120.5 nm|May 2024|
|$\textrm{F182M}$||[`jwst_f182m`](./jwst_f182m.par)|1842.4 nm|170.5 nm|May 2024|
|$\textrm{F210M}$||[`jwst_f210m`](./jwst_f210m.par)|2093.8 nm|146.3 nm|May 2024|
|$\textrm{F250M}$||[`jwst_f250m`](./jwst_f250m.par)|2502.0 nm|125.9 nm|May 2024|
|$\textrm{F300M}$||[`jwst_f300m`](./jwst_f300m.par)|2992.9 nm|231.7 nm|May 2024|
|$\textrm{F335M}$||[`jwst_f335m`](./jwst_f335m.par)|3358.6 nm|260.6 nm|May 2024|
|$\textrm{F360M}$||[`jwst_f360m`](./jwst_f360m.par)|3619.5 nm|276.0 nm|May 2024|
|$\textrm{F410M}$||[`jwst_f410m`](./jwst_f410m.par)|4078.7 nm|311.0 nm|May 2024|
|$\textrm{F430M}$||[`jwst_f430m`](./jwst_f430m.par)|4280.2 nm|162.9 nm|May 2024|
|$\textrm{F460M}$||[`jwst_f460m`](./jwst_f460m.par)|4628.5 nm|166.7 nm|May 2024|
|$\textrm{F480M}$||[`jwst_f480m`](./jwst_f480m.par)|4814.5 nm|235.1 nm|May 2024|
|$\textrm{F150W2}$||[`jwst_f150w2`](./jwst_f150w2.par)|1578.7 nm|909.4 nm|May 2024|
|$\textrm{F322W2}$||[`jwst_f322w2`](./jwst_f322w2.par)|3179.2 nm|1099.3 nm|May 2024|
|$\textrm{F560W}$|MIRI|[`jwst_f560w`](./jwst_f560w.par)|5.615 μm|791.9 nm||
|$\textrm{F770W}$|MIRI|[`jwst_f770w`](./jwst_f770w.par)|7.591 μm|1.426 μm||
|$\textrm{F1000W}$|MIRI|[`jwst_f1000w`](./jwst_f1000w.par)|9.923 μm|1.287 μm||
|$\textrm{F1130W}$|MIRI|[`jwst_f1130w`](./jwst_f1130w.par)|11.30 μm|558.2 nm||
|$\textrm{F1280W}$|MIRI|[`jwst_f1280w`](./jwst_f1280w.par)|12.77 μm|1.732 μm||
|$\textrm{F1500W}$|MIRI|[`jwst_f1500w`](./jwst_f1500w.par)|15.01 μm|2.151 μm||
|$\textrm{F1800W}$|MIRI|[`jwst_f1800w`](./jwst_f1800w.par)|17.94 μm|2.105 μm||
|$\textrm{F2100W}$|MIRI|[`jwst_f2100w`](./jwst_f2100w.par)|20.70 μm|3.296 μm||
|$\textrm{F2550W}$|MIRI|[`jwst_f2550w`](./jwst_f2550w.par)|25.28 μm|3.461 μm||

The default NIRCam filters are averages over the A and B modules. For every
JWST/NIRCam filter there are also module specific curves available as
`jwst_mod[a,b]_f*w.par`.  For the SW filters these correspond to detector '1' of
each module.


### Roman

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{F062}$|WFI|[`roman_wfi_f062`](./roman_wfi_f062.par)|617.8 nm|184.2 nm|03/24: Average of normalized SCA curves|
|$\textrm{F087}$|WFI|[`roman_wfi_f087`](./roman_wfi_f087.par)|867.1 nm|149.9 nm|03/24: Average of normalized SCA curves|
|$\textrm{F106}$|WFI|[`roman_wfi_f106`](./roman_wfi_f106.par)|1056.3 nm|179.4 nm|03/24: Average of normalized SCA curves|
|$\textrm{F129}$|WFI|[`roman_wfi_f129`](./roman_wfi_f129.par)|1277.2 nm|220.6 nm|03/24: Average of normalized SCA curves|
|$\textrm{F146}$|WFI|[`roman_wfi_f146`](./roman_wfi_f146.par)|1365.3 nm|719.9 nm|03/24: Average of normalized SCA curves|
|$\textrm{F158}$|WFI|[`roman_wfi_f158`](./roman_wfi_f158.par)|1571.6 nm|270.8 nm|03/24: Average of normalized SCA curves|
|$\textrm{F184}$|WFI|[`roman_wfi_f184`](./roman_wfi_f184.par)|1833.7 nm|218.2 nm|03/24: Average of normalized SCA curves|
|$\textrm{F213}$|WFI|[`roman_wfi_f213`](./roman_wfi_f213.par)|2126.4 nm|237.1 nm|03/24: Average of normalized SCA curves|


The default Roman filters are averages over the 18 available SCAs. However, by
asking for `trans_colname="sca??"` you can obtain the transmission curve for the
`??`th SCA (beginning with `sca01`)

### SPHEREx

There are 6 `SPHEREx` 'bands' beginning with `band1`, each having 17 channels per band accessed using `trans_colname='ch??'`, beginning with `ch01`.  So, e.g. ``filt=Filter("spherex_band2", trans_colname="ch14")`` gets you the 14th channel of the 2nd band.

### Spitzer Space Telescope

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{Band}\ 1$|IRAC|[`spitzer_irac_ch1`](./spitzer_irac_ch1.par)|3.538 μm|504.5 nm||
|$\textrm{Band}\ 2$|IRAC|[`spitzer_irac_ch2`](./spitzer_irac_ch2.par)|4.478 μm|667.1 nm||
|$\textrm{Band}\ 3$|IRAC|[`spitzer_irac_ch3`](./spitzer_irac_ch3.par)|5.696 μm|946.9 nm||
|$\textrm{Band}\ 4$|IRAC|[`spitzer_irac_ch4`](./spitzer_irac_ch4.par)|7.798 μm|1.937 μm||
|$24\ \textrm{μm}$|MIPS|[`spitzer_mips_24`](./spitzer_mips_24.par)|23.43 μm|4.497 μm||
|$70\ \textrm{μm}$|MIPS|[`spitzer_mips_70`](./spitzer_mips_70.par)|69.86 μm|19.50 μm||
|$160\ \textrm{μm}$|MIPS|[`spitzer_mips_160`](./spitzer_mips_160.par)|154.3 μm|30.53 μm||

### Wide-field Infrared Survey Explorer (WISE)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{W1}$|WISE|[`wise_w1`](./wise_w1.par)|3.346 μm|647.9 nm||
|$\textrm{W2}$|WISE|[`wise_w2`](./wise_w2.par)|4.595 μm|759.7 nm||
|$\textrm{W3}$|WISE|[`wise_w3`](./wise_w3.par)|11.55 μm|5.851 μm||
|$\textrm{W4}$|WISE|[`wise_w4`](./wise_w4.par)|22.08 μm|3.720 μm||

### SWIFT/UVOT

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$w1$|UVOT|[`uvot_w1`](./uvot_w1.par)|251.7 nm|93.8 nm||
|$w2$|UVOT|[`uvot_w2`](./uvot_w2.par)|201.0 nm|67.2 nm||
|$m2$|UVOT|[`uvot_m2`](./uvot_m2.par)|223.0 nm|45.9 nm||



## Ground-Based Telescopes

### Two Micron All-Sky Survey (2MASS) Telescope

|Filter|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:----------:|---------------------:|---------------------------:|:-------|
|$J$|[`twomass_J`](./twomass_J.par)|1.232 μm|156.4 nm||
|$H$|[`twomass_H`](./twomass_H.par)|1.642 μm|184.7 nm||
|$Ks$|[`twomass_Ks`](./twomass_Ks.par)|2.157 μm|206.9 nm||

### Víctor M. Blanco Telescope

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$u$|DECam|[`decam_u`](./decam_u.par)|380.9 nm|33.07 nm||
|$g$|DECam|[`decam_g`](./decam_g.par)|477.1 nm|99.45 nm||
|$r$|DECam|[`decam_r`](./decam_r.par)|639.2 nm|101.5 nm||
|$i$|DECam|[`decam_i`](./decam_i.par)|778.9 nm|102.0 nm||
|$z$|DECam|[`decam_z`](./decam_z.par)|914.9 nm|101.4 nm||
|$Y$|DECam|[`decam_Y`](./decam_Y.par)|988.5 nm|62.56 nm||

### Canada–France–Hawaii Telescope (CFHT)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$u_s\ \textrm{(9301)}$|MegaCam|[`cfht_megacam_us_9301`](./cfht_megacam_us_9301.par)|381.3 nm|55.39 nm||
|$g_s\ \textrm{(9401)}$|MegaCam|[`cfht_megacam_gs_9401`](./cfht_megacam_gs_9401.par)|483.3 nm|99.01 nm||
|$r_s\ \textrm{(9601)}$|MegaCam|[`cfht_megacam_rs_9601`](./cfht_megacam_rs_9601.par)|622.4 nm|86.20 nm||
|$i_s\ \textrm{(9701)}$|MegaCam|[`cfht_megacam_is_9701`](./cfht_megacam_is_9701.par)|765.1 nm|102.0 nm||
|$z_s\ \textrm{(9801)}$|MegaCam|[`cfht_megacam_zs_9801`](./cfht_megacam_zs_9801.par)|884.7 nm|108.7 nm||
|$J\ \textrm{(8101)}$|WIRCam|[`cfht_wircam_J_8101`](./cfht_wircam_J_8101.par)|1.251 μm|108.0 nm||
|$H\ \textrm{(8201)}$|WIRCam|[`cfht_wircam_H_8201`](./cfht_wircam_H_8201.par)|1.625 μm|196.5 nm||
|$K_s\ \textrm{(8302)}$|WIRCam|[`cfht_wircam_Ks_8302`](./cfht_wircam_Ks_8302.par)|2.154 μm|213.8 nm||

### W. M. Keck Telescope I (Keck I)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$g$|LRIS|[`keck_lris_g`](./keck_lris_g.par)|473.3 nm|66.87 nm||
|$R_s$|LRIS|[`keck_lris_Rs`](./keck_lris_Rs.par)|678.7 nm|103.6 nm||

### Nicholas U. Mayall Telescope

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$U\ \textrm{(k1001)}$|Mosaic|[`mayall_mosaic_U_k1001`](./mayall_mosaic_U_k1001.par)|358.2 nm|46.70 nm||
|$J_1$|NEWFIRM|[`mayall_newfirm_J1`](./mayall_newfirm_J1.par)|1.044 μm|102.8 nm||
|$J_2$|NEWFIRM|[`mayall_newfirm_J2`](./mayall_newfirm_J2.par)|1.193 μm|103.6 nm||
|$J_3$|NEWFIRM|[`mayall_newfirm_J3`](./mayall_newfirm_J3.par)|1.276 μm|98.84 nm||
|$H_1$|NEWFIRM|[`mayall_newfirm_H1`](./mayall_newfirm_H1.par)|1.559 μm|115.9 nm||
|$H_2$|NEWFIRM|[`mayall_newfirm_H2`](./mayall_newfirm_H2.par)|1.705 μm|120.9 nm||
|$K$|NEWFIRM|[`mayall_newfirm_K`](./mayall_newfirm_K.par)|2.164 μm|219.6 nm||

### MPG/ESO Telescope

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$U_{38}\ \textrm{(ESO 841)}$|WFI|[`mpgeso_wfi_U38_eso841`](./mpgeso_wfi_U38_eso841.par)|368.2 nm|29.64 nm||
|$B\ \textrm{(ESO 842)}$|WFI|[`mpgeso_wfi_B_eso842`](./mpgeso_wfi_B_eso842.par)|457.3 nm|68.16 nm||
|$V\ \textrm{(ESO 843)}$|WFI|[`mpgeso_wfi_V_eso843`](./mpgeso_wfi_V_eso843.par)|536.0 nm|61.65 nm||
|$R_c\ \textrm{(ESO 844)}$|WFI|[`mpgeso_wfi_Rc_eso844`](./mpgeso_wfi_Rc_eso844.par)|646.5 nm|110.6 nm||
|$I_c\ \textrm{(ESO 845)}$|WFI|[`mpgeso_wfi_Ic_eso845`](./mpgeso_wfi_Ic_eso845.par)|859.5 nm|139.9 nm||

### Sloan Digital Sky Survey (SDSS) Telescope

|Filter|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:----------:|---------------------:|---------------------------:|:-------|
|$u$|[`sdss_u0`](./sdss_u0.par)|354.6 nm|45.72 nm||
|$g$|[`sdss_g0`](./sdss_g0.par)|467.0 nm|92.79 nm||
|$r$|[`sdss_r0`](./sdss_r0.par)|615.6 nm|81.28 nm||
|$i$|[`sdss_i0`](./sdss_i0.par)|747.2 nm|89.08 nm||
|$z$|[`sdss_z0`](./sdss_z0.par)|891.7 nm|118.3 nm||

### Subaru Telescope

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$J$|MOIRCS|[`subaru_moircs_J`](./subaru_moircs_J.par)|1.250 μm|110.9 nm||
|$H$|MOIRCS|[`subaru_moircs_H`](./subaru_moircs_H.par)|1.631 μm|191.3 nm||
|$K_s$|MOIRCS|[`subaru_moircs_Ks`](./subaru_moircs_Ks.par)|2.154 μm|207.1 nm||
|$B$|Suprime-Cam|[`subaru_suprimecam_B`](./subaru_suprimecam_B.par)|442.7 nm|71.91 nm||
|$V$|Suprime-Cam|[`subaru_suprimecam_V`](./subaru_suprimecam_V.par)|545.5 nm|68.03 nm||
|$R_c$|Suprime-Cam|[`subaru_suprimecam_Rc`](./subaru_suprimecam_Rc.par)|649.0 nm|81.83 nm||
|$r^\prime$|Suprime-Cam|[`subaru_suprimecam_rp`](./subaru_suprimecam_rp.par)|624.9 nm|95.78 nm||
|$i^\prime$|Suprime-Cam|[`subaru_suprimecam_ip`](./subaru_suprimecam_ip.par)|764.6 nm|103.3 nm||
|$z^\prime$|Suprime-Cam|[`subaru_suprimecam_zp`](./subaru_suprimecam_zp.par)|901.1 nm|92.26 nm||
|$\textrm{IA427}$|Suprime-Cam|[`subaru_suprimecam_ia427`](./subaru_suprimecam_ia427.par)|425.9 nm|15.13 nm||
|$\textrm{IA445}$|Suprime-Cam|[`subaru_suprimecam_ia445`](./subaru_suprimecam_ia445.par)|444.2 nm|15.00 nm||
|$\textrm{IA464}$|Suprime-Cam|[`subaru_suprimecam_ia464`](./subaru_suprimecam_ia464.par)|463.2 nm|15.81 nm||
|$\textrm{IA484}$|Suprime-Cam|[`subaru_suprimecam_ia484`](./subaru_suprimecam_ia484.par)|484.6 nm|16.83 nm||
|$\textrm{IA505}$|Suprime-Cam|[`subaru_suprimecam_ia505`](./subaru_suprimecam_ia505.par)|506.0 nm|17.22 nm||
|$\textrm{IA527}$|Suprime-Cam|[`subaru_suprimecam_ia527`](./subaru_suprimecam_ia527.par)|525.8 nm|18.49 nm||
|$\textrm{IA550}$|Suprime-Cam|[`subaru_suprimecam_ia550`](./subaru_suprimecam_ia550.par)|549.4 nm|20.35 nm||
|$\textrm{IA574}$|Suprime-Cam|[`subaru_suprimecam_ia574`](./subaru_suprimecam_ia574.par)|576.2 nm|20.49 nm||
|$\textrm{IA598}$|Suprime-Cam|[`subaru_suprimecam_ia598`](./subaru_suprimecam_ia598.par)|600.6 nm|22.20 nm||
|$\textrm{IA624}$|Suprime-Cam|[`subaru_suprimecam_ia624`](./subaru_suprimecam_ia624.par)|622.9 nm|22.55 nm||
|$\textrm{IA651}$|Suprime-Cam|[`subaru_suprimecam_ia651`](./subaru_suprimecam_ia651.par)|649.7 nm|24.03 nm||
|$\textrm{IA679}$|Suprime-Cam|[`subaru_suprimecam_ia679`](./subaru_suprimecam_ia679.par)|678.0 nm|25.04 nm||
|$\textrm{IA709}$|Suprime-Cam|[`subaru_suprimecam_ia709`](./subaru_suprimecam_ia709.par)|707.2 nm|23.69 nm||
|$\textrm{IA738}$|Suprime-Cam|[`subaru_suprimecam_ia738`](./subaru_suprimecam_ia738.par)|735.8 nm|23.74 nm||
|$\textrm{IA767}$|Suprime-Cam|[`subaru_suprimecam_ia767`](./subaru_suprimecam_ia767.par)|767.9 nm|26.27 nm||
|$\textrm{IA797}$|Suprime-Cam|[`subaru_suprimecam_ia797`](./subaru_suprimecam_ia797.par)|796.4 nm|26.61 nm||
|$\textrm{IA827}$|Suprime-Cam|[`subaru_suprimecam_ia827`](./subaru_suprimecam_ia827.par)|824.5 nm|24.62 nm||
|$\textrm{IA856}$|Suprime-Cam|[`subaru_suprimecam_ia856`](./subaru_suprimecam_ia856.par)|856.3 nm|26.78 nm||
|$g$|Hyper Suprime-Cam|[`hsc_g`](./hsc_g.par)|475.5 nm|97.7 nm||
|$r$|Hyper Suprime-Cam|[`hsc_r`](./hsc_r.par)|618.4 nm|101.4 nm||
|$i$|Hyper Suprime-Cam|[`hsc_i`](./hsc_i.par)|766.1 nm|107.9 nm||
|$z$|Hyper Suprime-Cam|[`hsc_z`](./hsc_z.par)|889.7 nm|55.2 nm||
|$y$|Hyper Suprime-Cam|[`hsc_y`](./hsc_y.par)|976.2 nm|73.1 nm||

### United Kingdom Infrared Telescope (UKIRT)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$J$|WFCAM|[`ukirt_wfcam_J`](./ukirt_wfcam_J.par)|1.248 μm|112.2 nm||
|$H$|WFCAM|[`ukirt_wfcam_H`](./ukirt_wfcam_H.par)|1.631 μm|204.6 nm||
|$K$|WFCAM|[`ukirt_wfcam_K`](./ukirt_wfcam_K.par)|2.201 μm|246.5 nm||

### Visible and Infrared Survey Telescope for Astronomy (VISTA)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$\textrm{Z}$|VIRCAM|[`vista_vircam_Z`](./vista_vircam_Z.par)|878.3 nm|86.8 nm|Filter+QE+atm|
|$\textrm{Y}$|VIRCAM|[`vista_vircam_Y`](./vista_vircam_Y.par)|1.020 μm|78.3 nm|Filter+QE+atm|
|$\textrm{J}$|VIRCAM|[`vista_vircam_J`](./vista_vircam_J.par)|1.251 μm|122.7 nm|Filter+QE+atm|
|$\textrm{H}$|VIRCAM|[`vista_vircam_H`](./vista_vircam_H.par)|1.639 μm|201.9 nm|Filter+QE+atm|
|$\textrm{KS}$|VIRCAM|[`vista_vircam_Ks`](./vista_vircam_Ks.par)|2.145 μm|261.5 nm|Filter+QE+atm|
|$\textrm{Z}$|VIRCAM|[`vista_vircam_Z_FILTER`](./vista_vircam_Z_FILTER.par)|831.5 nm|734.6 nm|Filter only|
|$\textrm{Y}$|VIRCAM|[`vista_vircam_Y_FILTER`](./vista_vircam_Y_FILTER.par)|1027.2 nm|223.9 nm|Filter only|
|$\textrm{J}$|VIRCAM|[`vista_vircam_J_FILTER`](./vista_vircam_J_FILTER.par)|1256.8 nm|229.3 nm|Filter only|
|$\textrm{H}$|VIRCAM|[`vista_vircam_H_FILTER`](./vista_vircam_H_FILTER.par)|1639.7 nm|207.1 nm|Filter only|
|$\textrm{KS}$|VIRCAM|[`vista_vircam_Ks_FILTER`](./vista_vircam_Ks_FILTER.par)|2139.6 nm|267.4 nm|Filter only|

### Very Large Telescope (VLT) Array

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$J$|ISAAC|[`vlt_isaac_J`](./vlt_isaac_J.par)|1.231 μm|178.0 nm||
|$H$|ISAAC|[`vlt_isaac_H`](./vlt_isaac_H.par)|1.645 μm|199.2 nm||
|$K_s$|ISAAC|[`vlt_isaac_Ks`](./vlt_isaac_Ks.par)|2.164 μm|185.5 nm||
|$U$|VIMOS|[`vlt_vimos_U`](./vlt_vimos_U.par)|374.2 nm|39.07 nm||
|$R$|VIMOS|[`vlt_vimos_R`](./vlt_vimos_R.par)|641.7 nm|94.24 nm||

### BASS and MzLS (https://legacysurvey.org)

|Filter|Instrument|`sedpy` Name|$\lambda_\textrm{eff}$|$\Delta\lambda_\textrm{eff}$|Comments|
|:----:|:--------:|:----------:|---------------------:|---------------------------:|:-------|
|$g$|90-Prime|[`bass_g`](./bass_g.par)|472.1 nm|100.2 nm||
|$r$|90-Prime|[`bass_r`](./bass_r.par)|637.0 nm|97.76 nm||
|$z$|Mosaic-3|[`mzls_z`](./mzls_z.par)|917.9 nm|100.7 nm||
