; Script to simulate temperature image cubes of Didymos/Dimorphos based on the
; pre-calculated T_ref from T.N. Titus
;
; Models are pre-processed to break down to individual files each containing a
; temperature model in the format fully complying with the requirement of
; 'tpm_sph2plt.pro'.
;
; Input temperatur emodel are saved in subdirectory 'KRC/'.  Each individual
; model file has a name like 'tref_gamma{:d}_emis{:.1f}.fits'.format(ti, emis).


; switch for Juno shape model or triaxial ellipsoid
shape = 'gany'   ; use Juno shape model

if shape eq 'gany' then begin
    shapefile = '/Users/zouxd/Documents/2-Projects/2024-YORPD/Ganymed/Ganymed_261shape_tmp.plt';
    trefdir = '/Users/zouxd/Documents/2-Projects/2024-YORPD/Ganymed/KRC/Ganymed/'
    outdir = '/Users/zouxd/Documents/2-Projects/2024-YORPD/Ganymed/temp_maps2/'
    subelat = 35.9
    subslat = 18.8
    delta_lon = 27.1
    shape_a = 0.4250
    shape_b = 0.4250
    shape_c = 0.310
    pxlscl = 2.8 ; pixel scale in mas

endif 

; longitude step size for full lightcurve calculation
dlon = 90      ; longitudinal step size for rotation

; latitudes
rh = 1.31 ;solar distance
delta = 0.39 ; range from Earth

; simulated image parameters

xs = 256   ; x size of image
ys = 256   ; y size of image


;-------------------------------------------

; load observing geometry for ALMA data
;readcol, geofile, utc, rh, delta, phase, ra, dec, subelat, subelon, subslat, subslon, polepa, poleinc, sunpa, suninc, delim=',', skipline=1

; load shape model
 readshape_triplate, shapefile, vert, tri
 print, 'Dimensions of vert, tri: ', size(vert), size(tri)
 
; use ellipsoidal shape
;mesh_ellipsoid, 72, 36, shape_a, shape_b, shape_c, vert, tri


;Tfiles = file_search(trefdir + '/*krc.fits')  ; Ceres temperature model for test
Tfiles = file_search(trefdir + 'tref_gamma*_emis*.fits')  ; Juno temperature models

; fits keys to propagrate to output
keys = ['ti', 'emiss', 'rho', 'c', 'p_orb', 'p_rot', 'a_skin', 'd_skin']


; prepare longitudes
nlon = 360 / dlon
subelon = findgen(nlon) / nlon * 360
subelat = replicate(subelat, nlon)
subslon = (subelon + delta_lon) mod 360
subslat = replicate(subslat, nlon)
rh = replicate(rh, nlon)
delta = replicate(delta, nlon)


for j=0, n_elements(Tfiles) - 1 do begin

    print, 'processing file: '+(strsplit(Tfiles[j], '/', /ext))[-1]

    ; load reference temperature array
    tref = readfits(Tfiles[j], hdr, /silent)
    lst = readfits(Tfiles[j], ext=1, /silent)
    hour_angle = ((lst - 12) * 15 + 360) mod 360  ; hour angle from local noon in degrees
    zz = readfits(Tfiles[j], ext=2, /silent)
    lats = readfits(Tfiles[j], ext=3, /silent)

    ; prepare output directory
    dirname = outdir + strjoin((strsplit(file_basename(Tfiles[j]), '.', /ext))[0:-2], '.')
    if not file_test(dirname, /dir) then file_mkdir, dirname

    ; loop throughout longitudes
    for i=0, n_elements(subslon)-1 do begin ; nlon-1 do begin

        print, 'Sub-Earth longitude ', subelon[i]
        ; print, subslon[i], subslat[i], subelon[i], subelat[i], rh[i], delta[i]

        ; calculate observing geometry
        sunpos = vect_match_view(rd2xyz([subslon[i], subslat[i]]), subelat[i], subelon[i], 0) * rh[i] * 1.496e8
        vert1 = vect_match_view(vert, subelat[i], subelon[i], 0)
        res = pxlscl / 206265000 * delta[i] * 1.496e8   ; image resolution in km/s
        ;print, 'resolution: ', res

 ;      print, 'run mesh_geomap'
;        print, 'vert1 =', vert1, 'tri = ', tri
        mesh_geomap, vert1, tri, sunpos, delta[i] * 1.496e8, imap, emap, amap, mask, pltmap, xres=res, yres=res, xs=ys, ys=ys
  ;      print, 'mesh_geomap completed'

        ; save plate map file
        suffix = '_' + string(round(subelon[i]), format='(i3.3)') + '.fits'
        
        ; save plate files
        outfile = dirname + '/platemap' + suffix
        mkhdr, hdr1, pltmap, /ext
        fxaddpar, hdr1, 'long', subelon[i], 'longitude (deg)'
        writefits, outfile, pltmap, hdr1
        writefits, outfile, emap, /append
        writefits, outfile, imap, /append
        writefits, outfile, amap, /append

        ; calculate tpm for plate shape model
        temp_list = tpm_sph2plt(vert, tri, tref, subslon[i], time=[hour_angle, 360], lats=lats)
        ;writefits,outdir+'templist_'+string(i,format='(i3.3)')+'.fits',temp_list

        ; temperature image
        tempmap = tpm_mapping(temp_list, pltmap, depth=zz, zz=zz)
        
        ; save simulation images
        outfile = dirname + '/tempmap' + suffix
        ; prepare primary extension
        mkhdr, hdr1, tempmap, /ext
        fxaddpar, hdr1, 'bunit', 'K'
        for k=0, n_elements(keys) - 1 do begin
            val = sxpar(hdr, keys[k], comment=com)
            fxaddpar, hdr1, keys[k], val, com
        endfor
        fxaddpar, hdr1, 'long', subelon[i], 'longitude (deg)'
        writefits, outfile, tempmap, hdr1
        ; first extension: depth
        mkhdr, hdr1, zz, /ima
        fxaddpar, hdr1, 'bunit', 'm'
        writefits, outfile, zz, hdr1, /append
        ; second extension: emission angle
        mkhdr, hdr1, emap, /ima
        fxaddpar, hdr1, 'bunit', 'deg'
        writefits, outfile, emap, hdr1, /append

    endfor

endfor

end
