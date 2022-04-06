pro ds_merra2_accum, data, o3accum, itime=itime, height=height, hpa=hpa
if (not keyword_set(height)) and (not keyword_set(hpa)) then BEGIN
  message, 'height or hpa keyword need to be set'
endif
if keyword_set(height) and keyword_set(hpa) then BEGIN
  message, 'height and hpa keyword cannot be set simultaneously'
endif

if not keyword_set(itime) then BEGIN
  print, 'keyword itime is not set. itime is set to 1'
  itime = 1
endif

EPV = data.EPV
delp = data.delp
delp = reform(delp[*, *, *, itime])
H = data.H
H = reform(h[*, *, *, itime])
lat = data.lat
lev = data.lev
lon = data.lon
O3 = data.o3
o3 = reform(o3[*, *, *, itime])

OMEGA = data.OMEGA
PHIS = data.PHIS
PS = data.PS
PL = data.PL
pl = reform(pl[*, *, *, itime])
QI = data.QI
QL = data.QL
RH = data.RH
SLP = data.SLP
T = data.T
time = data.time
U = data.U
v = data.V

lat2d = data.lat2d
lon2d = data.lat2d

sz = size(o3)
dim1 = sz[1]
dim2 = sz[2]
dim3 = sz[3]

pres = 0.0

m2mm = 1000.
mm2du = 100.

;convert mass mxing ratio into vmr
vmr_O3 = 28.9644/47.99 * reform(o3[*,*,*])

vmr_o3_ds = ds_mmr2vmr(o3[*, *, *], /ppmv)

;vmr to No3 (number density , molec/cm3)
k = 1.3806505 * 10.^(-23); ![kg m^2  s^-2 K^-1] = [J K^-1] = [N m K^-1]

;O3 number density
NO3 = vmr_o3 * pl/(T*k)

;Air number density
Nair = 2.69*10.^(25) ; molec/m ^2

dh      = fltarr(sz[1],sz[2],sz[3])

o3du = fltarr(sz[1], sz[2], sz[3])
dh[*, *, 0:70] = h[*, *,0:70] - h[*, *, 1:71] 
dh[*, *, 71] = h[*, *, 70] - h[*, *, 71]

o3du = double(no3/nair*dh * m2mm * mm2du) ;m -> cm thickness

;o3du[where(o3du GE 1000 or o3du lt 0)] = !values.f_nan

;tco = total(o3du[*, *, 0:sz[3]-1],/nan,3)


trco = fltarr(dim1, dim2)


ds_xymake2d, lon, lat, lon2d, lat2d

;;----------------------------
;; MERRA2 vertical column layer for under 300 hPa
;;----------------------------
if keyword_set(hpa) then begin
  ;;----------------------------
  ;; construct MERRA2 boundary pressure for layers
  ;;----------------------------
  pboundary = fltarr(dim1, dim2, dim3+1)
  pboundary[*, *, 0:71] = pl - delp/2.
  pboundary[*, *, 72] = pl[*, * ,71] + delp[*, *, 71]/2.
  pboundary = pboundary / 100.
  for iy=0, dim2-1 do begin
  for ix=0, dim1-1 do begin
    for ilevel=dim3, 0, -1 do begin
      ilayer = ilevel-1
      ; if upper and lower boundary are both lower than the specific level
      if pboundary[ix, iy, ilevel] gt hpa $ ; lower boundary
          and pboundary[ix, iy,ilevel-1] gt hpa then begin ; upper boundar 
        trco[ix, iy] = trco[ix, iy] + o3du[ix, iy, ilayer]

      ; if upper boundary is higher and lower boundary is lohwer than the specific level
      endif else if pboundary[ix, iy, ilevel] gt hpa $
          and pboundary[ix, iy, ilevel-1] le hpa then begin

        trco[ix, iy] = trco[ix, iy] + o3du[ix, iy, ilayer]*$
          (pboundary[ix, iy, ilevel]-hpa)/$
          (pboundary[ix, iy, ilevel]-pboundary[ix, iy, ilevel-1])
        break
      endif
    endfor
  endfor
  endfor
  trco[where(trco lt 0)] = !values.f_nan
  o3accum = trco
  return
endif

; GEMS vertical column layer for under 10 km
if keyword_set(height) then begin
  ;;----------------------------
  ;; construct MERRA2 boundary pressure for layers
  ;;----------------------------
  hboundary = fltarr(dim1, dim2, dim3+1)
  hboundary[*, *, 0:70] = h[*, *, 0:70] + (h[*, *, 0:70] - h[*, *, 1:71])/2.
  hboundary[*, *, 0:72] = h[*, * ,71] - (h[*, *, 70] - h[*, *, 71])/2.
  hboundary = hboundary / 1000.

  for iy=0, dim2-1 do begin
  for ix=0, dim1-1 do begin
    for ilevel=dim3, 0, -1 do begin
      ilayer = ilevel-1
      if hboundary[ix, iy, ilevel] lt height $
          and hboundary[ix, iy, ilevel-1] lt height then begin
        trco[ix, iy] = trco[ix, iy] + o3[ix, iy, ilayer]
        ;print, ilayer

      endif else if hboundary[ix, iy, ilevel] lt height $
          and hboundary[ix, iy, ilevel-1] ge height then begin

        trco[ix, iy] = trco[ix, iy] + o3[ix, iy, ilayer]*$
          (height-hboundary[ix, iy, ilevel])/$
          (hboundary[ix, iy, ilevel-1]-hboundary[ix, iy, ilevel])
        ;print, ilayer
        break
      endif
    endfor
  endfor
  endfor
  trco[where(trco le 0)] = !values.f_nan
  o3accum = trco
  return
endif
return
end
