pro ds_plot_geographic, data, xx, yy, limit=limit, range=range, $
  pngfile=pngfile, title=title, cbtitle=cbtitle

if not keyword_set(range) then begin
  range=minmax(data)
endif
if not keyword_set(title) then begin
  title=''
endif
if not keyword_set(cbtitle) then begin
  cbtitle=''
endif
m = map('Geographic', /buffer, limit=limit)
c = contour(reform(data), xx, yy, /fill, /buffer, $
  n_levels=100, $
  rgb_table=colortable(70, /reverse), $
  ;zrange=[240, 320], $
  max_value=range[1], $
  min_value=range[0], $
  title=title, over=m)
map = mapcontinents()
cb = colorbar( $
  orientation=0, $
  position=[0.1, 0.10, 0.9, 0.15], $
  tickvalues=tickvalues,$
  ;range=[300, 320], $
  ;color='black', $
  title=cbtitle, $
  taper=1)
if not keyword_set(pngfile) then begin 
  pngfile = './plot/test.png'
endif

;c.title.font_size=16
c.save, pngfile
c.close
; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'

end
