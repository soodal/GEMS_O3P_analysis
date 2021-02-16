pro ds_xymake2d, x, y, x2d, y2d

xsize = size(x)
ysize = size(y)

if xsize[0] ne 1 then begin
  print, 'x is not 1 dimensional variable'
  return
endif

if ysize[0] ne 1 then begin
  print, 'y is not 1 dimensional variable'
  return

endif
x2d = rebin(x, xsize[1], ysize[1])
y2d = transpose(rebin(y, ysize[1], xsize[1]))

return
end
