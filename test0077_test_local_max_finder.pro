
;main level program
compile_opt idl2
ireset, /no_prompt

;make some test data
xvals = findgen(2000)
yvals = 2*sin( xvals*(2*!PI/500) ) + sin( xvals*(2*!PI/1000) )
;plot the test data
p = plot(xvals, yvals,$
    xtitle = 'X position', ytitle = 'Amplitude',$
    yminor = 3, xminor = 3, title = 'Neat Local Max/Min example')
    

;find the local max index
local_maxmin_index = local_max_finder(xvals, yvals)

;extract the x/y pairings of the local max/min
x_extrema = xvals[local_maxmin_index]
y_extrema = yvals[local_maxmin_index]

;overplot the max/min on the existing plot
p2 = scatterplot(x_extrema, y_extrema, /current, /overplot, $
    symbol = 'o', sym_color = 'r', sym_thick = 2)


;find the local min index with /minima keyword
local_maxmin_index = local_max_finder(xvals, yvals, /minima)

;extract the x/y pairings of the local max/min
x_extrema = xvals[local_maxmin_index]
y_extrema = yvals[local_maxmin_index]

;overplot the min on the existing plot
p3 = scatterplot(x_extrema, y_extrema, /current, /overplot, $
    symbol = 'o', sym_color = 'b', sym_thick = 2)

end
