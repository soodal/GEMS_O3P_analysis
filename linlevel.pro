function linlevel, start, final, nval
return, findgen(nval, increment=(float(final) - float(start) + 1)/float(nval)) + float(start)
end
