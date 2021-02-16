
;function plot_gems_o3p_profile, o3, p, color=color, overplot=overplot
;;if not keyword_set(color) then begin
  ;;color = 'black'
;;endif
;;if not keyword_set(overplot) then begin
  ;;overplot=1
;;endif

;p1 = plot(o3, p, /buffer, yrange=[1000, 0], $
  ;color=color, $
  ;;overplot=overplot)
  ;overplot=1)
;return, p1
;end

pos = [0.12, 0.1, 0.9, 0.85]
leg_pos = [0.85, 0.8]

path = '../GEMS_O3P_Yonsei/out/measurement_error_test/'
ci = [30, 50, 70, 90];102
cs = [100, 200, 300];378

for ici = 0, n_elements(ci)-1 do begin
  for icj = 0, n_elements(cs)-1 do begin
    cics = 'ci' + strtrim(string(ci[ici], format='(i3)'), 2) + 'cs' + string(cs[icj], format='(i3)')

    color = ['red', 'blue', 'green', 'cyan'] 
    name = ['ME 0.01', '0.001', '0.005', '0.04']
    o3p_total = fltarr(4, 24)

    for ime = 0, 3 do begin
      fp = path + cics + '/' + string(ime, format='(i03)') + '/'  
      o3fn = 'output_l2o3p_ozret.txt'
      presfn = 'output_l2o3p_p.txt'
      extrafn = 'output_l2o3p_extra_info.txt'
      ozapfn = 'output_l2o3p_ozap.txt'

      readcol, fp + o3fn, o3, format='F', numline=24
      readcol, fp + presfn, p, format='F', numline=24
      readcol, fp + extrafn, info, numline=4
      readcol, fp + ozapfn, ozap, format='F', numline=24

      o3p_total[ime,*] = o3
      if ime eq 0 then begin
        command = "p" + string(ime, format='(i1)') + $
          " = plot(o3, p, color=color[ime], overplot=ime, /buffer, " + $
          " title='GEMS O3PROFILE 2020-06-16T0345', " + $
          " xtitle='O3[DU]', " + $
          " ytitle='Pressure[hPa]', " + $
          " yrange=[1000, 0], /ylog, pos=pos, name=name[ime])"
        print, command
      endif else begin
        command = "p" + string(ime, format='(i1)') + $
          " = plot(o3, p, color=color[ime], overplot=ime, /buffer, " + $
          " yrange=[1000, 0], /ylog, pos=pos, name=name[ime])"
      endelse
      dumm = execute(command)

    endfor

    pap = plot(o3, p, color='black', /overplot, /buffer, $
      yrange=[1000, 0], /ylog, pos=pos, name='A priori')
      
    t0 = text(0.45, 0.3, 'Final Quality Flag:' + string(info[0], format='(I1)'))
    t1 = text(0.45, 0.25, 'Pressure Quality Flag:' + string(info[1], format='(I1)'))
    t2 = text(0.45, 0.20, 'Cloud Pressure:' + string(info[2]))
    t3 = text(0.45, 0.15, 'Effective Cloud Pressure:' + string(info[3], format='(F7.2)'))

    leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
    pngfile = './plot/gems_l2_o3_profile_' + cics + '.png'
    p1.save, pngfile
    p1.close
    ; send image to pc
    spawn, 'scp -P18742 -p ' + pngfile + $
      ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp/profiles'
    print, o3p_total
  endfor
endfor
end




