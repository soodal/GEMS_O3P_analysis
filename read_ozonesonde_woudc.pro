pro read_ozonesonde_woudc, filename, result_struct

openr, lun, filename, /get_lun
s = ''
while ~ eof(lun) do begin
  readf, lun, s
  ;print, s
  if strmid(s, 0, 1) eq '*' then begin
    print, 'header detected.'
    slen = strlen(s)
    line = strsplit(strmid(s, 2,slen-2),': ', /extract, /regex, /preserve_null)
    if line[0] eq 'SOFTWARE' then begin
      SOFTWARE = line[1]
    endif else if line[0] eq 'Ground equipment' then begin
      Ground_equipment = line[1]
    endif else if line[0] eq 'Sampling interval' then begin
      Sampling_interval = line[1]
    endif else if line[0] eq 'Smoothing interval' then begin
      Smoothing_interval = line[1]
    endif else if line[0] eq 'Amount of cathode solution' then begin
      Amount_of_cathode_solution = line[1]
    endif else if line[0] eq 'Concentration of cathode solution' then begin
      Concentration_of_cathod_solution = line[1]
    endif else if line[0] eq 'Background current correction' then begin
      Background_current_correction = line[1]
    endif else if line[0] eq 'Pump efficiency' then begin
      Pumb_efficiency = line[1]
    endif else if line[0] eq 'Balloon weight' then begin
      Balloon_weight = line[1]
    endif else if line[0] eq 'Tropopause' then begin
      Tropopause = line[1]
    endif
  endif else if strmid(s, 0, 1) eq '#' then begin
    slen = strlen(s)
    title = strmid(s, 1,slen-1)
    if title eq 'CONTENT' then begin
      read_or_not = 1
      while read_or_not do begin 
        readf, lun, s
        ;print, s
        if strlen(s) eq 0 then begin
          read_or_not = 0
        endif
      endwhile
    endif else if title eq 'DATA_GENERATION' then begin
      read_or_not = 1
      while read_or_not do begin
        readf, lun, s
        ;print, s
        if strlen(s) eq 0 then begin
          read_or_not = 0
        endif
      endwhile
    endif else if title eq 'PLATFORM' then begin
      read_or_not = 1
      while read_or_not do begin
        readf, lun, s
        ;print, s
        if strlen(s) eq 0 then begin
          read_or_not = 0
        endif
      endwhile
    endif else if title eq 'INSTRUMENT' then begin
      read_or_not = 1
      while read_or_not do begin
        readf, lun, s
        ;print, s
        if strlen(s) eq 0 then begin
          read_or_not = 0
        endif
      endwhile
    endif else if title eq 'LOCATION' then begin
      ;read_or_not = 1
      ;while read_or_not do begin 
        ;readf, lun, s
        ;print, s
        ;if strlen(s) eq 0 then begin
          ;read_or_not = 0
        ;endif
        ;strsplit(s, ',')
      ;endwhile
      readf, lun, s
      colname = strsplit(s, ',', /regex, /extract)
      readf, lun, s
      column = strsplit(s, ',', /regex, /extract)
      launch_latitude = float(column[0])
      launch_longitude = float(column[1])
      launch_height = float(column[2])
    endif else if title eq 'TIMESTAMP' then begin
      read_or_not = 1
      while read_or_not do begin 
        readf, lun, s
        ;print, s
        if strlen(s) eq 0 then begin
          read_or_not = 0
        endif
      endwhile
    endif else if title eq 'FLIGHT_SUMMARY' then begin
      read_or_not = 1
      while read_or_not do begin 
        readf, lun, s
        ;print, s
        if strlen(s) eq 0 then begin
          read_or_not = 0
        endif
      endwhile
    endif else if title eq 'AUXILIARY_DATA' then begin
      read_or_not = 1
      while read_or_not do begin 
        readf, lun, s
        ;print, s
        if strlen(s) eq 0 then begin
          read_or_not = 0
        endif
      endwhile
    endif else if title eq 'PUMP_CORRECTION' then begin
      read_or_not = 1
      while read_or_not do begin 
        readf, lun, s
        ;print, s
        if strlen(s) eq 0 then begin
          read_or_not = 0
        endif
      endwhile
    endif else if title eq 'PROFILE' then begin
      readf, lun, s
      read_or_not = 1
      column_name_profile = strsplit(s, ',', /extract, /preserve_null)
      ncol_profile = n_elements(column_name_profile)
      for i = 0, ncol_profile-1 do begin
        dummy = column_name_profile[i] + ' = []'
        dummy = execute(dummy)
      endfor
        
      while (not eof(lun)) and read_or_not do begin 
        readf, lun, s
        if strlen(s) eq 0 then begin
          read_or_not = 0
          break
        endif

        row = strsplit(s, ',', /extract, /preserve_null)
        nullidx = where(strlen(row) eq 0, /null)
        ;if n_elements(nullidx) ne 1 then begin
          ;print, nullidx, s
        ;endif
        row[nullidx] = 'nan'
        for i = 0, n_elements(row)-1 do begin
          dummy = column_name_profile[i] + ' = [' $
            + column_name_profile[i] + ', float(row[' + string(i, format='(i2)') + '])]'
          dummy = execute(dummy)
        endfor
      endwhile
    endif else if title eq 'PRELAUNCH' then begin
      ;print, 'prelaunch start'
      ;print, n_elements(pressure), n_elements(WindSpeed)
      readf, lun, s
      read_or_not = 1
      column_name_prelaunch = strsplit(s, ',', /extract, /preserve_null)
      column_name_prelaunch = column_name_prelaunch + '_prelaunch'
      ncol_prelaunch = n_elements(column_name_prelaunch)
      for i = 0, ncol_profile-1 do begin
        dummy = column_name_prelaunch[i] + ' = []'
        dummy = execute(dummy)
      endfor

      while read_or_not do begin 
        readf, lun, s
        if strlen(s) eq 0 then begin
          read_or_not = 0
          break
        endif

        column_prelaunch = strsplit(s, ',', /extract, /preserve_null)
        nullidx = where(strlen(column_prelaunch) eq 0, /null)
        column_prelaunch[nullidx] = 'nan'
        for i = 0, n_elements(column_prelaunch)-1 do begin
          dummy = column_name_prelaunch[i] + ' = [' $
            + column_name_prelaunch[i] + ', float(column_prelaunch[' + string(i, format='(i2)') + '])]'
          dummy = execute(dummy)
        endfor
      endwhile
    endif
  endif
endwhile

; add variable to struct

exist = size(column_name_profile, /dimension)
if exist ge 1  then begin
  for i = 0, ncol_profile-1 do begin
    struct_exist = size(result_struct, /dimension)
    if struct_exist then begin
      dummy = 'result_struct = create_struct(result_struct, "' + $
        column_name_profile[i] + '", ' + column_name_profile[i] + ')'
      endif else begin
      dummy = 'result_struct = create_struct("' + column_name_profile[i] + $
        '", ' + column_name_profile[i] + ')'
    endelse
    dummy_result = execute(dummy)
  endfor
endif

exist = size(column_name_prelaunch, /dimension)
if exist ge 1 then begin
  for i = 0, ncol_prelaunch-1 do begin
    struct_exist = size(result_struct, /dimension)
    if struct_exist then begin
      dummy = 'result_struct = create_struct(result_struct, "' + $
        column_name_prelaunch[i] + '_prelaunch", ' + column_name_prelaunch[i] + ')'
      endif else begin
      dummy = 'result_struct = create_struct("' + column_name_prelaunch[i] + $
        '_prelaunch", ' + column_name_prelaunch[i] + ')'
    endelse
    dummy_result = execute(dummy)
  endfor
endif

result_struct =create_struct(result_struct, 'launch_latitude', launch_latitude)
result_struct =create_struct(result_struct, 'launch_longitude', launch_longitude)
result_struct =create_struct(result_struct, 'launch_height', launch_height)

free_lun, lun
return

end
