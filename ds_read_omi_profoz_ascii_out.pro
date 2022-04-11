function ds_read_omi_profoz_ascii_out, fnames

silence = 0

read_omil2_file, fnames, nl, nf, nalb, ngas, naer,nprof, $
    omilon, omilat, omisza, omivza, omiaza, omirms, omiavgres, $
    omicfrac, omictp, omicldflg, omiai, omiutc,omintp,   $ 
    ominw, omisaa, omiexval, ominiter, ominspike, omiglint, omildwt, omisnow, $
    omimon, omiday, omiyear, omipix, omiline, omitaod, omitsca, omisaod, omialb, atmos, $
    ozprofs, omicol, omitrace,  omifitvar, omiavgk, omicorrel, omicovar, omicontri, omifitspec, $
    omisimrad, omiwaves, omiclmrad, omiactrad, omiwf, omisnr, omiring, $
    nfail, flon, flat, fsza, fvza, faza, fmon, fday, fyear, $
    fpix,  fline, fsaa, fexval, fniter, fnspike, fglint, fldwt, fsnow, $
    orbits=orbits, orbspix = orbspix, orbepix = orbepix, omiorb=omiorb,varname = varname, $
    get_trace=get_trace, get_fitvar=get_fitvar, get_correl=get_correl, get_covar=get_covar, $
    get_avgk = get_avgk, get_contri = get_contri, get_fitres=get_fitres, $
    get_radcal = get_radcal, get_wf = get_wf, get_snr = get_snr, get_ring=get_ring,  $
    omitraceprof=omitraceprof, omitraceavgk=omitraceavgk, omitracecontri=omitracecontri, omialb0=omialb0, $
    omiptrp = omiptrp, omiztrp=omiztrp, silence=silence



result = {nl:nl, nf:nf, nalb:nalb, ngas:ngas, naer:naer, nprof:nprof, $
  omilon:omilon, omilat:omilat, omisza:omisza, omiaza:omiaza, omirms:omirms, $
  omiavgres:omiavgres, omicfrac:omicfrac, omictp:omictp, omicldflg:omicldflg, $
  omiai:omiai, omiutc:omiutc, omintp:omintp, ominw:ominw, omisaa:omisaa, $
  omiexval:omiexval, ominiter:ominiter, ominspike:ominspike, omitsca:omitsca, $
  omisaod:omisaod, omialb:omialb, atmos:atmos, ozprofs:ozprofs, $
  omicol:omicol, omitrace:omitrace, omifitvar:omifitvar, omiavgk:omiavgk, $
  omicorrel:omicorrel, omicovar:omicovar, omicontri:omicontri, $
  omifitspec:omifitspec, omisimrad:omisimrad, omiwaves:omiwaves, $
  omiclmrad:omiclmrad, omiactrad:omiactrad, omiwf:omiwf, omisnr:omisnr, $
  omiring:omiring, nfail:nfail, flon:flon, flat:flat, fsza:fsza, fvza:fvza, $
  faza:faza, fmon:fmon, fday:fday, fyear:fyear, fpix:fpix, fline:fline, $
  fsaa:fsaa, fexval:fexval, fniter:fniter, fnspike:fnspike, fglint:fglint, $
  fldwt:fldwt, fsnow:fsnow, orbits:orbits, orbspix:orbspix, orbepix:orbepix, $
  omiorb:omiorb, varname:varname, $
  ;get_trace:get_trace, $
  ;get_fitvar:get_fitvar, $
  ;get_correl:get_correl, $
  ;get_covar:get_covar, $
  ;get_avgk:get_avgk, $
  ;get_contri:get_contri, $
  ;get_fitres:get_fitres, $
  ;get_radcal:get_radcal, $
  ;get_wf:get_wf, $
  ;get_snr:get_snr, $
  ;get_ring:get_ring, $
  omitraceprof:omitraceprof, omitraceavgk:omitraceavgk, $
  omitracecontri:omitracecontri, omialb0:omialb0, omiptrp:omiptrp, $
  omiztrp:omiztrp, silence:silence}



return, result
end

