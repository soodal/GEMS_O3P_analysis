pro sphergeom_intermediate, lat1, lon1, lat2, lon2, c0, c, lat, lon

    ;real (kind=r8),    intent (in) :: lat1, lat2, lon1, lon2, c0, c
    ;real (kind=r8),  intent (out)  :: lat, lon
    pi         = 3.14159265358979d
    pihalf     = 0.5 * pi
    twopi      = 2.0 * pi

    ;! ---------------
    ;! Local variables
    ;! ---------------
    ;real (kind=r8)  :: x, y, z, tmp1, tmp2, frc, gamsign, theta, gam0, gam

    lat = 0.0
    lon = 0.0
    gam0 = angle_minus_twopi(lon2 - lon1, pi) ;! lon difference
    gamsign = abs(gam0) / gam0
    gam0 = abs(gam0)
    if (gam0 eq 0) then gamsign = 1

    ;! Get straight line (AB) segment fraction frc intercepted by the line
    ;! from center to C
    ;! If frc < 0, extrapolation, but it is limited to |c| < (180-c0)/2.0
    tmp1 = sin(c)
    frc = tmp1 / (sin(c0 - c) + tmp1)

    ;! Work in Cartesian Coordinate
    tmp1 = frc * cos(lat2)
    tmp2 = 1.0 - frc

    x = tmp2 * cos(lat1) + tmp1 * cos(gam0)
    y = tmp1 * sin(gam0)
    z = tmp2 * sin(lat1) + frc * sin(lat2)

    gam = atan(y/x)                          ;! -90 < gam < 90
    if (frc >= 0) then begin
       if (gam < 0) then gam = gam + pi           ;! 0 <= gam <= 180
    endif else begin
       if (gam > 0) then gam = gam - pi           ;! -180 <= gam <= 0
    endif
    gam = gamsign * gam                      ;! get correct sign
    lon = gam + lon1

    theta = atan (sqrt(x**2 + y**2) / z)     ;! -90 < theta < 90
    if (theta < 0) then theta = theta + pi        ;! 0 <= theta <= 180
    lat = pihalf - theta                     ;! convert to latitude
    return
end

