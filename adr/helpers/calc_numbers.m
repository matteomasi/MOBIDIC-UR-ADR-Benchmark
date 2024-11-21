%% Calculation of Cr, Pe and F numbers
function out = calc_numbers(dt, dx, v, D)
   out.Cr = v.*dt./dx;
   out.Pe = v./D.*dx;
   out.F = D.*dt./(dx.^2);
end