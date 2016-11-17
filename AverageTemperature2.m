function [radius,Tavg, Tc, q, qtotal, qavg] = AverageTemperature2(source)
%
% returns the average radial temperature from picked_mode file
%
%
 load picked_mode.dat
 
 
% retrieve radial levels
nr = 97;
ilevel = picked_mode(:,3);
nentry = length(ilevel);    % total number of entries in file

% retrieve radius
radius = picked_mode(1:nr,4);   % get the radial grid

% retrieve temperature array
temp = picked_mode(:,7);       % get temperature
flux = picked_mode(:,18);      % get heat flux

% initialize and evaluate average temperature (use ilevel to define radial position)
Tavg(1:nr) = 0;
qavg(1:nr) = 0;
for k = 1 : nentry
    Tavg(ilevel(k)) = Tavg(ilevel(k)) + temp(k);
    qavg(ilevel(k)) = qavg(ilevel(k)) + flux(k);
end

navg = nentry / nr
Tavg = Tavg / navg;
qavg = qavg / navg;

% evaluate temperature gradient  (you probably don’t need this)
dTdr = Derivative(radius,Tavg);
for k = 1 : nr
   q(k) = -radius(k)^2 * dTdr(k);    (convected heat flux)
end

% total heat flux  (heat flux at bottom minus sink (or negative source)
qtotal(1) = q(1);
for k = 2:nr
    qtotal(k) = qtotal(1) + source*(radius(k)^3 - radius(1)^3)/3;
end

% Trms
Trms = 0;
for i = 1 : nr-1
    ravg = 0.5 * ( radius(i) + radius(i+1) );
    dr = radius(i+1) - radius(i);
    Tbar = 0.5 * (Tavg(i) + Tavg(i+1));
    Trms = Trms + 4 * pi * ravg^2 * dr * Tbar^2;
end
V = 4 * pi * (radius(end)^3 - radius(1)^3)/3;
Trms/V

% conductive temperature
cdim = radius(1) / radius(nr);
qb2 = source * radius(nr)^2 / 6
r2 = radius.^2;
tmp = (1 - cdim);
c1 = radius(1)*(1 - qb2*(1 - cdim^2))/(1 - cdim);
c2 = qb2 - c1/radius(nr);

Tc = c1./radius + c2 - (source/6).*r2;

end

    



 
 