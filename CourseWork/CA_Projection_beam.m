clear all;
clc;

% This programm calculates projection beam width for circular aperture

a = 10;     % Aperture radius, m
sc = 0.01; % Scale parameter
b = a*sc;   % Zond aperture radius

lam = a/100;    % Wavelength in free space
R0 = 8*a^2/lam; % Standard far zone griteria, m
R = R0*0.1:0.1*R0:2*R0; % Distance to zond

C = a*R./(R.^2 + b^2);

Th_pr = asin(C + sqrt(C.^2 - (a^2 - b^2)./(R.^2 + b^2))); % Upper boundary of projection beam

figure(1);  % Projection beam in degrees
plot(R/R0, Th_pr*180/pi); grid
xlabel('R/R0');
ylabel('Projection beam width, deg.');


figure(2);  % Normalized projection beam
plot(R/R0, 2*a/lam*sin(Th_pr)); grid
xlabel('R/R0');
ylabel('Projection beam width, 2*a/\lambda*sin(\Theta_p_r)');

