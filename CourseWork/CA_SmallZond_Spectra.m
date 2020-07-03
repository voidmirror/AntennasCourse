clear all;
clc;

% Programm of spectra from time response calculation
% for circular aperture measured with infinitely small zond

tic
% Design parameters
a = 10;         % Aperture radius, m


R0 = 8*a^2/1    % Reference distance
R = 10000;     % Distance from analyzed aperture to zond
c = 3e+8;       % speed of light, m/sec

Th = [0, 2.5]*pi/180; % Theta angle. First value - 0, second value - considered angle
N_Th = length(Th);
N_FFT = 8192*32;  % Number of FFT points (should be decreased for faster calculation with accuracy degradation)
T = 100/c;    % Time interval for FFT
d_t = T/(N_FFT - 1); % Sample time, sec
t = R/c - T/2:d_t:R/c + T/2; % Time in seconds normalized to R/c
N_i = length(t);
E_e = zeros(1, N_i);

for k = 1:N_Th
    ro = R*sin(Th(k));   % Vector R projection to aperture plane, m
    z = R*cos(Th(k));    % Distance from observation point to aperture plane, m
    B = sqrt((c*t).^2 - z^2); % Radius of G curve, m
    
    if abs(ro) <= a
        i2 = find((c*t >= z) & (c*t < sqrt(z^2 + (a - abs(ro))^2)));
        E_e(i2) = z^2./((c*t(i2)).^2);
    end
    i3 = find((c*t >= sqrt(z^2 + (a - abs(ro))^2)) & (c*t < sqrt(z^2 + (a + abs(ro))^2)));
    E_e(i3) = z^2/pi./((c*t(i3)).^2).*acos((-a^2 + abs(ro)^2 + B(i3).^2)./(2*abs(ro)*B(i3)));

    E_f(k, :) = fft(E_e, N_FFT); % Wideband radiation pattern
    k
end

% figure(1); % Time response at a given angle
% plot(t*1e+9, E_e); grid
% xlabel('time, nsec');
% ylabel('E field amplitude');
% title('E field TR for CA. R = 1000 m, \Theta = 10^o, a = 10 m');

figure(2); % Spectra at a given angle
f = 1/T*(1:1:N_FFT);
k1 = find((f*2*a/c >= 10) & (f*2*a/c <= 100));
plot(f(k1)*2*a/c, 20*log10(abs(E_f(2, k1))./abs(E_f(1, k1)))); grid
xlabel('Normalized frequency, f*2*a/c');
ylabel('E field spectrum, dB');
title('E field TR for CA. R = 1000 m, \Theta = 10^o, a = 10 m');

toc



