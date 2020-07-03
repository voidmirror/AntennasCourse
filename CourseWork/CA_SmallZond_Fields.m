clear all;
clc;

% Programm of time-space electic fields or radiation pattern calculation
% for circular aperture measured with infinitely small zond

tic
% Design parameters
a = 10;         % Aperture radius, m


lam = a/10;     % Wavelength in free space
R0 = 8*a^2/lam; % Reference distance
R = 10*R0/1;     % Distance from analyzed aperture to zond
c = 3e+8;       % speed of light, m/sec

Th = (0:0.1:80)*pi/180; % Theta angle. Can be set to one value (result will be a time response on zond output)
                 % or to array (result will be 2D RP at a given distance and frequency)
N_Th = length(Th);
N_FFT = 8192*64;  % Number of FFT points (should be decreased for faster calculation with accuracy degradation)
T = 100/c*lam;    % Time interval for FFT
d_t = T/(N_FFT - 1); % Sample time, sec
t = R/c - T/2:d_t:R/c + T/2; % Time in seconds normalized to R/c
N_i = length(t);
E_e = zeros(1, N_i); % Time response array

d_f = 1/N_FFT/d_t;
m = 1:1:N_FFT/2;
for k = 1:N_Th
        
    ro = R*sin(Th(k));   % Vector R projection to aperture plane, m
    z = R*cos(Th(k));     % Distance from observation point to aperture plane, m
    B = sqrt((c*t).^2 - z^2); % Radius of G curve, m
    if abs(ro) <= a
        i2 = find((c*t >= z) & (c*t < sqrt(z^2 + (a - abs(ro))^2)));
        E_e(i2) = z^2./((c*t(i2)).^2);
    end
    i3 = find((c*t >= sqrt(z^2 + (a - abs(ro))^2)) & (c*t < sqrt(z^2 + (a + abs(ro))^2)));
    E_e(i3) = z^2/pi./((c*t(i3)).^2).*acos((-a^2 + abs(ro)^2 + B(i3).^2)./(2*abs(ro)*B(i3)));
    
    E_f = fft(E_e, N_FFT); % Wideband radiation pattern
    [F, p] = min(abs(d_f*m - 3e+8/lam));
    E_ff(k) = E_f(p);
    k
end

if (length(Th) == 1)
    figure(1); % Time response at a given angle
    plot((t - R/c)*1e+9, E_e/max(E_e)); grid
    xlabel('time, nsec');
    ylabel('E field amplitude');
    title('E field TR for a/\lambda = 10. R_0 = 8*a^2/\lambda. \Theta = 13.9^o');
%     figure(2); % Spectra at a given angle
%     d_f = 1/N_FFT/d_t;
%     k = 1:1:N_FFT/2;
%     plot(d_f*k, 20*log10(abs(E_f(1:end/2))/max(abs(E_f)))); grid
else
    
    figure(4); % Radiation pattern at a given frequency
    
    plot(2*a/lam*sin(Th), 20*log10(abs(E_ff)/max(abs(E_ff)))); grid
    xlabel('2*a/\lambda*sin(\Theta)');
    ylabel('Amplitude, dB');
    title('Radiation pattern of circular aperture for a/\lambda = 10. R_0 = 8*a^2/\lambda');
end
toc




