clear all;
clc;

% Programm of spectra of time response calculation
% for circular aperture measured with iaperture zond

tic
% Design parameters
a = 10;     % Analyzed aperture radius, m
b = 5;     % Zond aperture radius, m


R0 = 8*a^2/1
R = 10000;    % Distance from analyzed aperture to zond
Th = [0, 2.5]*pi/180; % Theta angle. Can be set to one value (result will be a time field on zond output)
                        % or to array (result will be 2D RP at a given distance and frequency, and
                        % errors in RP compared with far field pattern)
r = 0:0.5:b;        % Integration parameter r
phi = 0:0.2:2*pi;   % Integration parameter phi
c = 3e+8;   % Speed of light, m/sec

N_Th = length(Th);
N_FFT = 8192*32;       % Number of FFT points (should be varied for faster calculation without accuracy degradation)
T = 100/c;       % Time interval for FFT
d_t = T/(N_FFT - 1);% Sample time, sec

t = R/c - T/2:d_t:R/c + T/2; % Time in seconds normalized to R/c
N_i = length(t);

N_r = length(r);
N_phi = length(phi);

E_b1 = zeros(N_phi, N_i);
E_b = zeros(1, N_i);
E_f = zeros(N_Th, N_FFT);
% figure(5);
% hold on;
for j = 1:N_Th
    ro = R*sin(Th(j));   % Vector R projection to aperture plane, m
    z = R*cos(Th(j));     % Distance from observation point to aperture plane, m
    zb = z - r'*sin(phi)*sin(Th(j));
    rob = ro + r'*sin(phi)*cos(Th(j));
    for i = 1:N_phi
        E_e = zeros(N_r, N_i);
        for k = 1:N_r
            B = sqrt((c*t).^2 - zb(k, i)^2); % Radius of G curve, m
            if abs(rob(k, i)) <= a
                i2 = find((c*t >= zb(k, i)) & (c*t < sqrt(zb(k, i)^2 + (a - abs(rob(k, i)))^2)));
                E_e(k, i2) = zb(k, i)^2./((c*t(i2)).^2);
                i3 = find((c*t >= sqrt(zb(k, i)^2 + (a - abs(rob(k, i)))^2)) & (c*t < sqrt(zb(k, i)^2 + (a + abs(rob(k, i)))^2)));
                E_e(k, i3) = zb(k, i)^2/pi./((c*t(i3)).^2).*acos((-a^2 + abs(rob(k, i))^2 + B(i3).^2)./(2*abs(rob(k, i))*B(i3)));
            else
                i2 = find((c*t >= sqrt(zb(k, i)^2 + (a - abs(rob(k, i)))^2)) & (c*t <= sqrt(zb(k, i)^2 + (a + abs(rob(k, i)))^2)));
                E_e(k, i2) = zb(k, i)^2/pi./((c*t(i2)).^2).*acos((-a^2 + abs(rob(k, i))^2 + B(i2).^2)./(2*abs(rob(k, i))*B(i2)));
            end
%             plot((t - R/c)*1e+9, E_e(k, :));
        end
        E_b1(i, :) = trapz(r, E_e.*(r'*ones(1,N_i)));
    end
    
    E_b(:) = trapz(phi, E_b1);
    E_f(j, :) = fft(E_b, N_FFT);
    j
end
% grid;

% figure(1);
% plot(t*1e+9, E_e); grid
% xlabel('time, nsec');
% ylabel('E field amplitude');
% title('E field TR for CA. R = 1000 m, \Theta = 10^o, a = 10 m');

figure(2);
f = 1/T*(1:1:N_FFT);
k1 = find((f*2*a/c >= 10) & (f*2*a/c <= 100));
plot(f(k1)*2*a/c, 20*log10(abs(E_f(2, k1))./abs(E_f(1, k1)))); grid
xlabel('Normalized frequency, f*2*a/c');
ylabel('E field spectrum, dB');
title('E field TR for CA. R = 1000 m, \Theta = 10^o, a = 10 m');


toc


