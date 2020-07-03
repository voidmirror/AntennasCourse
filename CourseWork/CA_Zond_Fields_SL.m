clear all;
clc;

% Programm of Side lobes error calculation for circular aperture
% measured with aperture zond
tic
% Design parameters
a = 10;     % Aperture radius, m
b = 0.55*a;  % Zond radius


lam = a/10;     % Wave length in free space
R0 = 8*a^2/lam; % Reference distance
R = 1*R0/2;       % Distance to zond
n_SL = [1, 49, 79, 110, 140, 170, 201, 233, 265, 298, 333, 368,...
   406, 445, 488, 534, 586, 646, 721];  % Indexes of side lobes in far
%    field RP for a/lam = 10

% [1, 48, 78, 107, 136, 165, 194, 223, 252, 281, 310, 338, 367, 396, 425, 454,...
%     483, 512, 540, 569, 598, 627, 656, 685, 714, 743, 772, 801, 830, 859, 888, 917, 946, 976];

Th = (n_SL - 1)*0.1*pi/180;         % Corresponding side lobe angle
N_Th = length(Th);
r = 0:0.5:b;        % Integration parameter r
phi = 0:0.2:2*pi;   % Integration parameter phi
c = 3e+8;   % speed of light, m/sec

N_R = length(R);
N_FFT = 8192*4;   % FFT size
T = 100/c*lam;   % Time interval
d_t = T/(N_FFT - 1);    % Sample time

E_f = zeros(N_Th, N_R, N_FFT);
for n = 1:N_Th
for j = 1:N_R
    ro = R(j)*sin(Th(n));   % Vector R projection to aperture plane, m
    z = R(j)*cos(Th(n));     % Distance from observation point to aperture plane, m
    zb = z - r'*sin(phi)*sin(Th(n));
    rob = ro + r'*sin(phi)*cos(Th(n));
    t = R(j)/c - T/2:d_t:R(j)/c + T/2; % Time in seconds
    N_i = length(t);
    N_r = length(r);
    N_phi = length(phi);
    
    E_b1 = zeros(N_phi, N_i);
    E_b = zeros(1, N_i);
    
    for i = 1:N_phi
        E_e = zeros(N_r, N_i);
        for k = 1:N_r
            B = sqrt((c*t).^2 - zb(k, i)^2); % Radius of G curve, m
            if abs(rob(k, i)) <= a
                i2 = find((c*t >= zb(k, i)) & (c*t < sqrt(zb(k, i)^2 + (a - abs(rob(k, i)))^2)));
                E_e(k, i2) = zb(k, i)^2./((c*t(i2)).^2);
            end
                i3 = find((c*t >= sqrt(zb(k, i)^2 + (a - abs(rob(k, i)))^2)) & (c*t < sqrt(zb(k, i)^2 + (a + abs(rob(k, i)))^2)));
                E_e(k, i3) = zb(k, i)^2/pi./((c*t(i3)).^2).*acos((-a^2 + abs(rob(k, i))^2 + B(i3).^2)./(2*abs(rob(k, i))*B(i3)));
        end
        E_b1(i, :) = trapz(r, E_e.*(r'*ones(1,N_i)));
    end
    
    E_b(:) = trapz(phi, E_b1);
    E_f(n, j, :) = fft(E_b, N_FFT);
%     j
end
n
end

d_f = 1/N_FFT/d_t;
k = 1:1:N_FFT/2;
[F, p] = min(abs(d_f*k - 3e+8/lam));
    
load CA_Ideal_Field_Th_80_a_lam10.mat;
       
Er_SL_dB = (20*log10(abs(E_f(:, 1, p))./abs(E_f(1, 1, p)))).' - 20*log10(abs(E_i_CA_Th_80(n_SL))./max(abs(E_i_CA_Th_80)));
Er_SL_dB(1) = 0;

figure(2);
plot(0:length(n_SL) - 1, Er_SL_dB); grid
xlabel('Side lobe number');
ylabel('Amplitude error, dB');
title('Error of side lobes of CA RP for a/\lambda = 10. R = R_0. R_0 = 8*a^2/\lambda.');

toc


