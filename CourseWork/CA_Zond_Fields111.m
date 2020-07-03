clear all;
clc;

% Programm of time-space electic fields or radiation pattern calculation
% for circular aperture measured with aperture zond

tic
% Design parameters
a = 10;     % Analyzed aperture radius, m
b = 0.1;     % Zond aperture radius, m

lam = a/10;     % Wavelength in free space
R0 = 8*a^2/lam; % Reference distance
R = 10*R0/1;    % Distance from analyzed aperture to zond


Th = (0:0.05:10)*pi/180; % Theta angle. Can be set to one value (result will be a time field on zond output)
                        % or to array (result will be 2D RP at a given distance and frequency, and
                        % errors in RP compared with far field pattern)
r = 0:0.1:b;        % Integration parameter r
phi = 0:0.2:2*pi;   % Integration parameter phi
c = 3e+8;   % Speed of light, m/sec

N_Th = length(Th);
N_FFT = 8192*2;       % Number of FFT points (should be varied for faster calculation without accuracy degradation)
T = 100/c*lam;       % Time interval for FFT
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
                E_e(k, i2)=1.;  %zb(k, i)^1./((c*t(i2)).^1);
            end
                i3 = find((c*t >= sqrt(zb(k, i)^2 + (a - abs(rob(k, i)))^2)) & (c*t < sqrt(zb(k, i)^2 + (a + abs(rob(k, i)))^2)));
                E_e(k, i3) = acos((-a^2 + abs(rob(k, i))^2 + B(i3).^2)./(2*abs(rob(k, i))*B(i3)));
%             plot((t - R/c)*1e+9, E_e(k, :));
        end
        E_b1(i, :) = trapz(r, E_e.*(r'*ones(1,N_i)));
    end
    
    E_b(:) = trapz(phi, E_b1);
    E_f(j, :) = fft(E_b, N_FFT);
    j
end
% grid;
if (length(Th) == 1)
    figure(1);
    plot((t - R/c)*1e+9, E_b/max(E_b)); grid
    xlabel('time, nsec');
    ylabel('E field amplitude');
    title('E field time response for a/\lambda = 10. R_0 = 8*a^2/\lambda. R = R_0/4');
%     figure(2);
%     d_f = 1/N_FFT/d_t;
%     k = 1:1:N_FFT/2;
%     plot(d_f*k, 20*log10(abs(E_f(1:end/2))/max(abs(E_f)))); grid
else
%     figure(3);
%     mesh(t*1e9, Th*180/pi, E_b);
%     xlabel('time, nsec');
%     ylabel('\Theta, deg.');
%     zlabel('E field amplitude');
    
    d_f = 1/N_FFT/d_t;
    k = 1:1:N_FFT/2;
    figure(3);
    [F, p] = min(abs(d_f*k - 3e+8/lam));
    plot(2*a/lam*sin(Th), 20*log10(abs(E_f(:, p))/max(abs(E_f(:, p))))); grid
    xlabel('2*a/\lambda*sin(\Theta)');
    ylabel('Amplitude, dB');
    title('Radiation pattern of circular aperture for a/\lambda = 10. R_0 = 8*a^2/\lambda');
%     figure(5);
%     load CA_Ideal_Field_Th_a_lam40.mat;
%     Er = (abs(E_f(:, p))/max(abs(E_f(:, p)))).^2 - (abs(E_i_CA(:))/max(abs(E_i_CA))).^2;
% %     i9 = find(Th >= atan((a + b)/R));
%     plot(2*a/lam*sin(Th), Er); grid
%     xlabel('2*a/\lambda*sin(\Theta)');
%     ylabel('Amplitude error, E_b^2 - E_id^2');
%     title('Error of circular aperture RP for a/\lambda = 10. R = R_0/4. R_0 = 8*a^2/\lambda');
%     figure(6);
%     Er_dB = 20*log10(abs(E_f(:, p))/max(abs(E_f(:, p)))) - 20*log10(abs(E_i_CA(:))/max(abs(E_i_CA)));
%     plot(2*a/lam*sin(Th), Er_dB); grid
%     xlabel('2*a/\lambda*sin(\Theta)');
%     ylabel('Amplitude error, dB');
%     title('Error of circular aperture RP for a/\lambda = 10. R = R_0/4. R_0 = 8*a^2/\lambda');
%     
end

toc


