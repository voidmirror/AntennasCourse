clear all;
clc;

% Programm of side lobe error calculation for circular aperture for
% different reference distance
tic
% Design parameters
a = 10;     % Aperture radius, m
b = 5;


lam = a/10;
R0 = 8*a^2/lam;
R = [R0*0.1:0.05*R0:1*R0, (R0 + R0*0.2):0.2*R0:2*R0];
Th = [0, 4.8]*pi/180; %, 7.8, 10.9, 13.9, 16.9]*pi/180; % set angle of required side lobe
N_Th = length(Th);
r = 0:0.5:b;
phi = 0:0.2:2*pi;
c = 3e+8;   % speed of light, m/sec

N_R = length(R);
N_FFT = 8192*2;
T = 100/c*lam;
d_t = T/(N_FFT - 1);

d_f = 1/N_FFT/d_t;
m = 1:1:N_FFT/2;
[F, p] = min(abs(d_f*m - 3e+8/lam));

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
%             plot((t - R/c)*1e+9, E_e(k, :));
        end
        E_b1(i, :) = trapz(r, E_e.*(r'*ones(1,N_i)));
    end
    
    E_b(:) = trapz(phi, E_b1);
    E_f = fft(E_b, N_FFT);
    E_ff(n, j) = E_f(p);
    j
end
n
end
% grid;
if (length(R) == 1)
    figure(1);
    plot((t - R/c)*1e+9, E_b); grid
    xlabel('time, nsec');
    ylabel('E field amplitude');
    title('E field time response for a/\lambda = 10. R_0 = 8*a^2/\lambda');
%     figure(2);
%     d_f = 1/N_FFT/d_t;
%     k = 1:1:N_FFT/2;
%     plot(d_f*k, 20*log10(abs(E_f(1:end/2))/max(abs(E_f)))); grid
else
%     figure(3);
%     plot(d_f*k, 20*log10(abs(E_f(10, 1:end/2))/max(abs(E_f(10, 1:end/2))))); grid
%     figure(4);
%     plot(R/R0, 20*log10(abs(E_f(2, :, p))/max(abs(E_f(1, :, p))))); grid
%     xlabel('R/R0');
%     ylabel('Amplitude, dB');
%     title('Radiation pattern of circular aperture for a/\lambda = 10. R_0 = 8*a^2/\lambda');
    
    load CA_Ideal_Field_Th_80_a_lam10.mat;
    
    Er_1_dB = 20*log10(abs(E_ff(2, :))./abs(E_ff(1, :))) - 20*log10(abs(E_i_CA_Th_80(49))/max(abs(E_i_CA_Th_80)));
%     Er_2_dB = 20*log10(abs(E_ff(3, :))./abs(E_ff(1, :))) - 20*log10(abs(E_i_CA_Th_80(80))/max(abs(E_i_CA_Th_80)));
%     Er_3_dB = 20*log10(abs(E_ff(4, :))./abs(E_ff(1, :))) - 20*log10(abs(E_i_CA_Th_80(110))/max(abs(E_i_CA_Th_80)));
%     Er_4_dB = 20*log10(abs(E_ff(5, :))./abs(E_ff(1, :))) - 20*log10(abs(E_i_CA_Th_80(140))/max(abs(E_i_CA_Th_80)));
%     Er_5_dB = 20*log10(abs(E_ff(6, :))./abs(E_ff(1, :))) - 20*log10(abs(E_i_CA_Th_80(170))/max(abs(E_i_CA_Th_80)));
    
    figure(10);
    plot(R/R0, Er_1_dB); grid
    xlabel('R/R0');
    ylabel('Amplitude error, dB');
    title('Error of 1^s^t side lobe of CA RP for a/\lambda = 10. R_0 = 8*a^2/\lambda. \Theta  = 4.8^0');
    
%     figure(11);
%     plot(R/R0, Er_2_dB); grid
%     xlabel('R/R0');
%     ylabel('Amplitude error, dB');
%     title('Error of 2^s^t side lobe of CA RP for a/\lambda = 10. R_0 = 8*a^2/\lambda. \Theta  = 7.85^0');
%     
%     figure(12);
%     plot(R/R0, Er_3_dB); grid
%     xlabel('R/R0');
%     ylabel('Amplitude error, dB');
%     title('Error of 3^s^t side lobe of CA RP for a/\lambda = 10. R_0 = 8*a^2/\lambda. \Theta  = 10.88^0');
%     
%     figure(13);
%     plot(R/R0, Er_4_dB); grid
%     xlabel('R/R0');
%     ylabel('Amplitude error, dB');
%     title('Error of 4^s^t side lobe of CA RP for a/\lambda = 10. R_0 = 8*a^2/\lambda. \Theta  = 13.9^0');
%     
%     figure(14);
%     plot(R/R0, Er_5_dB); grid
%     xlabel('R/R0');
%     ylabel('Amplitude error, dB');
%     title('Error of 5^s^t side lobe of CA RP for a/\lambda = 10. R_0 = 8*a^2/\lambda. \Theta  = 16.95^0');
end

toc


