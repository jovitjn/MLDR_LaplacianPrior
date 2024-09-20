% MLDR beamformer with complex Gaussian prior

clc; clearvars; close all;

%Read the speech files-----------
[s1, fs] = audioread('61-70968-0000.flac');
s1 = s1(1:64000,:);

% s2 = audioread('237-126133-0000.flac');
s2 = audioread('61-70968-0003.flac');
s2 = s2(1:64000,:);
%--------------------------------------------------

N = 7; % number of microphonrs
uca_design = uca_2d(N, 0.1);
uca_design.element_positions(3,:) = 0; 
L_dim = [10 10 4];
r_positions = uca_design.element_positions' + L_dim/2;
[s1_pos(1), s1_pos(2), s1_pos(3)] = sph2cart(45*pi/180,(90-60)*pi/180,3);
s1_pos = s1_pos + L_dim/2;
[s2_pos(1), s2_pos(2), s2_pos(3)] = sph2cart(135*pi/180,(90-60)*pi/180,3);
s2_pos = s2_pos + L_dim/2;
c_vel = 340;
s = [s1_pos; s2_pos];

% generating room impulse responses
for s_itr = 1:2
    beta = 0.2;                 % Reverberation time (s)
    n = 2048;                   % Number of samples
    h(:,:, s_itr) = 50*rir_generator(c_vel, fs, r_positions, s(s_itr, :), L_dim, beta, n);
end


SNRdB = 30; 
noise_std = 10^(-SNRdB/20); 

%Convolving the speech signals with the room impulse responses
for n = 1:N
    Noise = laprnd(1, length(conv(h(1,:,1), s1)));
    %Noise = randn(1, length(conv(h(1,:,1), s1)));
    xn(n,:) = conv(h(n,:,1), s1) + conv(h(n,:,2), s2);
    xn_norm = norm(xn(n,:));
    pn(n,:) = xn(n,:)+ noise_std*xn_norm*Noise/norm(Noise);

end

%Microphone signal in the STFT domain------------------------------
OverlapLength = 256;
fftlength = 512;
window_length = 512; %corresponds to 640/16000 = 40 ms of data

for n = 1:N
[x(:,:,n), f, t] = stft(pn(n,:),fs,'Window',hanning(window_length),'OverlapLength',OverlapLength,'FFTLength',fftlength, 'FrequencyRange','onesided');
end
x = permute(x,[3 2 1]);

% Adding noise--------------------
% SNRdB = 30; 
% noise_std = 10^(-SNRdB/20);
% 
% for i = 1:length(f)
%     Noise = (randn(N, length(t)) + 1i * randn(N, length(t)))/sqrt(2);
%     x_norm = norm(x(:,:,i));
%     p(:,:,i) = x(:,:,i) + noise_std*x_norm*Noise/norm(Noise); %(figure out a way to add diffuse noise)
% end

H = [];
for n = 1:N
[Hh, fh, th] = stft(h(n,:,1),fs,'Window',hanning(window_length),'OverlapLength',OverlapLength,'FFTLength',fftlength, 'FrequencyRange','onesided');
H = [H Hh(:,1)];
end
H = H.';

delta = 1e-6;
epsilon = 1e-8;

% MPDR Algorithm ------------------------------------------
for ff = 1:length(f)
    ff
    y = x(:,:,ff); 
    R = (y*y')/length(t);
    a = H(:,ff)/H(1,ff);
    w_zero = inv(R + delta*eye(N))*a/(a'*inv(R + delta*eye(N))*a);
    for sub_itr = 1:4
        R_t = zeros(N,N);
        for t_itr = 1:length(t)
            R_t = R_t + (y(:,t_itr)*y(:,t_itr)');
        end
        w = inv(R_t + delta*eye(N))*a/(a'*inv(R_t + delta*eye(N))*a);
    end
    x_est_MLDR(ff, :) = w'*y;
end
s_est_MLDR = istft(x_est_MLDR,fs,'Window',hanning(window_length),'OverlapLength',OverlapLength,'FFTLength',fftlength, 'FrequencyRange','onesided');

%----Plots-----
%STFT of the first speaker
figure
stft(s1,fs,'Window',hanning(window_length),'OverlapLength',OverlapLength,'FFTLength',fftlength);
set(gcf, 'color', 'w')
ylim([0 8])

%STFT of the second speaker
figure
stft(s2,fs,'Window',hanning(window_length),'OverlapLength',OverlapLength,'FFTLength',fftlength);
set(gcf, 'color', 'w')
ylim([0 8])

%STFT at the signal received at the first microphone
figure
stft(xn(1, 1:64000),fs,'Window',hanning(window_length),'OverlapLength',OverlapLength,'FFTLength',fftlength);
set(gcf, 'color', 'w')
ylim([0 8])

%STFT at the output of the beamformer (should match with the desired speech STFT)
figure
stft(real(s_est_MLDR(1:64000)),fs,'Window',hanning(window_length),'OverlapLength',OverlapLength,'FFTLength',fftlength);
set(gcf, 'color', 'w')
ylim([0 8])


%Uncomment these to play the sound

%speech of first speaker
%sound(s1, fs) 
% 
% %speech of second speaker
% sound(s2, fs) 
% 
% %signal received at the first microphone
% sound(xn(1, 1:64000),fs) 
% 
% %signal at the output of the beamformer (should match the desired speaker)
% sound(real(s_est_MLDR(1:64000)),fs) 


