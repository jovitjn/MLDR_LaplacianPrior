% MLDR beamformer with circular Laplacian prior

clc; clearvars; close all;

%Read the speech files-----------
[s1, fs] = audioread('61-70968-0000.flac');
s1 = s1(1:64000,:);

% s2 = audioread('237-126133-0000.flac');
s2 = audioread('61-70968-0003.flac');
s2 = s2(1:64000,:);

s = [s1, s2];
%--------------------------------------------------

N = 7; % number of microphonrs
uca_design = uca_2d(N, 0.1);
uca_design.element_positions(3,:) = 0; 
L_dim = [10 10 4];
r_positions = uca_design.element_positions';
[s1_pos(1), s1_pos(2), s1_pos(3)] = sph2cart(45*pi/180,(90-60)*pi/180,3);
% s1_pos = s1_pos + L_dim/2;
[s2_pos(1), s2_pos(2), s2_pos(3)] = sph2cart(135*pi/180,(90-60)*pi/180,3);
% s2_pos = s2_pos + L_dim/2;
c_vel = 340;
s_pos = [s1_pos; s2_pos];

figure
scatter3(r_positions(:,1), r_positions(:,2), r_positions(:,3)); xlim([0 L_dim(1)]); ylim([0 L_dim(2)]); zlim([0 L_dim(3)]); hold on;
scatter3(s1_pos(:,1), s1_pos(:,2), s1_pos(:,3))
scatter3(s2_pos(:,1), s2_pos(:,2), s2_pos(:,3))


%Microphone signal in the STFT domain------------------------------
OverlapLength = 256;
fftlength = 512;
window_length = 512; %corresponds to 640/16000 = 40 ms of data

for i = 1:2
[S(:,:,i), f, t] = stft(s(:,i),fs,'Window',hanning(window_length),'OverlapLength',OverlapLength,'FFTLength',fftlength, 'FrequencyRange','onesided');
end

S = permute(S,[3 2 1]);

SNRdB = 20; 
noise_std = 10^(-SNRdB/20);

for i = 1:length(f)
    Noise = (randn(N, length(t)) + 1i * randn(N, length(t)))/sqrt(2);
    s_norm = norm(S(:,:,i));
    p(:,:,i) = exp(-1i*2*pi/340*f(i)*s_pos*r_positions').' * S(:,:,i); % + noise_std*s_norm*Noise/norm(Noise); %(figure out a way to add diffuse noise)
end


delta = 1e-6;
epsilon = 1e-8;


% MLDR Algorithm ------------------------------------------
for ff = 1:length(f)
    ff
    y = p(:,:,ff); 
    R = (y*y')/length(t);
    a = exp(-1i*2*pi/340*f(ff)*s_pos(1,:)*r_positions').';
    w_zero = inv(R + delta*eye(N))*a/(a'*inv(R + delta*eye(N))*a);
    lambda = (abs(w_zero'*y)).^2;
    lambda = max(lambda, epsilon);
    %MLDR_LP beamformer
    for sub_itr = 1:3
        cvx_begin quiet
            variables bb(2*length(t), 1) w(N, 1)
            minimize( norm(bb, 1) )
            subject to
                bb >= 0;
                (abs(real(w'*y))).' <= sqrt(lambda')/2 .*(bb(1:2:end));
                (abs(imag(w'*y))).' <= sqrt(lambda')/2 .*(bb(2:2:end));
                w'*a == 1;
        cvx_end
        lambda = (abs(real(w'*y)) + abs(imag(w'*y))).^2;
        lambda = max(lambda, epsilon);
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

% %speech of first speaker
% sound(s1, fs) 
% 
% %speech of second speaker
% sound(s2, fs) 
% 
% %signal received at the first microphone
% sound(xn(1, 1:64000),fs) 
% 
% %signal at the output of the beamformer (should match the desired speaker)
% sound(real(s_est_MLDR(1:64000)),fs) 