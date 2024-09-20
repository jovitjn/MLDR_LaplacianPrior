% Assuming s1 and s_est_MLDR are already loaded and pre-processed as per your initial script

% Find the time delay
[crossCorr, lag] = xcorr(s_est_MLDR(1:64000), s1);
[~, idx] = max(abs(crossCorr));
delay = lag(idx);

% Align the signals based on the delay
if delay > 0
    aligned_s_est_MLDR = s_est_MLDR(delay+1:end);
    aligned_s1 = s1(1:end-delay);
elseif delay < 0
    aligned_s1 = s1(-delay+1:end);
    aligned_s_est_MLDR = s_est_MLDR(1:end+delay);
else
    aligned_s1 = s1;
    aligned_s_est_MLDR = s_est_MLDR;
end

% Ensure the signals are the same length for PESQ calculation
minLength = min(length(aligned_s1), length(aligned_s_est_MLDR));
aligned_s1 = aligned_s1(1:minLength);
aligned_s_est_MLDR = aligned_s_est_MLDR(1:minLength);

% Ensure the signals are real before saving
real_aligned_s1 = real(aligned_s1);
real_aligned_s_est_MLDR = real(aligned_s_est_MLDR);

% Save the real parts of the aligned signals to disk
audiowrite('reference.wav', real_aligned_s1, fs);
audiowrite('enhanced.wav', real_aligned_s_est_MLDR,fs)

% Display PESQ, cepstral distance, WSS, and LLR measures
fprintf('PESQ: %.4f\n', pesq('reference.wav', 'enhanced.wav'));
fprintf('Cepstral Distance: %.4f\n', comp_cep('reference.wav', 'enhanced.wav'));
fprintf('WSS: %.4f\n', comp_wss('reference.wav', 'enhanced.wav'));
fprintf('LLR: %.4f\n', comp_llr('reference.wav', 'enhanced.wav'));