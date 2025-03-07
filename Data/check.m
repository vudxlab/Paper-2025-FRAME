%% plot check
%% Plot diagram
load('TrussO0.mat');
for i=1:length(acceleration(:,1))
     figure;
     plot(acceleration(i,:));
     title(["acc" i]);
     xlabel("samples");
     ylabel("acc");
end
Number of data points
acc = acceleration(1,:);
Fs = 1/0.02;
N = length(acc);



% Apply the FFT to the acceleration data
fft_acc = fft(acc);

% Calculate the two-sided spectrum and then the single-sided spectrum
P2 = abs(fft_acc / N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Create a frequency vector
f = Fs*(0:(N/2))/N;

% Plot the single-sided amplitude spectrum
plot(f, P1) 
title('Single-Sided Amplitude Spectrum of Acceleration')
xlabel('Frequency (Hz)')
ylabel('|P1(f)|')