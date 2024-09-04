clc;
clear;
close all;

% Parameters
bit_rates = [2, 3, 4, 5]; % bits/sec
amplitudes = [2, 3, 4, 5]; % Volts
sampling_frequencies = [12, 16, 20, 24]; % Hz
num_bits_per_packet = 10; % bits
rng(7,"twister") %% To define seed, here 7 is seed and twister is the algorithm

% Fixed parameters
fixed_bit_rate = 3; % bits per second
fixed_amplitude = 3; % Volts

% Number of packets for BER calculation
num_packets = 10;
% Number of bits per packet
num_bits = num_packets * num_bits_per_packet; % Ans 100

% Generate random bits once
data = randi([0 1], 1, num_bits); 
% Uniformly distributed pseudorandom integers
% [0 1] is interval here and 1*num_bits array

% Function to calculate theoretical BER for FSK using Q-function
theoretical_BER = @(SNR) qfunc(sqrt(2 * SNR));

% Preallocate matrices
measured_BER_bit_rate = zeros(length(bit_rates), length(sampling_frequencies));
theoretical_BER_bit_rate = zeros(length(bit_rates), length(sampling_frequencies));
PER_bit_rate = zeros(length(bit_rates), length(sampling_frequencies));

% Loop for varying bit rates and sampling frequencies
fprintf('\n--- Varying Bit Rates and Sampling Frequencies ---\n');
for i = 1:length(bit_rates)
    for j = 1:length(sampling_frequencies)
        
        fs = sampling_frequencies(j);
        freq_sep = fs / 4; % Ensure freq_sep is less than or equal to fs/2
        nsamp = round(fs / bit_rates(i));
        
        % FSK Modulation
        fsk_signal = fskmod(data, 2, freq_sep, nsamp, fs);
        fsk_signal = fixed_amplitude * fsk_signal;
        
        % Add AWGN noise
        % At AWGN channel keeping SNR fixed
        SNR = 5; % in dB
        noisy_signal = awgn(fsk_signal, SNR, 'measured');
        
        % FSK Demodulation
        demodulated_data = fskdemod(noisy_signal, 2, freq_sep, nsamp, fs);
        
        threshold = .5;
        demodulated_data = (demodulated_data > threshold);
        
        % Calculate BER
        [num_errors, ber] = biterr(data, demodulated_data);
        measured_BER_bit_rate(i, j) = ber;
        SNR_linear = 10^(SNR/10);
        theoretical_BER_bit_rate(i, j) = theoretical_BER(SNR_linear);
        
        % Calculate PER
        packet_errors = sum(sum(reshape(data, num_bits_per_packet, num_packets) ~= ...
                         reshape(demodulated_data, num_bits_per_packet, num_packets), 1) > 0);
        PER_bit_rate(i, j) = packet_errors / num_packets;
        
        % Print results
        fprintf('Bit Rate: %d bits/sec, Sampling Frequency: %d Hz, Measured BER: %.4e, Theoretical BER: %.4e, PER: %.4e\n', ...
            bit_rates(i), sampling_frequencies(j), measured_BER_bit_rate(i, j), theoretical_BER_bit_rate(i, j), PER_bit_rate(i, j));
        
        fprintf('    Total Bit Errors: %d, Total Packet Errors: %d\n', ...
            num_errors, packet_errors);
    end
end

% Plot BER vs Bit Rate for different sampling frequencies
figure(1);
subplot(221)
hold on;
grid on;
colors = 'bgrm'; % Define different colors for different sampling frequencies
markers = 'ox+*'; % Define different markers for different sampling frequencies
legend_entries = cell(1, 2 * length(sampling_frequencies));
k = 1;

for j = 1:length(sampling_frequencies)
    plot(bit_rates, measured_BER_bit_rate(:, j), [colors(j) markers(j) '-'], 'LineWidth', 2);
    legend_entries{k} = sprintf('Measured BER (Fs = %d Hz)', sampling_frequencies(j));
    k = k + 1;
    plot(bit_rates, theoretical_BER_bit_rate(:, j), [colors(j) markers(j) '--'], 'LineWidth', 2);
    legend_entries{k} = sprintf('Theoretical BER (Fs = %d Hz)', sampling_frequencies(j));
    k = k + 1;
end

xlabel('Bit Rate (bits/sec)');
ylabel('BER');
title('BER vs Bit Rate for Different Sampling Frequencies');
legend(legend_entries, 'Location', 'best');
hold off;



% Preallocate matrices
measured_BER_amplitude = zeros(length(amplitudes), length(sampling_frequencies));
theoretical_BER_amplitude = zeros(length(amplitudes), length(sampling_frequencies));
PER_amplitude = zeros(length(amplitudes), length(sampling_frequencies));
SNR_amplitude= zeros(length(amplitudes), length(sampling_frequencies));
noise_amplitude = 2;

fprintf('\n--- Varying Amplitudes ---\n');
for i = 1:length(amplitudes)
    for j = 1:length(sampling_frequencies)
        bit_rate = fixed_bit_rate;
        fs = sampling_frequencies(j);
        freq_sep = fs / 4; % Ensure freq_sep is less than or equal to fs/2
        nsamp = round(fs / bit_rate);
        
        % FSK Modulation
        fsk_signal = fskmod(data, 2, freq_sep, nsamp, fs);
        fsk_signal = amplitudes(i) * fsk_signal;
        SNR_amplitude(i, j) = 20 * log10(amplitudes(i) / noise_amplitude);
        
        % Add AWGN noise
        noisy_signal = awgn(fsk_signal, SNR_amplitude(i, j), 'measured');
        
        % FSK Demodulation
        demodulated_data = fskdemod(noisy_signal, 2, freq_sep, nsamp, fs);

        threshold = .5;
        demodulated_data = (demodulated_data > threshold);
        
        % Calculate BER
        [num_errors, ber] = biterr(data, demodulated_data);
        measured_BER_amplitude(i, j) = ber;
        SNR_linear = 10^(SNR_amplitude(i, j) / 10);
        theoretical_BER_amplitude(i, j) = theoretical_BER(SNR_linear);
        
        % Calculate PER
        packet_errors = sum(sum(reshape(data, num_bits_per_packet, num_packets) ~= ...
                         reshape(demodulated_data, num_bits_per_packet, num_packets), 1) > 0);
        PER_amplitude(i, j) = packet_errors / num_packets;
        
        % Print results
        fprintf('Amplitude: %d V, Sampling Frequency: %d Hz, SNR: %.3f, Measured BER: %.4e, Theoretical BER: %.4e, PER: %.4e\n', ...
            amplitudes(i), sampling_frequencies(j), SNR_amplitude(i, j), measured_BER_amplitude(i, j), theoretical_BER_amplitude(i, j), PER_amplitude(i, j));
        
        fprintf('    Total Bit Errors: %d, Total Packet Errors: %d\n', ...
            num_errors, packet_errors);
    end
end

% Plot BER vs Amplitude for different sampling frequencies
subplot(222)
hold on;
grid on;
colors = 'bgrm'; % Define different colors for different sampling frequencies
markers = 'ox+*'; % Define different markers for different sampling frequencies
legend_entries = cell(1, 2 * length(sampling_frequencies));
k = 1;

for j = 1:length(sampling_frequencies)
    plot(amplitudes, measured_BER_amplitude(:, j), [colors(j) markers(j) '-'], 'LineWidth', 2);
    legend_entries{k} = sprintf('Measured BER (Fs = %d Hz)', sampling_frequencies(j));
    k = k + 1;
    plot(amplitudes, theoretical_BER_amplitude(:, j), [colors(j) markers(j) '--'], 'LineWidth', 2);
    legend_entries{k} = sprintf('Theoretical BER (Fs = %d Hz)', sampling_frequencies(j));
    k = k + 1;
end

xlabel('Amplitude (V)');
ylabel('BER');
title('BER vs Amplitude for Different Sampling Frequencies');
legend(legend_entries, 'Location', 'best');
hold off;


rng(6,"twister");
data = randi([0 1], 1, num_bits); 

% Preallocate matrices
measured_BER_fs = zeros(length(sampling_frequencies), length(amplitudes));
theoretical_BER_fs = zeros(length(sampling_frequencies), length(amplitudes));
PER_fs = zeros(length(sampling_frequencies), length(amplitudes));

fprintf('\n--- Varying Sampling Frequencies and Amplitudes ---\n');
for i = 1:length(sampling_frequencies)
    for j = 1:length(amplitudes)
        bit_rate = fixed_bit_rate;
        fs = sampling_frequencies(i);
        freq_sep = fs / 4; % Ensure freq_sep is less than or equal to fs/2
        nsamp = round(fs / bit_rate);
        
        % FSK Modulation
        fsk_signal = fskmod(data, 2, freq_sep, nsamp, fs);
        fsk_signal = amplitudes(j) * fsk_signal;
        
        % Add AWGN noise
        SNR = 2; % in dB
        noisy_signal = awgn(fsk_signal, SNR, 'measured');
        
        % FSK Demodulation
        demodulated_data = fskdemod(noisy_signal, 2, freq_sep, nsamp, fs);

        threshold = .5;
        demodulated_data = (demodulated_data > threshold);
        
        % Calculate BER
        [num_errors, ber] = biterr(data, demodulated_data);
        measured_BER_fs(i, j) = ber;
        SNR_linear = 10^(SNR / 10);
        theoretical_BER_fs(i, j) = theoretical_BER(SNR_linear);
        
        % Calculate PER
        packet_errors = sum(sum(reshape(data, num_bits_per_packet, num_packets) ~= ...
                         reshape(demodulated_data, num_bits_per_packet, num_packets), 1) > 0);
        PER_fs(i, j) = packet_errors / num_packets;
        
        % Print results
        fprintf('Sampling Frequency: %d Hz, Amplitude: %d V, Measured BER: %.4e, Theoretical BER: %.4e, PER: %.4e\n', ...
            sampling_frequencies(i), amplitudes(j), measured_BER_fs(i, j), theoretical_BER_fs(i, j), PER_fs(i, j));
        
        fprintf('    Total Bit Errors: %d, Total Packet Errors: %d\n', ...
            num_errors, packet_errors);
    end
end

% Plot BER vs Sampling Frequency for different amplitudes
subplot(223)
hold on;
grid on;
colors = 'bgrm'; % Define different colors for different amplitudes
markers = 'ox+*'; % Define different markers for different amplitudes
legend_entries = cell(1, 2 * length(amplitudes));
k = 1;

for j = 1:length(amplitudes)
    plot(sampling_frequencies, measured_BER_fs(:, j), [colors(j) markers(j) '-'], 'LineWidth', 2);
    legend_entries{k} = sprintf('Measured BER (Amplitude = %d V)', amplitudes(j));
    k = k + 1;
    plot(sampling_frequencies, theoretical_BER_fs(:, j), [colors(j) markers(j) '--'], 'LineWidth', 2);
    legend_entries{k} = sprintf('Theoretical BER (Amplitude = %d V)', amplitudes(j));
    k = k + 1;
end

xlabel('Sampling Frequency (Hz)');
ylabel('BER');
title('BER vs Sampling Frequency for Different Amplitudes');
legend(legend_entries, 'Location', 'best');
hold off;


rng(1,"twister");
data = randi([0 1], 1, num_bits); 
% Preallocate matrices
SNR_varied = [0, 2, 4, 6, 8, 10, 12]; % SNR in dB
measured_BER_SNR = zeros(length(SNR_varied), length(sampling_frequencies));
theoretical_BER_SNR = zeros(length(SNR_varied), length(sampling_frequencies));
PER_SNR = zeros(length(SNR_varied), length(sampling_frequencies));

fprintf('\n--- Varying SNR and Sampling Frequencies ---\n');
for i = 1:length(SNR_varied)
    for j = 1:length(sampling_frequencies)
        bit_rate = fixed_bit_rate; 
        fs = sampling_frequencies(j);
        freq_sep = fs / 4; % Ensure freq_sep is less than or equal to fs/2
        nsamp = round(fs / bit_rate);
        
        % FSK Modulation
        fsk_signal = fskmod(data, 2, freq_sep, nsamp, fs);
        fsk_signal = fixed_amplitude * fsk_signal;
        
        % Add AWGN noise
        noisy_signal = awgn(fsk_signal, SNR_varied(i), 'measured');
        
        % FSK Demodulation
        demodulated_data = fskdemod(noisy_signal, 2, freq_sep, nsamp, fs);

        threshold = .5;
        demodulated_data = (demodulated_data > threshold);
        
        % Calculate BER
        [num_errors, ber] = biterr(data, demodulated_data);
        measured_BER_SNR(i, j) = ber;
        SNR_linear = 10^(SNR_varied(i) / 10);
        theoretical_BER_SNR(i, j) = theoretical_BER(SNR_linear);
        
        % Calculate PER
        packet_errors = sum(sum(reshape(data, num_bits_per_packet, num_packets) ~= ...
                         reshape(demodulated_data, num_bits_per_packet, num_packets), 1) > 0);
        PER_SNR(i, j) = packet_errors / num_packets;
        
        % Print results
        fprintf('SNR: %d dB, Sampling Frequency: %d Hz, Measured BER: %.4e, Theoretical BER: %.4e, PER: %.4e\n', ...
            SNR_varied(i), sampling_frequencies(j), measured_BER_SNR(i, j), theoretical_BER_SNR(i, j), PER_SNR(i, j));
        
        fprintf('    Total Bit Errors: %d, Total Packet Errors: %d\n', ...
            num_errors, packet_errors);
    end
end

% Plot BER vs SNR for different sampling frequencies
subplot(224)
hold on;
grid on;
colors = 'bgrm'; % Define different colors for different sampling frequencies
markers = 'ox+*'; % Define different markers for different sampling frequencies
legend_entries = cell(1, 2 * length(sampling_frequencies));
k = 1;

for j = 1:length(sampling_frequencies)
    plot(SNR_varied, measured_BER_SNR(:, j), [colors(j) markers(j) '-'], 'LineWidth', 2);
    legend_entries{k} = sprintf('Measured BER (Fs = %d Hz)', sampling_frequencies(j));
    k = k + 1;
    plot(SNR_varied, theoretical_BER_SNR(:, j), [colors(j) markers(j) '--'], 'LineWidth', 2);
    legend_entries{k} = sprintf('Theoretical BER (Fs = %d Hz)', sampling_frequencies(j));
    k = k + 1;
end

xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR for Different Sampling Frequencies');
legend(legend_entries, 'Location', 'best');
hold off;

bit_duration = 1 / fixed_bit_rate;
t = 0:1/sampling_frequencies(4):(bit_duration * num_bits - 1/sampling_frequencies(4));


fc=9;
figure(2)
subplot(311)
stairs(0:bit_duration:bit_duration*length(data), [data data(end)], 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Amplitude');
title('Input Generated Random Data Signal');
ylim([-0.2 1.2]);
xlim([0,30]);


subplot(313)
stairs(0:bit_duration:bit_duration*length(data), [demodulated_data demodulated_data(end)], 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Amplitude');
title('Output Data Signal');
ylim([-0.2 1.2]);
xlim([0,30]);

subplot(312)
fsk_signal_mod = real(fsk_signal) .* cos(2*pi*fc*t) - imag(fsk_signal) .* sin(2*pi*fc*t);
plot(t, fsk_signal_mod);
xlabel('Time (s)');
ylabel('Amplitude');
title('FSK Signal');
ylim([-3.5 3.5]);
xlim([0,30]);
grid on;
