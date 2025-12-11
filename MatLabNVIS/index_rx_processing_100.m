
clc;
clear;

addpath('src');

SAMPLE_RATE_Hz = 192000;

% Параметры зондирующего ЛЧМ-сигнала:

chirp_bandwidth_Hz = 100000;
reverse_chirp_divider = 32; % доля времени на "обратный" chirp

% Обработка без децимации

mas_CSS_etalon_IQ = read_cfloat32('test_iq_etalon/etalon_100kHz');

fft_len = length(mas_CSS_etalon_IQ);

fft_len_slow = 4096;

tau_range_ms = chirp_bandwidth_Hz * reverse_chirp_divider / (reverse_chirp_divider-1);

mas_tau_ms = (-fft_len/2+1 : 1 : fft_len/2) / tau_range_ms * 1000;
mas_t_sec = (1:1:fft_len_slow) * fft_len / SAMPLE_RATE_Hz;
mas_f_doppler_Hz = ((1:fft_len_slow).' - fft_len_slow/2 - 1) * SAMPLE_RATE_Hz / fft_len / fft_len_slow;

% mas_window = ones(fft_len, 1);
% mas_window = window(@blackmanharris, fft_len);
mas_window = [zeros(fft_len/8,1); window(@blackmanharris, 3*fft_len/4); zeros(fft_len/8,1)];

% *** Цифровой эксперимент:

% mas_Rx_IQ = read_cfloat32('data_rx\digital\css100_rx');

% *** Хохольский поворот, WinRadio:

% load('data_rx\2024-10-04_SZV\2024-10-04_15-25_LCHM_100k_11960_sr192kHz.mat'); % сигнал есть, но разрывы IQ на приёме
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% load('data_rx\2024-10-04_SZV\2024-10-04_21-05_LCHM_100k_6540_sr192kHz.mat'); % брызги
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% load('data_rx\2024-10-04_SZV\2024-10-04_21-45_LCHM_100k_5075_sr192kHz.mat'); % сигнал есть, но разрывы IQ на приёме
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% load('data_rx\2024-10-04_SZV\2024-10-04_22-12_LCHM_100k_5075_sr192kHz.mat'); % сигнал есть, но разрывы IQ на приёме
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% *** Малышево, 8 Rx:

load('data_rx\2024-10-04_Mal\2024-10-04_21-06_LCHM_100k_sr192kHz.mat'); % !!! интересно
symbol_sync_ratio = 1;
freq_shift_Hz = 0;

% load('data_rx\2024-10-04_Mal\2024-10-04_21-45_LCHM_100k_sr192kHz.mat'); % есть сигнал 
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% load('data_rx\2024-10-04_Mal\2024-10-04_22-12_LCHM_100k_sr192kHz.mat'); % есть сигнал 
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

i_Rx = 4;
mas_Rx_IQ = matr_Rx_IQ(:,i_Rx);
clear matr_Rx_IQ;


mas_Rx_IQ = my_freq_shift(mas_Rx_IQ, 0, freq_shift_Hz, SAMPLE_RATE_Hz);

% my_oscillogram(mas_Rx_IQ, SAMPLE_RATE_Hz);
my_spectr(mas_Rx_IQ, SAMPLE_RATE_Hz);
my_spectrogram_db(mas_Rx_IQ, 1024, SAMPLE_RATE_Hz);

% Оставляем целое число периодов:
mas_Rx_IQ = mas_Rx_IQ(1 : floor(length(mas_Rx_IQ) / fft_len) * fft_len);

% *** Rx: -----------------------------------------------------------------

% FIR Rx:

ripple_dB =     0.001;
suppr_dB  =   120;
f_pass_Hz = 54000;
f_stop_Hz = 60000;

[taps_FIR_Rx, nse_sigma_mul] = LPF_Remez(ripple_dB, suppr_dB, f_pass_Hz, f_stop_Hz, SAMPLE_RATE_Hz);
% taps_FIR_Rx = 1;

% save_to_float32(taps_FIR_Rx, 'data_for_cpp_maket/taps_FIR_Rx_48k_56k_120_SR192k');
% float32_to_txt(taps_FIR_Rx, 'data_for_cpp_maket/taps_FIR_Rx_48k_56k_120_SR192k');

mas_Rx_IQ_fir = conv(mas_Rx_IQ, taps_FIR_Rx, 'same');
clear mas_Rx_IQ;

% my_oscillogram(mas_Rx_IQ_fir, SAMPLE_RATE_Hz);
my_spectr(mas_Rx_IQ_fir, SAMPLE_RATE_Hz);
% my_spectrogram_db(mas_Rx_IQ_fir, 1024, SAMPLE_RATE_Hz);

% Time cyclic sync:
[sync_point, mas_resp_sync] = search_sync_point(mas_Rx_IQ_fir, mas_CSS_etalon_IQ);
if sync_point > 0
   mas_Rx_IQ_fir = mas_Rx_IQ_fir(sync_point : end);
end % if
if sync_point < 0
   mas_Rx_IQ_fir = mas_Rx_IQ_fir(fft_len+sync_point : end);
end % if

% Оставляем целое число периодов:

N_repeats = floor(length(mas_Rx_IQ_fir) / fft_len);
mas_Rx_IQ_fir = mas_Rx_IQ_fir(1 : N_repeats * fft_len);

% my_oscillogram(mas_Rx_IQ_fir, SAMPLE_RATE_Hz);
% my_spectrogram_db(mas_Rx_IQ_fir, 1024, SAMPLE_RATE_Hz);

% Считаем диаграмму "задержка-время":

matr_Delay_Time_diagram = zeros(fft_len_slow, fft_len);
for n = 1 : 1 : fft_len_slow
   mas_indexes = round((n-1)*fft_len*symbol_sync_ratio) + (1:fft_len);

   mas_resp = mas_Rx_IQ_fir(mas_indexes) .* conj(mas_CSS_etalon_IQ);
   matr_Delay_Time_diagram(n,:) = fftshift( ifft(mas_resp .* mas_window) );
end % n
clear mas_Rx_IQ_fir;

% Считаем Power-Delay Profile:
mas_PDP = mean(abs(matr_Delay_Time_diagram).^2).';

% Считаем диаграмму "задержка-Допплер":

mas_window_slow = window(@blackmanharris, fft_len_slow);

mas_indexes_slow = (1:fft_len_slow).';

matr_Delay_Doppler = zeros(fft_len_slow, fft_len);
for i_tau = 1 : 1 : fft_len
   matr_Delay_Doppler(:, i_tau) = fftshift( fft(matr_Delay_Time_diagram(mas_indexes_slow, i_tau) .* mas_window_slow) );
end % i_tau

% Оценка ОСШ:

P_all = sum(mas_PDP);
P_all_signal_fragment = sum(mas_PDP(3*fft_len/8 + (1 : 1 : fft_len/4)));

% Оценка шумовой полки в "сигнальном фрагменте" задержек:
P_nse_fragment = sum(mas_PDP(fft_len/4 + (1 : 1 : fft_len/2))) - P_all_signal_fragment;
P_nse_per_bin = P_nse_fragment / (fft_len/4);

P_signal = P_all_signal_fragment - P_nse_per_bin * fft_len/4;
if P_signal < 0
   P_signal = P_all_signal_fragment;
end % if

SNR_estim_dB = pow2db( P_signal / (P_all - P_signal) );

% Графики:

figure();
imagesc(mas_tau_ms, mas_t_sec, mag2db(abs(matr_Delay_Time_diagram)));
grid on;
aid = gca;
aid.set('CLim',[-120, -50]);
xlabel('Задержка, мс');
ylabel('Время, сек');
title('Диаграмма "Задержка-Время"');

% figure();
% mesh(mag2db(abs(matr_resp)));

figure();
plot(mas_tau_ms, pow2db( mas_PDP ), '.-');
% plot(pow2db( mas_PDP ), '.-');
grid on;
xlabel('Задержка, мс');
title('Power-Delay Profile');

figure();
imagesc(mas_tau_ms, mas_f_doppler_Hz, mag2db(abs(matr_Delay_Doppler)));
% mesh(mas_tau_ms, mas_f_doppler_Hz, mag2db(abs(matr_Delay_Doppler)));
grid on;
aid = gca;
aid.set('CLim',[-90, -30]);
xlabel('Задержка, мс');
ylabel('Допплер, Гц');
title('Диаграмма "Задержка-Допплер"');

% figure();
% plot(mas_f_doppler_Hz, mag2db(abs(matr_Delay_Doppler)).');
% grid on;
% xlabel('Допплер, Гц');
