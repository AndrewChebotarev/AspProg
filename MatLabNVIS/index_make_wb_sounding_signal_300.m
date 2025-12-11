
% Формирование ЛЧМ-сигнала для зондирования ДКМВ-канала в полосе 300 кГц.

% Причины:
% 1) Ширина полосы ограничена 300 кГц.
% 2) Хотим создать сигнал, с помощью которого можно увидеть лучи с разносом до 9 мс,
%    а также Doppler shift/spread от -5 до +5 Гц.
% 3) При SAMPLE_RATE 384 кГц мы отфильтруем на приёме 300 кГц (от -150 до +150 кГц).
%    Из-за этого рабочей областью Power_Delay_Profile (p.d.p.) будет половина выхода ifft.
%    При длине ifft 16384 мы получим длительность символа ЛЧМ 42.(6) мс,
%    а также рабочую часть p.d.p. от -5.(3) до +5.(3) мс. Это 16384 точек. Это нас устраивает.
% 4) Набрав 4096 ЛЧМ-символов ( 3 минуты ), мы сможем сделать fft-4096
%    по оси "медленного" времени, и получить диаграмму "задержка-Допплер".
%    Тогда у нас получится диапазон Допплера от -11.71875 до +11.71875 Гц.
%    Это 4096 точек. Это нас устраивает.
%
% 11.71875 Гц = 1/2/(42.(6) мс)

clc;
clear;

SAMPLE_RATE_Hz = 384000;

fft_loglen = 14;
fft_len = 2^fft_loglen;

fprintf('period %5.3f ms (%u samples)\n', fft_len/SAMPLE_RATE_Hz*1000, fft_len);

% Параметры зондирующего ЛЧМ-сигнала:

chirp_bandwidth_Hz = 300000;
reverse_chirp_divider = 32; % доля времени на "обратный" chirp
att_dB = 0.25;

att_mul = 10^(-att_dB/20);

mas_CSS_300k_IQ = att_mul * make_css_iq(chirp_bandwidth_Hz, SAMPLE_RATE_Hz, fft_len, fft_len/reverse_chirp_divider );

save_to_cfloat32(mas_CSS_300k_IQ, 'test_iq_etalon/etalon_300kHz');

% cfloat32_to_txt(mas_CSS_300k_IQ, 'data_for_cpp_maket/etalon_300kHz');

%%

N_repeats = 10;

mas_Tx_300k_IQ = mas_CSS_300k_IQ * ones(1,N_repeats);
mas_Tx_300k_IQ = mas_Tx_300k_IQ(:);

my_oscillogram(mas_Tx_300k_IQ, SAMPLE_RATE_Hz);
my_spectr(mas_Tx_300k_IQ, SAMPLE_RATE_Hz);
my_spectrogram_db(mas_Tx_300k_IQ, 1024, SAMPLE_RATE_Hz);

%%

save_to_cfloat32(mas_Tx_300k_IQ, 'test_iq_tx/sounding_300kHz');
