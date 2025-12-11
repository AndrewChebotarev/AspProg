% Очистка командного окна MATLAB
clc;
% Очистка рабочей области (удаление всех переменных)
clear;

% Добавление папки 'src' в путь поиска функций MATLAB
addpath('src');

% Установка частоты дискретизации 48 кГц для полосы 12 кГц
SAMPLE_RATE_Hz = 48000;

% Параметры зондирующего ЛЧМ-сигнала:

% Ширина полосы ЛЧМ-сигнала 12 кГц
chirp_bandwidth_Hz = 12000;
% Коэффициент для обратного chirp (1/32 времени)
reverse_chirp_divider = 32; % доля времени на "обратный" chirp

% Обработка без децимации (без уменьшения частоты дискретизации)

% Загрузка эталонного ЛЧМ-сигнала из файла формата cfloat32
mas_CSS_etalon_IQ = read_cfloat32('test_iq_etalon/etalon_12kHz');

% Определение длины ЛЧМ-сигнала (2048 отсчетов)
fft_len = length(mas_CSS_etalon_IQ);

% Количество "медленных" отсчетов для анализа Допплера
fft_len_slow = 4096;

% Расчет диапазона задержек в миллисекундах
% Формула: полоса * коэффициент / (коэффициент - 1)
tau_range_ms = chirp_bandwidth_Hz * reverse_chirp_divider / (reverse_chirp_divider-1);

% Создание масштабных осей для графиков:
% Ось задержек в миллисекундах (-N/2+1 до N/2 точек)
mas_tau_ms = (-fft_len/2+1 : 1 : fft_len/2) / tau_range_ms * 1000;
% Ось времени в секундах для диаграммы "задержка-время"
mas_t_sec = (1:1:fft_len_slow) * fft_len / SAMPLE_RATE_Hz;
% Ось допплеровских частот в Герцах для диаграммы "задержка-Допплер"
mas_f_doppler_Hz = ((1:fft_len_slow).' - fft_len_slow/2 - 1) * SAMPLE_RATE_Hz / fft_len / fft_len_slow;

% Выбор оконной функции для обработки сигнала:
% mas_window = ones(fft_len, 1);                              % Прямоугольное окно (не используется)
% mas_window = window(@blackmanharris, fft_len);              % Окно Блэкмана-Харриса (не используется)
% Кастомное окно с нулевыми зонами по краям для уменьшения краевых эффектов
mas_window = [zeros(fft_len/8,1); window(@blackmanharris, 3*fft_len/4); zeros(fft_len/8,1)];


% *** Цифровой эксперимент (моделирование канала):

% Загрузка данных цифрового эксперимента (закомментировано)
% mas_Rx_IQ = read_cfloat32('data_rx\digital\css12_rx');


% *** Хохольский поворот, WinRadio (разные файлы экспериментальных данных):

% load('data_rx\2024-10-04_SZV\2024-10-04_15-00_LCHM_12k_11960_sr48kHz.mat');             % сигнал есть
% symbol_sync_ratio = fft_len / (fft_len + (1.45+0.89)*SAMPLE_RATE_Hz/1000/fft_len_slow);
% freq_shift_Hz = 0;

% mas_Rx_IQ = read_cfloat32('data_rx\2024-10-04_SZV\LCHM_12k_2500_sr48kHz'); % nothing
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% mas_Rx_IQ = read_cfloat32('data_rx\2024-10-04_SZV\LCHM_12k_3100_sr48kHz'); % nothing
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% mas_Rx_IQ = read_cfloat32('data_rx\2024-10-04_SZV\LCHM_12k_3540_sr48kHz'); % nothing
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% mas_Rx_IQ = read_cfloat32('data_rx\2024-10-04_SZV\LCHM_12k_4490_sr48kHz'); % брызги
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% mas_Rx_IQ = read_cfloat32('data_rx\2024-10-04_SZV\LCHM_12k_5075_1_sr48kHz'); % сигнал есть, но разрывы на приёме
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% mas_Rx_IQ = read_cfloat32('data_rx\2024-10-04_SZV\LCHM_12k_5075_sr48kHz'); % сигнал есть, но разрывы на приёме
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% mas_Rx_IQ = read_cfloat32('data_rx\2024-10-04_SZV\LCHM_12k_6500_sr48kHz'); % сигнал есть, но разрывы на приёме
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% mas_Rx_IQ = read_cfloat32('data_rx\2024-10-04_SZV\LCHM_12k_6535_sr48kHz'); % сигнал есть, но разрывы на приёме
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% mas_Rx_IQ = read_cfloat32('data_rx\2024-10-04_SZV\LCHM_12k_6540_sr48kHz'); % сигнал есть, но разрывы на приёме
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% mas_Rx_IQ = read_cfloat32('data_rx\2024-10-04_SZV\LCHM_11960_sr48kHz');
% symbol_sync_ratio = fft_len / (fft_len + (1.45+0.89)*SAMPLE_RATE_Hz/1000/fft_len_slow); % сигнал есть
% freq_shift_Hz = 0;


% *** Малышево, 8 Rx (активный вариант загрузки данных):

% Загрузка данных эксперимента из Малышево с 8 приемниками
load('data_rx\2024-10-04_Mal\2024-10-04_21-33_LCHM_12k_sr48kHz.mat'); % !!! сигнал есть
% Коэффициент синхронизации символов
symbol_sync_ratio = 1;
% Частотный сдвиг, подобранный по лучу земной волны для компенсации допплера
freq_shift_Hz = 10.322; % подобрано по лучу земной волны

% Выбор второго приемного канала из массива данных
i_Rx = 2;
% Извлечение данных выбранного приемного канала
mas_Rx_IQ = matr_Rx_IQ(:,i_Rx);
% Очистка неиспользуемой переменной из памяти
clear matr_Rx_IQ_rec;


% *** Тамбов, KiwiSDR (закомментированные варианты загрузки данных):

% load('data_rx\2024-10-04_Tambov\2024-10-04_14-25_LCHM_12k_12000_sr48kHz.mat'); % сигнал есть
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;
% fft_len_slow = 128;

% load('data_rx\2024-10-04_Tambov\2024-10-04_15-02_LCHM_12k_11960_sr48kHz.mat'); % сигнал есть, нужен фильтр помехи
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% load('data_rx\2024-10-04_Tambov\2024-10-04_15-13_LCHM_12k_11960_sr48kHz.mat'); % сигнал есть, нужен фильтр помехи
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;
% fft_len_slow = 128;

% load('data_rx\2024-10-04_Tambov\2024-10-04_19-54_LCHM_12k_3540_sr48kHz.mat'); % сигнал есть, нужен фильтр помехи
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;
% fft_len_slow = 512;

% load('data_rx\2024-10-04_Tambov\2024-10-04_20-42_LCHM_12k_6540_sr48kHz.mat'); % сигнала нет, нужен фильтр помехи
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% load('data_rx\2024-10-04_Tambov\2024-10-04_20-48_LCHM_12k_6540_sr48kHz.mat'); % сигнала нет, нужен фильтр помехи
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;

% load('data_rx\2024-10-04_Tambov\2024-10-04_20-53_LCHM_12k_6540_sr48kHz.mat'); % сигнала нет, нужен фильтр помехи
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;
% fft_len_slow = 128;

% load('data_rx\2024-10-04_Tambov\2024-10-04_20-55_LCHM_12k_6540_sr48kHz.mat'); % сигнала нет, нужен фильтр помехи
% symbol_sync_ratio = 1;
% freq_shift_Hz = 0;
% fft_len_slow = 128;

% load('data_rx\2024-10-04_Tambov\2024-10-04_20-56_LCHM_12k_6540_sr48kHz.mat'); % !!! сигнал есть
% symbol_sync_ratio = fft_len / (fft_len + (16)*SAMPLE_RATE_Hz/1000/fft_len_slow); % 16 мс подобрано вручную
% freq_shift_Hz = 8; % подобрано вручную


% Применение частотного сдвига к принятому сигналу для компенсации допплера
mas_Rx_IQ = my_freq_shift(mas_Rx_IQ, 0, freq_shift_Hz, SAMPLE_RATE_Hz);

% Визуализация сигнала на различных этапах обработки:
% Построение осциллограммы принятого сигнала
my_oscillogram(mas_Rx_IQ, SAMPLE_RATE_Hz);
% Построение спектра принятого сигнала
my_spectr(mas_Rx_IQ, SAMPLE_RATE_Hz);
% Построение спектрограммы первых 180 секунд сигнала
my_spectrogram_db(mas_Rx_IQ(1:180*SAMPLE_RATE_Hz), 1024, SAMPLE_RATE_Hz);

%%

% Обрезка сигнала до целого числа периодов ЛЧМ-сигнала:
% Это необходимо для корректной обработки периодического сигнала
mas_Rx_IQ = mas_Rx_IQ(1 : floor(length(mas_Rx_IQ) / fft_len) * fft_len);

% *** Приемная обработка (Rx): ----------------------------------------------------------------

% Параметры КИХ-фильтра нижних частот для подавления внеполосных помех:

% Неравномерность в полосе пропускания [дБ]
ripple_dB = 0.001;
% Подавление в полосе задерживания [дБ]
suppr_dB  = 120;
% Граничная частота полосы пропускания [Гц]
f_pass_Hz = 6000;
% Граничная частота полосы задерживания [Гц]
f_stop_Hz = 8000;

% Синтез КИХ-фильтра методом Ремеза с заданными параметрами
[taps_FIR_Rx, nse_sigma_mul] = LPF_Remez(ripple_dB, suppr_dB, f_pass_Hz, f_stop_Hz, SAMPLE_RATE_Hz);
% taps_FIR_Rx = 1;  % Альтернатива: без фильтрации (закомментировано)

% Сохранение коэффициентов фильтра для использования в других приложениях:
% save_to_float32(taps_FIR_Rx, 'data_for_cpp_maket/taps_FIR_Rx_24k_28k_120_SR96k');
% float32_to_txt(taps_FIR_Rx, 'data_for_cpp_maket/taps_FIR_Rx_6k_8k_120_SR48k');

% Фильтрация принятого сигнала КИХ-фильтром
% conv с параметром 'same' сохраняет исходную длину сигнала
mas_Rx_IQ_fir = conv(mas_Rx_IQ, taps_FIR_Rx, 'same');
% Очистка исходного сигнала для экономии памяти
clear mas_Rx_IQ;

% Визуализация отфильтрованного сигнала:
% my_oscillogram(mas_Rx_IQ_fir, SAMPLE_RATE_Hz);  % Осциллограмма (закомментировано)
% Построение спектра отфильтрованного сигнала
my_spectr(mas_Rx_IQ_fir, SAMPLE_RATE_Hz);
% my_spectrogram_db(mas_Rx_IQ_fir, 1024, SAMPLE_RATE_Hz);  % Спектрограмма (закомментировано)

% Временная синхронизация по циклической корреляции:
% Поиск точки начала сигнала путем корреляции с эталоном
[sync_point, mas_resp_sync] = search_sync_point(mas_Rx_IQ_fir, mas_CSS_etalon_IQ);
% Обрезка сигнала в зависимости от найденной точки синхронизации
if sync_point > 0
   % Если точка синхронизации положительная - обрезаем начало
   mas_Rx_IQ_fir = mas_Rx_IQ_fir(sync_point : end);
end % if
if sync_point < 0
   % Если точка синхронизации отрицательная - обрезаем с конца
   mas_Rx_IQ_fir = mas_Rx_IQ_fir(fft_len+sync_point : end);
end % if

% Обрезка до целого числа периодов после синхронизации:
% Расчет количества полных периодов ЛЧМ в сигнале
N_repeats = floor(length(mas_Rx_IQ_fir) / fft_len);
% Обрезка до целого числа периодов
mas_Rx_IQ_fir = mas_Rx_IQ_fir(1 : N_repeats * fft_len);

% Дополнительная визуализация синхронизированного сигнала:
% my_oscillogram(mas_Rx_IQ_fir, SAMPLE_RATE_Hz);        % Осциллограмма (закомментировано)
% my_spectrogram_db(mas_Rx_IQ_fir, 1024, SAMPLE_RATE_Hz);  % Спектрограмма (закомментировано)

% Построение диаграммы "задержка-время":
% Инициализация матрицы для хранения результатов
matr_Delay_Time_diagram = zeros(fft_len_slow, fft_len);
% Обработка каждого "медленного" отсчета
for n = 1 : 1 : fft_len_slow
   % Выбор сегмента сигнала для обработки с учетом синхронизации
   mas_indexes = round((n-1)*fft_len*symbol_sync_ratio) + (1:fft_len);

   % Корреляция выбранного сегмента с эталонным сигналом
   mas_resp = mas_Rx_IQ_fir(mas_indexes) .* conj(mas_CSS_etalon_IQ);
   % Обратное БПФ для получения импульсной характеристики канала
   % fftshift сдвигает нулевую задержку в центр
   matr_Delay_Time_diagram(n,:) = fftshift( ifft(mas_resp .* mas_window) );
end % n
% Очистка отфильтрованного сигнала для экономии памяти
clear mas_Rx_IQ_fir;

% Расчет Power-Delay Profile (PDP):
% Усреднение мощности по времени для получения профиля задержек
mas_PDP = mean(abs(matr_Delay_Time_diagram).^2).';

% Построение диаграммы "задержка-Допплер":

% Выбор оконной функции для допплер-анализа:
% mas_window_slow = window(@blackmanharris, fft_len_slow);  % Окно Блэкмана-Харриса (закомментировано)
% Окно Тейлора для анализа Допплера (лучшее подавление боковых лепестков)
mas_window_slow = window(@taylorwin, fft_len_slow);

% Индексы для "медленного" времени (оси Допплера)
mas_indexes_slow = (1:fft_len_slow).';

% Инициализация матрицы для диаграммы "задержка-Допплер"
matr_Delay_Doppler = zeros(fft_len_slow, fft_len);
% Обработка каждой задержки для анализа Допплера
for i_tau = 1 : 1 : fft_len
   % БПФ по "медленному" времени для анализа допплеровских сдвигов
   % fftshift сдвигает нулевую частоту Допплера в центр
   matr_Delay_Doppler(:, i_tau) = fftshift( fft(matr_Delay_Time_diagram(mas_indexes_slow, i_tau) .* mas_window_slow) );
end % i_tau

% Оценка отношения сигнал-шум (ОСШ):

% Общая мощность в профиле задержек
P_all = sum(mas_PDP);
% Мощность в сигнальной области (центральная часть задержек)
P_all_signal_fragment = sum(mas_PDP(3*fft_len/8 + (1 : 1 : fft_len/4)));

% Оценка шумовой полки в области задержек вне сигнала:
% Мощность в области, где ожидается только шум
P_nse_fragment = sum(mas_PDP(fft_len/4 + (1 : 1 : fft_len/2))) - P_all_signal_fragment;
% Мощность шума на один бин задержки
P_nse_per_bin = P_nse_fragment / (fft_len/4);

% Расчет мощности сигнала путем вычитания шумовой составляющей
P_signal = P_all_signal_fragment - P_nse_per_bin * fft_len/4;
% Защита от отрицательных значений мощности сигнала
if P_signal < 0
   P_signal = P_all_signal_fragment;
end % if

% Расчет ОСШ в децибелах
SNR_estim_dB = pow2db( P_signal / (P_all - P_signal) );

% Построение графиков результатов:

% График 1: Диаграмма "задержка-время"
figure(1);
% Отображение диаграммы в виде изображения с цветовой кодировкой амплитуды
imagesc(mas_tau_ms, mas_t_sec, mag2db(abs(matr_Delay_Time_diagram)));
% Включение сетки на графике
grid on;
% Получение handle текущих осей
aid = gca;
% Установка динамического диапазона цветовой шкалы от -120 до -35 дБ
aid.set('CLim',[-120, -35]);
% Подпись оси X
xlabel('Задержка, мс');
% Подпись оси Y
ylabel('Время, сек');
% Заголовок графика
title('Диаграмма "Задержка-Время"');

% Альтернативная 3D визуализация (закомментировано):
% figure();
% mesh(mag2db(abs(matr_resp)));

% График 2: Power-Delay Profile (PDP)
figure(2);
% Построение графика профиля задержек в логарифмическом масштабе
plot(mas_tau_ms, pow2db( mas_PDP ), '.-');
% Альтернативный вариант без масштабирования по задержкам (закомментировано):
% plot(pow2db( mas_PDP ), '.-');
% Включение сетки
grid on;
% Подпись оси X
xlabel('Задержка, мс');
% Заголовок графика
title('Power-Delay Profile');

% График 3: Диаграмма "задержка-Допплер"
figure(3);
% Отображение диаграммы в виде изображения с цветовой кодировкой
imagesc(mas_tau_ms, mas_f_doppler_Hz, mag2db(abs(matr_Delay_Doppler)));
% Альтернативная 3D визуализация (закомментировано):
% mesh(mas_tau_ms, mas_f_doppler_Hz, mag2db(abs(matr_Delay_Doppler)));
% Включение сетки
grid on;
% Получение handle текущих осей
aid = gca;
% Установка динамического диапазона цветовой шкалы от -90 до -15 дБ
aid.set('CLim',[-90, -15]);
% Подпись оси X
xlabel('Задержка, мс');
% Подпись оси Y
ylabel('Допплер, Гц');
% Заголовок графика
title('Диаграмма "Задержка-Допплер"');

% Альтернативное представление диаграммы Допплера (закомментировано):
% figure();
% plot(mas_f_doppler_Hz, mag2db(abs(matr_Delay_Doppler)).');
% grid on;
% xlabel('Допплер, Гц');