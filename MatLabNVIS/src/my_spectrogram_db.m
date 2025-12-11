function figure_id = my_spectrogram_db(IQ_mas, N_freqs, SAMPLE_RATE_Hz)
% Построение спектрограммы сигнала в децибелах
% Физика: визуализация изменения спектра во времени (время-частота-интенсивность)
% Программирование: использование встроенной функции spectrogram с дополнительной обработкой

   % Вычисление спектрограммы с оконной функцией Блэкмана-Харриса
   SpecGram = spectrogram(IQ_mas, blackman(N_freqs), N_freqs/2, N_freqs);

   % Определение типа сигнала (вещественный или комплексный)
   if (sum(abs(real(IQ_mas))) == sum(abs(IQ_mas)))
      % Вещественный сигнал - только положительные частоты
      f_mas = (0:N_freqs/2)/N_freqs*SAMPLE_RATE_Hz/1000;  % Частота в кГц
   else
      % Комплексный сигнал - полный диапазон частот
      SpecGram = fftshift(SpecGram, 1);  % Сдвиг нуля частоты в центр
      f_mas = ((0:N_freqs)/N_freqs - 0.5)*SAMPLE_RATE_Hz/1000;  % ± частоты в кГц
   end

   dim = length(SpecGram);  % Количество временных отсчетов в спектрограмме
   % Временная ось в миллисекундах
   t_mas = (0:dim)/dim * length(IQ_mas)/SAMPLE_RATE_Hz*1000;

   % Нормировка спектрограммы (эмпирический коэффициент)
   SpecGram = SpecGram/106.73;

   % Перевод амплитуды в децибелы
   SpecGram_dB = 20*log10(abs(SpecGram));

   % Построение изображения спектрограммы
   figure_id = figure();
   imagesc(t_mas, f_mas, SpecGram_dB);
   % Создание градиентной цветовой карты (от белого к черному)
   map = (1 : -0.0033 : 0).^0.5';
   colormap([map map map]);  % Серошкальная карта
   grid on;
   xlabel('t, мсек');
   ylabel('f, кГц');

end % function