function [mas_out, phase_out] = my_freq_shift(mas_in, phase_in, freq_shift_Hz, SAMPLE_RATE_Hz)
% Частотный сдвиг сигнала
% Физика: компенсация допплеровского смещения частоты или настройка приемника
% Программирование: умножение на комплексную экспоненту

   N = length(mas_in);
   % Создание фазового набега для частотного сдвига
   mas_out = mas_in .* exp(1i*(2*pi*freq_shift_Hz*(0:1:N-1)/SAMPLE_RATE_Hz + phase_in)).';
   
   % Расчет накопленной фазы для следующего блока данных
   phase_out = phase_in + 2*pi*freq_shift_Hz*N/SAMPLE_RATE_Hz;
end % function