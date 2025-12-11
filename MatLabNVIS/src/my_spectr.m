function [figure_id, mas_spectr] = my_spectr(mas_IQ, SAMPLE_RATE_Hz, N_repeats)
% Построение спектра сигнала
% Физика: анализ частотного состава сигнала
% Программирование: вычисление и визуализация БПФ

   if (nargin <= 2)
      N_repeats = 1;  % Количество повторений для улучшения оценки спектра
   end % if

   dims = size(mas_IQ);
   if (length(dims) == 1)
      dims = [dims 1];  % Преобразование вектора в матрицу-столбец
   end
   if (dims(1) <= dims(2))
      fprintf('ERROR in my_spectr()! mas_IQ N_rows must be less than N_cols! 1 col is 1 graph!\n');
   end % if
   assert (dims(1) > dims(2));  % Проверка ориентации данных
  
   % Создание удлиненного сигнала для улучшения оценки спектра
   mas_IQ_long = [];
   for n = 1 : 1 : N_repeats
      mas_IQ_long = [mas_IQ_long; exp(2i*pi*randi(512,1,1)/512)* mas_IQ];  % Случайный фазовый сдвиг
   end % n
   
   % Определение размера БПФ
   if (2^nextpow2(dims(1)) == dims(1))
      n_fft = dims(1);  % Использовать длину сигнала если она степень 2
   else
      n_fft = 2^(nextpow2(dims(1))-1);  % Ближайшая меньшая степень 2
   end

   % Оконная функция для уменьшения эффектов спектральной утечки
   w = window(@blackmanharris, n_fft);
   
   mas_spectr = [];
   for n = 1 : 1 : dims(2)
      % Вычисление БПФ с оконной функцией
      y_mas = fft(mas_IQ_long(1:n_fft, n) .* w, n_fft);
      y_mas = 20*log10(fftshift(abs(y_mas)));  % Перевод в дБ и сдвиг нуля частоты в центр

      mas_spectr = [mas_spectr, y_mas];
   end % for n
   
   % Нормировка спектра
   max_val = max(mas_spectr(:));
   
   % Построение графика
   x_mas = (-n_fft/2 : 1 : n_fft/2-1)'/n_fft*SAMPLE_RATE_Hz/1000;  % Частота в кГц
   figure_id = figure();
   hold all;
   for n = 1 : 1 : dims(2)
      plot(x_mas, mas_spectr(:,n) - max_val);  % Нормированный спектр
   end % for n
   hold off;
   title('Спектр');
   grid on;
   xlabel('\deltaf, кГц');
   ylabel('P_н_о_р_м, дБ');  
end
