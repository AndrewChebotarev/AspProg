function fid = my_oscillogram(matr_IQ, SAMPLE_RATE_Hz)
% Построение осциллограммы сигнала
% Физика: визуализация временной формы сигнала
% Программирование: отображение действительной, мнимой части и огибающей

   [N_IQ, N_channels] = size(matr_IQ);
   assert(N_IQ > N_channels);  % Проверка ориентации массива
   
   if (nargin > 1)
      mas_time_ms = (0:1:N_IQ-1)/SAMPLE_RATE_Hz*1000;  % Время в миллисекундах
   else
      mas_time_ms = (1:1:N_IQ);  % Без привязки к реальному времени
   end % else
   
   fid = figure();
   hold all;
   for n = 1 : 1 : N_channels
      plot(mas_time_ms, real(matr_IQ(:,n)), 'DisplayName', ['Re '  num2str(n)]);
      plot(mas_time_ms, imag(matr_IQ(:,n)), 'DisplayName', ['Im '  num2str(n)]);
      plot(mas_time_ms,  abs(matr_IQ(:,n)), 'DisplayName', ['Abs ' num2str(n)]);
   end % for n
   hold off;
   
   if (nargin > 1)
      xlabel('t, мсек');  % Подпись оси при известной частоте дискретизации
   end % if
   
   title('Осциллограмма');
   grid on;
end % function