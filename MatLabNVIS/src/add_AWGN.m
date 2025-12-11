function mas_out = add_AWGN(mas_in, mas_nse_sigma)
% Добавление аддитивного белого гауссовского шума (AWGN) к сигналу
% Физика: моделирует тепловой шум в приемнике и атмосферные помехи
% Программирование: создает комплексный шум с заданным СКО

   [N, N_ch] = size(mas_in);
   assert(N>N_ch);  % Проверка: строк должно быть больше чем столбцов
   
   N_ch_nse = length(mas_nse_sigma);
   if (N_ch_nse>1)
      assert(N_ch_nse == N_ch);  % Проверка соответствия размерностей
   else
      mas_nse_sigma = mas_nse_sigma * ones(N_ch,1);  % Расширение скаляра до вектора
   end % else
   
   mas_out = 0 * mas_in; % Выделение памяти для выходного массива
   for i_ch = 1 : 1 : N_ch
      % Генерация комплексного гауссовского шума
      % Физика: моделирует шум в квадратурных каналах приемника
      mas_out(:,i_ch) = mas_in(:,i_ch) + 0.707 * mas_nse_sigma(i_ch) * (randn(N,1) + 1i * randn(N,1));
   end % i_ch
end % function