function save_to_cint32(IQ_mas, filename, extension)
% Сохранение комплексного массива в формат cint32
% Физика: сохранение данных с АЦП или для обмена с другими системами
% Программирование: запись комплексных чисел в 32-битном целочисленном формате

   if nargin < 3
      extension = 'cint32';  % Расширение файла по умолчанию
   end % if

   % Переводим в столбец:
   if (isrow(IQ_mas))
      IQ_mas = IQ_mas.';
   end % if

   % Чередование реальной и мнимой частей
   mas_tmp = [real(IQ_mas) imag(IQ_mas)].'; % чередование IQ

   fid = fopen([filename '.' extension], 'w+', 'ieee-le');  % Открытие файла
   fwrite(fid, mas_tmp(:), 'int32', 0, 'ieee-le');  % Запись 32-битных целых
   fclose(fid);
   
   clear mas_tmp;
   
end