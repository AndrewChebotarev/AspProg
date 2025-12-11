function save_to_cfloat32(IQ_mas, filename)
% Сохранение комплексного массива в файл формата cfloat32
% Физика: сохранение эталонных сигналов и данных для передачи
% Программирование: запись комплексных чисел в бинарный формат float32

   % Переводим в столбец для единообразия:
   if (isrow(IQ_mas))
      IQ_mas = IQ_mas.';
   end % if

   % Создание массива с чередованием реальной и мнимой частей
   mas_tmp = [real(IQ_mas) imag(IQ_mas)].'; % чередование IQ

   fid = fopen([filename '.cfloat32'], 'w+', 'ieee-le');  % Открытие файла
   fwrite(fid, mas_tmp(:), 'float32', 0, 'ieee-le');  % Запись float32
   fclose(fid);
   
   clear mas_tmp;  % Очистка временной переменной
   
end % function