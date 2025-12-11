function save_to_cint8_HackRF(IQ_mas, filename)
% Сохранение комплексного массива в формат cint8 для SDR HackRF
% Физика: подготовка сигналов для передачи через SDR-устройство HackRF
% Программирование: преобразование в 8-битный целочисленный формат, поддерживаемый HackRF

   % Переводим в столбец:
   if (isrow(IQ_mas))
      IQ_mas = IQ_mas.';
   end % if

   % Чередование реальной и мнимой частей
   mas_tmp = [real(IQ_mas) imag(IQ_mas)].'; % чередование IQ

   fid = fopen([filename '.bin'], 'w+');  % Открытие бинарного файла
   fwrite(fid, mas_tmp(:), 'int8');  % Запись 8-битных целых чисел
   fclose(fid);
   
   clear mas_tmp;
   
end % function