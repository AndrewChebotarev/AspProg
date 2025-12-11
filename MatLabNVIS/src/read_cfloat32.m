function mas_IQ = read_cfloat32(filename, N_IQ_to_read, N_IQ_to_skip)
% Чтение комплексных чисел из файла формата cfloat32
% Физика: загрузка записанных сигналов с приемников
% Программирование: чтение чередующихся float32 (real, imag)

   fid = fopen([filename '.cfloat32'], 'r');  % Открытие файла

   if (nargin < 3)
      N_IQ_to_skip = 0;  % Пропуск отсчетов по умолчанию
   end % if
   
   if (nargin < 2)
      fseek(fid,0,'eof');
      N_IQ_to_read = ftell(fid)/4/2;  % Автоматическое определение длины файла
   end % if
   
   if (N_IQ_to_read == 0)
      fseek(fid,0,'eof');
      N_IQ_to_read = ftell(fid)/4/2;  % Чтение всего файла
   end % if

   % Позиционирование и чтение данных
   fseek(fid, N_IQ_to_skip *4*2, 'bof');
   mas_data = fread(fid, N_IQ_to_read*2, 'float32', 0, 'ieee-le');  % Чтение float32
   fclose(fid);

   % Формирование комплексного массива из чередующихся real, imag
   mas_IQ = complex(mas_data(1:2:end-1), mas_data(2:2:end));
end % function 