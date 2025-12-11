function mas_IQ = read_cint32(filename, N_IQ_to_read, N_IQ_to_skip, endianness)
% Чтение комплексных чисел из файла формата cint32 (32-битные целые)
% Физика: загрузка данных с АЦП, которые часто сохраняются в целочисленном формате
% Программирование: чтение чередующихся int32 (real, imag)

   fid = fopen([filename '.cint32'], 'r');  % Открытие файла

   if (nargin < 4)
      endianness = 'ieee-le';  % Порядок байт по умолчанию (little-endian)
   end % if

   if (nargin < 3)
      N_IQ_to_skip = 0;  % Пропуск отсчетов по умолчанию
   end % if
   
   if (nargin < 2)
      % Автоматическое определение размера файла
      fseek(fid,0,'eof');
      N_IQ_to_read = ftell(fid)/4/2;  % Размер в комплексных отсчетах
   end % if
   
   if (N_IQ_to_read == 0)
      % Чтение всего файла
      fseek(fid,0,'eof');
      N_IQ_to_read = ftell(fid)/4/2;
   end % if

   % Позиционирование в файле и чтение данных
   fseek(fid, N_IQ_to_skip *4*2, 'bof');  % Пропуск указанного количества отсчетов
   mas_data = fread(fid, N_IQ_to_read*2, 'int32', 0, endianness);  % Чтение int32
   fclose(fid);

   % Формирование комплексного массива из чередующихся real, imag
   mas_IQ = complex(mas_data(1:2:end-1), mas_data(2:2:end));

end % function