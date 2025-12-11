function [sync_point, mas_resp_sync] = search_sync_point(masIn_Rx_IQ, masIn_etalon)

   N = length(masIn_etalon);

   taps_etalon = flip(conj(masIn_etalon));
   mas_resp = conv(masIn_Rx_IQ, taps_etalon, 'same');

   N_repeats = floor(length(masIn_Rx_IQ)/N);
   
   matr_resp = reshape(mas_resp(1:N*N_repeats), N, N_repeats);

   mas_resp_sync = mean( abs(matr_resp.').^2 );

   [vmax, imax] = max(mas_resp_sync);
%   sync_point = imax - N/2;
   
   % Ищем sync_point как центр между первым и последним лучом.
   % Всё, что выше -12 дБ от max считаем лучом.
   
   sync_lo = 0;
   sync_hi = 0;
   for i = 1 : 1 : length(mas_resp_sync)
      
      if (mas_resp_sync(i) > vmax/4/4)
         if (sync_lo == 0)
            sync_lo = i;
         end % if
         
         sync_hi = i;
      end % if
   end % i
   
   sync_point = round((sync_lo + sync_hi)/2) - N/2;
   
end % function
