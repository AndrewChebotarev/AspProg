function [masOut_IQ, masOut_IQ_noFIR] = make_css_iq(BW_sounding_Hz, ...
                                                    SAMPLE_RATE_Hz, ...
                                                    period_len, ...
                                                    reverse_len, ...
                                                    gaussfir_bt, ...
                                                    gaussfir_span, ...
                                                    gaussfir_sps)
   if nargin < 7
      gaussfir_sps = 4;
   end % if
   if nargin < 6
      gaussfir_span = 8;
   end % if
   if nargin < 5
      gaussfir_bt = 0.1;
   end % if
   if nargin < 4
      reverse_len = floor(period_len/32);
   end % if
   
   % Main sweep:
   min_freq_Hz = -BW_sounding_Hz/2;
   max_freq_Hz =  BW_sounding_Hz/2;

   dl = period_len - reverse_len;
   rl = reverse_len;
   mas_freqs_noFIR_Hz = [0           + (0:1:dl/2-1).' * max_freq_Hz/(dl/2-1); ...
                         max_freq_Hz + (0:1:rl-1).' * (min_freq_Hz-max_freq_Hz)/(rl-1); ...
                         min_freq_Hz - (0:1:dl/2-1).' * min_freq_Hz/(dl/2-1)];

   taps_FIR_freq = gaussdesign(gaussfir_bt, gaussfir_span, gaussfir_sps).';
   mas_freqs_Hz = cyclic_conv(mas_freqs_noFIR_Hz, taps_FIR_freq);

   masOut_IQ_noFIR = freqs2iq(mas_freqs_noFIR_Hz/SAMPLE_RATE_Hz);
   masOut_IQ_noFIR = fftshift(masOut_IQ_noFIR);
   
   masOut_IQ = freqs2iq(mas_freqs_Hz/SAMPLE_RATE_Hz);
   masOut_IQ = fftshift(masOut_IQ);

   ripple_dB = 0.001;
   suppr_dB = 60;
   f_pass_Hz = 1.2 * BW_sounding_Hz/2;
   f_stop_Hz = 1.4 * BW_sounding_Hz/2;
   
   if f_pass_Hz >= SAMPLE_RATE_Hz/2
      f_pass_Hz = 0.45*SAMPLE_RATE_Hz;
   end % if
   if f_stop_Hz >= SAMPLE_RATE_Hz/2
      f_stop_Hz = 0.49*SAMPLE_RATE_Hz;
   end % if
   
   taps_FIR_Tx = LPF_Remez(ripple_dB, suppr_dB, f_pass_Hz, f_stop_Hz, SAMPLE_RATE_Hz);
   
   masOut_IQ = cyclic_conv(masOut_IQ, taps_FIR_Tx);
   
end % function
