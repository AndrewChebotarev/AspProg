
2024-10-04 Описание эксперимента:

ПРД Камаз КУНГ, Тамбовская трасса

ПРМ Камаз КУНГ, Хохольский поворот (data_rx/2024-10-04_SZV)
ПРМ Стационарная ЦАР ПРМ 8 Rx, Малышево (data_rx/2024-10-04_Mal)
ПРМ KiwiSDR по Интернету, Тамбов (data_rx/2024-10-04_Tambov)



Пояснение скриптов MatLab:

1) Формирование тестового сигнала:

index_make_wb_sounding_signal_12.m    для полосы  12 кГц
index_make_wb_sounding_signal_100.m   для полосы 100 кГц
index_make_wb_sounding_signal_300.m   для полосы 300 кГц

index_upsample_for_HackRF.m   интерполяция до sample rate 1536 кГц для SDR HackRF

IQ сигнала будет сохранено в папки:

test_iq_etalon
test_iq_tx

2) Цифровой/калибровочный эксперимент:

index_digital_experiment_12.m    для полосы  12 кГц
index_digital_experiment_100.m   для полосы 100 кГц
index_digital_experiment_300.m   для полосы 300 кГц

Можно задать известное заранее ОСШ, число лучей, задержку, допплеровские смещения лучей.

3) Обработка на приёме:

index_WinRadio_Rx_resample.m  децимация IQ c SDR WinRadio Excalibur Pro (делается 1 раз, уже выполнена)

исходные IQ берутся из папки data_rx/2024-10-04_SZV/ddc
IQ сохраняются в mat-файлы в папку data_rx/2024-10-04_SZV

index_KiwiSDR_Rx_resample.m   интерполяция IQ c KiwiSDR из Интернета (делается 1 раз, уже выполнена)

исходные IQ берутся из папки data_rx/2024-10-04_Tambov/kiwisdr
IQ сохраняются в mat-файлы в папку data_rx/2024-10-04_Tambov

Основная обработка ЛЧМ-сигнала

index_rx_processing_12.m    для полосы  12 кГц
index_rx_processing_100.m   для полосы 100 кГц
index_rx_processing_300.m   для полосы 300 кГц
