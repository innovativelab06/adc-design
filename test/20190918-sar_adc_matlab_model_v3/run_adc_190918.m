clear all;
clc;
%% Run Input Signal
AMP=1;
Ns=65536;
fin=0.1e5;
fs=1e6;
input=AMP*sin(2*pi*fin/fs*(1:Ns));
fft_hann(input,0,fs/2,fs,'i','input');

%% 9b sar adc - diff
[output_9b]=sar_adc_9b(Ns,fin,fs,AMP); 
real_output_9b=ideal_dac(output_9b,9);
fft_hann(real_output_9b,0,fs/2,fs,'i','output_9b');

%% ideal sar adc - diff
opt_comp.offset=0;
[output_ideal]=sar_adc(9,0,opt_comp,Ns,fin,fs,AMP,1); 
real_output_ideal=ideal_dac(output_ideal,9); 
fft_hann(real_output_ideal,0,fs/2,fs,'i','output_ideal');