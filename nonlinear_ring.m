%% nonlinear ring spectrum fitting
close all
%% Experimental data
inspect_lambda_min = 1543;
inspect_lambda_max = 1546;

[wl,GC,zero_bias, minus_half_bias, minus_one_bias, minus_one_half_bias,  minus_two_bias] = import_excel('ring spectrum.xlsx','0dbm',inspect_lambda_min,inspect_lambda_max);

ploting(wl,{-GC},'GC',14)

ploting(wl,{ -zero_bias-GC },'ring',14)

[wl,GC,zero_bias, minus_half_bias, minus_one_bias, minus_one_half_bias,  minus_two_bias] = import_excel('ring spectrum.xlsx','10dbm',inspect_lambda_min,inspect_lambda_max);

ploting(wl,{ -zero_bias-GC },'ring',14)