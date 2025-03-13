% This is the script file for calibrating camera nonlinear response curve
% Author: Geyou Zhang, University of Electronic Science and Technology of China, 2024/08/19 
clc; clear all; close all;
N = 80;
vIc = zeros( N, 1 ); 
vE=[ 299.0620 355.1360 411.2100 467.2840 523.3580 579.4320 635.5060 672.8890 728.9630 785.0370 841.1110 897.1850 953.2590 1009.3330 1065.4070 1121.4810 1177.5560 1233.6300 1271.0120 1327.0860 1383.1600 1439.2350 1495.3090 1551.3830 1607.4570 1663.5310 1719.6050 1775.6790 1831.7530 1869.1360 1925.2100 1981.2840 2037.3580 2093.4320 2149.5060 2205.5800 2261.6540 2317.7280 2373.8020 2429.8770 2485.9510 2523.3330 2579.4070 2635.4810 2691.5560 2747.6300 2803.7040 2859.7780 2915.8520 2971.9260 3028.0000 3084.0740 3121.4570 3177.5310 3233.6050 3289.6790 3345.7530 3401.8270 3457.9010 3513.9750 3570.0490 3626.1230 3682.1980 3738.2720 3775.6540 3831.7280 3887.8020 3943.8770 3999.9510 4056.0250 4112.0990 4168.1730 4224.2470 4280.3210 4336.3950 4373.7780 4429.8520 4485.9260 4542.0000 4598.0740 ];
vE = vE';
% vE = (1:N)';
for j = 1:N
    n = j - 1;
    Im = double( imread( sprintf( './CameraResponseCalibration/%04d.bmp', n ) ) ); 
%     Im = Im(230,310);
%     Im( 255 == Im ) = nan;
    vIc( j ) = mean( Im( : ), 'omitnan' );
end

% Set up fittype and options.
ft = 'linearinterp';
[Interp, gof] = fit( vIc, vE, ft, 'Normalize', 'on' );
vE = [0;vE;Interp(255)];
vIc = [0;vIc;255];
figure;grid on;hold on;
plot(vE, vIc,'r*');
vE = 255.*(vE - min(vE)) / (max(vE) - min(vE));
% vIc = 255.*(vIc - min(vIc)) / (max(vIc) - min(vIc));
[Interp, gof] = fit( vIc, vE, ft, 'Normalize', 'on' );


FontSize = 20;
figure;grid on;hold on;set(gcf, 'Position',[0  0  900  500]);xlim([0 255]);ylim([0 255]);
plot( vIc, vE,'.-', 'LineWidth', 2,'MarkerSize',15, 'Color', [1,0.5,0] );
plot(1:255,1:255, 'b-', 'LineWidth', 2, 'Color', [0,0.5,1] );
set( gca, 'FontSize', FontSize, 'FontWeight', 'bold' );
xlabel('Image Intensity');
ylabel('Normalized Irradiance');
legend('Camera Response Function','Reference Linear Response');
save Interp Interp


