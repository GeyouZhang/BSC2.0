% This is the script file that demonstrates the waveform generation
% procedure of OPWM technique and visualizes its mechanism via simulation
% Author: Geyou Zhang, University of Electronic Science and Technology of China, 2024/08/19 
close all; clear all;
vHarmonics = [3,5,7];
Mdesire = 0.8;
iHarmonics = length( vHarmonics );
iNotch = iHarmonics + 1;
%% The random search algorithm provides an initial value for iteration
iRandomNum = 100000;
errBest = 1000;
dM = 0.1;
while(1)
    valphaInit = pi/2*rand(iHarmonics + 1, 1);
    valphaInit = sort(valphaInit);
    err = norm( f(valphaInit, dM, vHarmonics) );
    if( err < errBest )
        errBest = err;
        valphaInitBest = valphaInit;
    end
    if( errBest < 0.3 )
        errBest = err;
        valphaInitBest = valphaInit;
        break;
    end
end
vM = dM:0.01:1;
% solve the notch position in the case of preset M
malpha = nan(length(vM),iNotch);
options=optimoptions('fsolve','Algorithm','levenberg-marquardt');
iCount = 1;
alpha = valphaInitBest;
for M = vM
    alpha = fsolve(@(a)f(a,M,vHarmonics), alpha, options);  
    malpha(iCount,:) = alpha';
    iCount = iCount + 1;
end
% generate OPWM waveform
iT = 36;
iStep = 10000;
iTtemp = iT * iStep;
vT = 1:iTtemp;
vNotch = malpha(round((Mdesire - dM)/0.01),:);
vNotch = round(vNotch./(pi/2).*(iTtemp/4-1)) + 1;
vNotch = [1,vNotch,iTtemp/4];     
iPolar = 1;
vY = nan(iTtemp/4,1);
for n = 1:length(vNotch)-1
    if(vNotch(n)==vNotch(n+1))
        iPolar = -iPolar;
        continue;
    end
    vY(vNotch(n):vNotch(n+1)) = iPolar;
    iPolar = -iPolar;
end
vY = vY(1:iStep:end);
vY = [vY;flip(vY)];
vY = [vY;-vY] + 1;
vYFFT = abs(fft(vY));
% Show the results
figure;set(gcf,'Position',[0 0 1800 400]);
subplot(131); hold on;xlabel('Preset value of M'); ylabel('Notch position');
for i = 1:iNotch
    plot(vM,malpha(:,i),'LineWidth',2);
end
subplot(132);plot(vY,'LineWidth',2);title('OPWM waveform');xlabel('Pixel'); ylabel('Intensity');
subplot(133);hold on; stem(0:18,vYFFT(1:19)); plot(vHarmonics,0,'r^','MarkerSize',6,'LineWidth',3);title('Spectrum of OPWM waveform');xlabel('Frequency'); ylabel('Amplitude');


function F = f(alpha, M, vHarmonics)
% F(1) = 1 - pi*M/4 - cos(alpha(1)) + cos(alpha(2)) - cos(alpha(3)) + cos(alpha(4));
% F(2) = 1 - cos(3*alpha(1)) + cos(3*alpha(2)) - cos(3*alpha(3)) + cos(3*alpha(4));
% F(3) = 1 - cos(5*alpha(1)) + cos(5*alpha(2)) - cos(5*alpha(3)) + cos(5*alpha(4));
% F(4) = 1 - cos(7*alpha(1)) + cos(7*alpha(2)) - cos(7*alpha(3)) + cos(7*alpha(4));
iHarmonics = length( vHarmonics );
iNotch = iHarmonics + 1;
F(1) = 1 - M;
for n = 1:iNotch
    F(1) = F(1) + 2*(-1)^n*cos(alpha(n));
end
for h = 2:iNotch
    F(h) = 1;
    for n = 1:iNotch
        F(h) = F(h) + 2*(-1)^n*cos(vHarmonics(h-1)*alpha(n));
    end
end
end

function [ Phase, Bc ] = Func_GetPhase( vvI, N )
Count = 0;
Bc = 0;
PSin = 0;
PCos = 0;
for j = 1 : N
    n = j - 1;
    vI = vvI(:,j);
    Count = Count + 1;
    PSin = PSin + vI .* sin( 2 * pi * n / N );
    PCos = PCos + vI .* cos( 2 * pi * n / N );
end
Phase = pi + atan2( - PSin, - PCos );
Bc = 2 / N * sqrt( PSin .* PSin + PCos .* PCos );
end