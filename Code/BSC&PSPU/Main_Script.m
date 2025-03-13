% This is the script file for motion-error-free 3D reconstruction according to the Binomial Self-compensation (BSC) technique. 
% Line 12: Switch between different sets of data
% Line 46: Switch between I-BSC and P-BSC
% Author: Geyou Zhang, University of Electronic Science and Technology of
% China, 2025/02/27
clc; clear all; close all;
%% Parameters Setting
addpath('./Package/');
load('../../Data/mCamera1Rectified.mat');
load('../../Data/mCamera2Rectified.mat');
load('../../Data/mProjector.mat');
% Select data
sFolderL = '../../Data/Statue/1_Rectified/';
sFolderR = '../../Data/Statue/2_Rectified/';
% sFolderL = '../../Data/Tissue/1_Rectified/';
% sFolderR = '../../Data/Tissue/2_Rectified/';
NSet = 4; FSet = 912/36;
iCameraWidth = 640; iCameraHeight = 480;
iFrameTotal = 50;
% Set Depth Range:
Zmin = -150; Zmax = 0;
% Set binomial order of motion error compensation, binomial order = 0 represents no compensation.
iBinomialOrder = 2;
iImageNum = iBinomialOrder + 4;
iBcThresh = 5;
% Compute Disparity Range
[ mDispMin, mDispMax ] = Func_DispartiyRange( Zmin, Zmax, mCamera1Rectified, mCamera2Rectified, iCameraHeight, iCameraWidth );

%% Image Sequence Loading
vmIL = nan( iCameraHeight, iCameraWidth, iFrameTotal );
vmIR = nan( iCameraHeight, iCameraWidth, iFrameTotal );
for i = 1:iFrameTotal  
    vmIL(:,:,i) = double( imread( sprintf( '%s%04d.bmp', sFolderL, i - 1 ) ) );
    vmIR(:,:,i) = double( imread( sprintf( '%s%04d.bmp', sFolderR, i - 1 ) ) );
end
%% Binomial Self-Compensation VS Traditional Four-step phase shifting for Dynamic 3D Scanning
% The computation speed can be accelerated by setting a larger number of workers in the MATLAB parallel pool
fig = figure;
set(gcf, 'Position', [0 0 1100 1000]); %tiledlayout(2, 2, 'TileSpacing', 'none', 'Padding', 'none');
for i = 1:iFrameTotal - iImageNum + 1
    %% 3D reconstruction with our P-BSC/I-BSC
    % Binomial Self-Compemsation for High Frequency Wrapped Phase
    tic
    % Test PBSC
%     [ mPhaseWrapLeft, mBcLeft ] = Func_PBSC( vmIL(:, :, i:i+ iImageNum - 1), iBcThresh );
%     [ mPhaseWrapRight, mBcRight ] = Func_PBSC( vmIR(:, :, i:i+ iImageNum - 1), iBcThresh );
    % Test IBSC
    [ mPhaseWrapLeft, mBcLeft ] = Func_IBSC( vmIL(:, :, i:i+ iImageNum - 1), iBcThresh );
    [ mPhaseWrapRight, mBcRight ] = Func_IBSC( vmIR(:, :, i:i+ iImageNum - 1), iBcThresh );
    dT1 = toc;
    % Correct Inherent Phase Shift
    dOffset = pi/2*mod(i - 1,4);
    mPhaseWrapLeft = mod(mPhaseWrapLeft + dOffset,2*pi);
    mPhaseWrapRight = mod(mPhaseWrapRight + dOffset,2*pi);
    
    % Paraxial Stereo Phase Unwarpping and 3D Reconstruction
    tic
    [ mX, mY, mZ, mXRaw, mYRaw, mZRaw, mPhase] = Func_Compute3D_PSPU( mPhaseWrapLeft, mPhaseWrapRight, mDispMin, mDispMax, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, 0 );
    dT2 = toc;

    % Remove the Outliers
    mInvalid = isnan( mZ ); mZ( mInvalid ) = 0; mZFilted = medfilt2( mZ, [9,9]); mOutlier = abs(mZ - mZFilted) > 2; mZ( mOutlier|mInvalid ) = nan;
    %% 3D reconstruction with traditional four-step phase shifting
    tic
    [ mPhaseWrapLeftFourStep, mBcLeftFourStep ] = Func_PBSC( vmIL(:, :, i:i + 3), iBcThresh );
    [ mPhaseWrapRightFourStep, mBcRightFourStep ] = Func_PBSC( vmIR(:, :, i:i + 3), iBcThresh );
    dT3 = toc;
    
    % Correct Inherent Phase Shift
    dOffset = pi/2*mod(i - 1,4);
    mPhaseWrapLeftFourStep = mod(mPhaseWrapLeftFourStep + dOffset,2*pi);
    mPhaseWrapRightFourStep = mod(mPhaseWrapRightFourStep + dOffset,2*pi);
    
    % Stereo Phase Unwarpping and 3D Reconstruction
    tic
    [ mXFourStep, mYFourStep, mZFourStep, mXFourStepRaw, mYFourStepRaw, mZFourStepRaw, mPhaseFourStep] = Func_Compute3D_PSPU( mPhaseWrapLeftFourStep, mPhaseWrapRightFourStep, mDispMin, mDispMax, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, 0 );
    dT4 = toc;
    disp(['Frame no.',num2str(i), '-----------------------------------------------------------------------------------------------------------']);
    disp(['BSC takes ', num2str(dT1), 's', ', PSPU and 3D reconstruction take ', num2str(dT2), 's']);
    disp(['Traditional four-step takes ', num2str(dT3), 's', ', PSPU and 3D reconstruction take ', num2str(dT4), 's']);
    
    % Remove the Outliers
    mInvalid = isnan( mZFourStep ); mZFourStep( mInvalid ) = 0; mZFilted = medfilt2( mZFourStep, [9,9]); mOutlier = abs(mZFourStep - mZFilted) > 2; mZFourStep( mOutlier|mInvalid ) = nan;
    
    % Compute the depth range for visualization
    mu = mean( mZFourStep(:), 'omitnan' ); sigma = std( mZFourStep(:), 'omitnan' ); CMin = mu - 2*sigma; CMax = mu + 2*sigma;
   
    %% Draw the point clouds
    subplot(221);
    pcshow( [mXRaw(:),mYRaw(:),mZRaw(:)] ); axis image; zlim([Zmin, Zmax]); xlim([-100, 40]); ylim([-120, 20]); 
    view([0 -90]);colormap(jet); title('BSC w/o PSPU', 'FontSize', 20, 'FontWeight', 'bold'); caxis([CMin, CMax]);     
    
    subplot(222);
    pcshow( [mX(:),mY(:),mZ(:)] ); axis image; zlim([Zmin, Zmax]); xlim([-100, 40]); ylim([-120, 20]); 
    view([0 -90]);colormap(jet); title('BSC w/ PSPU (Ours)', 'FontSize', 20, 'FontWeight', 'bold'); caxis([CMin, CMax]);    
    
    subplot(223);
    pcshow( [mXFourStepRaw(:),mYFourStepRaw(:),mZFourStepRaw(:)] ); axis image; zlim([Zmin, Zmax]); xlim([-100, 40]); ylim([-120, 20]); 
    view([0 -90]); colormap(jet);title('Traditional four-step w/o PSPU', 'FontSize', 20, 'FontWeight', 'bold'); caxis([CMin, CMax]); 
    
    subplot(224);
    pcshow( [mXFourStep(:),mYFourStep(:),mZFourStep(:)] ); axis image; zlim([Zmin, Zmax]); xlim([-100, 40]); ylim([-120, 20]); 
    view([0 -90]); colormap(jet);title('Traditional four-step w/ PSPU', 'FontSize', 20, 'FontWeight', 'bold'); caxis([CMin, CMax]); 
    
    h = axes(fig,'visible','off'); 
    cb = colorbar(h,'Position',[0.92 0.168 0.016 0.7], 'FontSize', 20, 'FontWeight', 'bold'); 
    set(get(cb,'Title'),'string','mm','FontWeight', 'bold', 'FontSize', 20, 'Color', 'white'); set(cb, 'FontWeight', 'bold', 'FontSize', 20, 'Color', 'white'); 
    caxis(h, [CMin, CMax]);
    drawnow;
end