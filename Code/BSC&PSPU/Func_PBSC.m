% This function conduct phase-sequential binomial self-compensation (PBSC) for motion error. It
% first compute motion-affected phase from image sequence, then compensate
% the motion error in the phase map by BSC.
%
% Input parameters:
% vmI      -- captured image sequence, format: H×W×(K+4), where H×W is the image resolution, K is the binomial order
% BcThresh -- pixels with modulation below this threshold will be disabled
%
% Output:
% mPhi     -- motion-error-free phase obtained through P-BSC
% mBc      -- modulation
% Author: Geyou Zhang, University of Electronic Science and Technology of China, 2025/02/27
function[ mPhi, mBc ] = Func_PBSC( vmI, BcThresh )
[iCameraHeight, iCameraWidth, iFrameNum] = size(vmI);
vmPhi = zeros(iCameraHeight, iCameraWidth, iFrameNum - 3);
vmBc = zeros(iCameraHeight, iCameraWidth, iFrameNum - 3);
iBinomialOrder = iFrameNum - 4;

%% Compute motion-affected phase frames and correct the inherent phase shifting
for i = 1:iBinomialOrder + 1
    [vmPhi(:, :, i), vmBc(:, :, i)] = Func_FourStepPhaseShifting( vmI(:, :, i), vmI(:, :, i + 1), vmI(:, :, i + 2), vmI(:, :, i + 3) );
    vmPhi(:, :, i) = mod(vmPhi(:, :, i) + pi/2*( i - 1 ),2*pi);
end

%% Binomial self-compensation implented by pairwise summation layer by layer 
for k = 1:iBinomialOrder
    for s = 1:iBinomialOrder - k + 1
        vmPhi(:, :, s) = Func_AddTwoPhase( vmPhi(:, :, s), vmPhi(:, :, s + 1) );
    end
end

%% Compute modulation
mBc = 0;
for i = 1:iBinomialOrder + 1
    mBc = nchoosek(iBinomialOrder, iBinomialOrder - ( i - 1 ) ) .* vmBc(:, :, i) + mBc;
end
mBc = mBc ./ (2^iBinomialOrder);
%% Eliminate the areas that are too dark
vDark = find(mBc<BcThresh);
mPhi = mod( vmPhi(:, :, 1), 2*pi);
mPhi(vDark) = nan;
end

function[ mPhi ] = Func_AddTwoPhase( mPhi0, mPhi1 )
    mLeap0 = abs(mPhi1 - mPhi0)>pi;
    mPhi = ( mPhi0 + mPhi1 + 2*pi.*mLeap0 ) ./ 2;
end