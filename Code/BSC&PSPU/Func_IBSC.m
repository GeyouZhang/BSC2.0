% This function conduct image-sequential binomial self-compensation (I-BSC) for motion error. It
% first compute four compensated fringe images, then compute motion-error-free phase.
%
% Input parameters:
% vmI      -- captured image sequence, format: H×W×(K+4), where H×W is the image resolution, K is the binomial order
% BcThresh -- pixels with modulation below this threshold will be disabled
%
% Output:
% mPhi     -- motion-error-free phase obtained through I-BSC
% mBc      -- modulation
% Author: Geyou Zhang, University of Electronic Science and Technology of China, 2025/02/27
function[ mPhi, mBc ] = Func_IBSC( vmI, BcThresh )
[h, w, iImageNum] = size(vmI);
iBinomialOrder = iImageNum - 4;
vmI_tilde = zeros(h, w, 4);
%% Compute four compensated fringe images
for i = 0:3
    for k = 0:iBinomialOrder
        V = k + 3 - mod(k + 3 - i,4);
        vmI_tilde(:, :, i + 1) = vmI_tilde(:, :, i + 1) + nchoosek(iBinomialOrder, k ) .* vmI(:, :, V + 1);
    end
end
%% Compute motion-error-free phase
mSin = vmI_tilde(:, :, 2) - vmI_tilde(:, :, 4);
mCos = vmI_tilde(:, :, 1) - vmI_tilde(:, :, 3);
mPhi = pi + atan2( - mSin, - mCos );
%% Compute motion-error-free modulation
mBc = sqrt(mSin.^2 + mCos.^2) ./ 2^(1 + iBinomialOrder);
%% Eliminate the areas that are too dark
vDark = find(mBc<BcThresh);
mPhi(vDark) = nan;
end