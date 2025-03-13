% This function computes the 3D point clouds from the wrapped phase map of
% both left and right cameras. This function corresponds to PSPU in section 5.
% Author: Geyou Zhang, University of Electronic Science and Technology of China, 2024/09/05
%
% Input parameters:
% mPhaseWrapLeft       -- wrapped phase map of left camera
% mPhaseWrapRight      -- wrapped phase map of right camera
% mDispMin, mDispMax   -- disparity range in auxiliary camera corresponding to each main camera pixel (x^m, y^m)
% mCamera1Rectified    -- perspective projection matrix of the main camera
% mCamera2Rectified    -- perspective projection matrix of the auxiliary camera
% mProjector           -- perspective projection matrix of the projector
% FSet                 -- frequency of fringe pattern
% bBlockMatchingFlag   -- equals to 1 indicates using block matching algorithm
%                         equals to 0 indicates using single-pixel matching algorithm (recommended)
%
% Output:
% mX,mY,mZ             -- reconstructed point cloud with PSPU
% mXRaw,mYRaw,mZRaw    -- reconstructed point cloud without PSPU
% mPhase               -- unwrapped phase map of left camera
function[ mX, mY, mZ, mXRaw, mYRaw, mZRaw, mPhase ] = Func_Compute3D_PSPU( mPhaseWrapLeft, mPhaseWrapRight, mDispMin, mDispMax, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, bBlockMatchingFlag )
% conduct stereo phase unwrapping to have a noisy absolute phase map
[ mPhase ] = Func_StereoPhaseUnwrapping( mPhaseWrapLeft, mPhaseWrapRight, mDispMin, mDispMax, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, bBlockMatchingFlag );
mXpUnwrapRaw = mPhase ./ (2*pi) .* 1280;
% compute raw 3D point cloud without outlier rectification
[ mXRaw, mYRaw, mZRaw ] = Func_Compute3DXp( mCamera1Rectified, mProjector, mXpUnwrapRaw, 0, 0 );

% conduct hierarchical adjustment for phase outlier
iThreshReliable = 3000;
[ mPhase ] = Func_RectifyOutlierPhase( mPhase, FSet, 3, iThreshReliable );
mXpUnwrap = mPhase ./ (2*pi) .* 1280;
% compute 3D point cloud
[ mX, mY, mZ ] = Func_Compute3DXp( mCamera1Rectified, mProjector, mXpUnwrap, 0, 0 );
end

% This function conducts stereo phase unwrapping to have a noisy absolute phase map of left camera
function[ mPhase ] = Func_StereoPhaseUnwrapping( mPhaseWrapLeft, mPhaseWrapRight, mDispMin, mDispMax, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, bBlockMatchingFlag )
if( bBlockMatchingFlag )
    [ mXTarget ] = Func_PhaseMatchSAD( mPhaseWrapLeft, mPhaseWrapRight, mDispMin, mDispMax,3);
else
    [ mXTarget ] = Func_PhaseMatchPoint( mPhaseWrapLeft, mPhaseWrapRight, mDispMin, mDispMax);
end
[ mXMatch, mYMatch, mZMatch ] = Func_Compute3DXp( mCamera1Rectified, mCamera2Rectified, mXTarget, 0, 0 );
vP = mProjector * [ mXMatch(:), mYMatch(:), mZMatch(:), ones( numel( mZMatch(:) ), 1 ) ]';
vP = vP ./ vP( 3, : );
vxp = vP( 1, : );
mXp = reshape( vxp, size(mPhaseWrapLeft) );
mPhiReference = mXp ./ 1280 .*(2*pi);
mPhaseOrder = round( ( FSet(end) .* mPhiReference - mPhaseWrapLeft ) ./ (2*pi) );
mPhase = ( mPhaseOrder .* (2*pi) + mPhaseWrapLeft ) ./ FSet(end);
end

% This function employs hierarchical adjustment strategy to rectify the outliers in the noisy absolute phase map
%
% Input parameters:
% mPhaseInput          -- noisy absolute phase map
% FSet                 -- frequency of fringe pattern
% iteNum               -- iterations
% iThreshReliable      -- area threshold of reliable connected domain, 
% connected domains that exceed the threshold size will be marked as 
% reliable connected domains
%
% Output:
% mPhase               -- absolute phase map after hierarchical adjustment
function [ mPhase ] = Func_RectifyOutlierPhase( mPhaseInput, FSet, iteNum, iThreshReliable )
mPhaseInputOrg = mPhaseInput;   
[H, W] = size(mPhaseInput);
mDiffer = 6.28.*ones(H, W);
for k = 1:iteNum
    mLabel = zeros(H, W);
    mLabel( isnan(mPhaseInput) ) = -1;
    % divide the connected domains by phase continuity
    [ mLabel ] = Func_PhaseFloodFillSegmentation( mLabel, mPhaseInput, FSet );
    % relabel connected domains in descending area order
    stats = regionprops(mLabel,'Area');
    areas = [stats.Area];
    [~, sortedIdx] = sort(areas, 'descend');
    temp = -1 .* ones(H, W);
    for i = 1:length(sortedIdx)
        temp(mLabel == sortedIdx(i)) = i;
    end
    mLabel = temp;
    stats = regionprops(mLabel,'Area','PixelIdxList');
    areas = [stats.Area];
    % connected domains that exceed the threshold size will be marked as reliable connected domains
    for i = 1:length(areas)
        if(areas(i)<iThreshReliable)
            break;
        end
    end
    iBaseBlockNum = i - 1;
    % direction vector
    vDirection = [1,-1,H,-H,H+1,H-1,-H+1,-H-1];
    for i = 1:iBaseBlockNum
        region1 = i == mLabel;
        se = strel('square', 3);
        dilatedRegion1 = imdilate(region1, se);
        mAdjacent = dilatedRegion1 - region1 > 0;
        mAdjacent(1,:) = 0; mAdjacent(H,:) = 0; mAdjacent(:,1) = 0; mAdjacent(:,W) = 0;
        vAdjacentIdx = find( mAdjacent );
        vRemoveIdx = [];
        for j = iBaseBlockNum + 1:length( stats )
            % compute the pixel index of outlier blocks in contact with reliable blocks
            vTempIdx = intersect(vAdjacentIdx,stats(j).PixelIdxList);
            if(~isempty(vTempIdx))
                %  Search for the direction in which the reliable block and outlier block contact, and obtain the reference phase          
                vPhaseOutlier = mPhaseInputOrg( vTempIdx );
                vPhaseBase = zeros(length(vTempIdx), 1);
                for m = 1:length(vTempIdx)
                    for direction = 1:8
                        if(i == mLabel( vTempIdx(m) + vDirection(direction) ))
                            vPhaseBase(m) = mPhaseInput( vTempIdx(m) + vDirection(direction) );
                            break;
                        end
                    end
                end
                % rectify the outlier block
                vPhaseDiffer = vPhaseBase - vPhaseOutlier;
                vOrder = vPhaseDiffer / (2*pi/FSet);
                order = round(mean(vOrder));
                OrderDiffer = mean(abs(vOrder - round(vOrder)));
                vContactRatio = length(vTempIdx) / stats(j).Area;
                if( OrderDiffer < mDiffer(vTempIdx(1)) && OrderDiffer < 0.2 && vContactRatio>0.05)
                    mPhaseInput(stats(j).PixelIdxList) = mPhaseInputOrg(stats(j).PixelIdxList) + (2*pi/FSet) * round( order );
                    mDiffer(stats(j).PixelIdxList) = OrderDiffer;
                end
            end
        end
    end
end
mPhase = mPhaseInput;
end

% This function employs flood fill algorithm to divide the connected domains by phase continuity
%
% Input parameters:
% mLabel          -- image mask, mark the valid area with 0 and mark the invalid area with -1
% mPhase          -- absolute phase map with noise
% F               -- frequency of fringe pattern
%
% Output:
% mLabel               -- absolute phase map after hierarchical adjustment
function [ mLabel ] = Func_PhaseFloodFillSegmentation( mLabel, mPhase, F )
[H, W] = size(mLabel);
dThresh = 0.5 / F;
mDirection = [0, 1;1, 0;0, -1; -1, 0;-1, -1;-1, 1; 1, -1; 1, 1];
vStack = nan( 300000, 2 );
    
iLabel = 1;
while(1)
    [iY, iX] = find( 0 == mLabel, 1 );
    % fill is complete if pixels without assigned number do not exist
    if(isempty(iY))
        break;
    end
    % assign stacks to points to be filled
    iEnd = 1;
    vStack( iEnd, : ) = [iX, iY];
    mLabel( iY, iX ) = iLabel;
    % fill the area numbered ilabel
    while( iEnd )
        iX = vStack(iEnd, 1);
        iY = vStack(iEnd, 2);
        vStack( iEnd, : ) = nan;
        iEnd = iEnd - 1;
        for i = 1:4
            iXnext = iX + mDirection( i, 1 );
            iYnext = iY + mDirection( i, 2 );
            % If the next phase continuous point has no number, it will be assigned a number
            if( iXnext > 0 && iXnext < W && iYnext > 0 && iYnext < H && 0 == mLabel( iYnext, iXnext ) )
                dPhaseGrad = abs(mPhase(iYnext, iXnext) - mPhase(iY, iX));
                bCriterion = min( dPhaseGrad, 2*pi - dPhaseGrad ) < dThresh;
%                 bCriterion = dPhaseGrad < dThresh;
                if( bCriterion )
                    mLabel( iYnext, iXnext ) = iLabel;
                    iEnd = iEnd + 1;
                    vStack( iEnd, : ) = [iXnext, iYnext];
                end
            end
        end
    end
    iLabel = iLabel + 1;
end
end

% This function conducts stereo matching based on the wrapped phase maps. This is a simplified version for improving speed. 
%
% Input parameters:
% mPhaseLeft           -- wrapped phase map of left camera
% mPhaseRight          -- wrapped phase map of right camera
% mDispMin, mDispMax   -- disparity range in auxiliary camera corresponding to each main camera pixel (x^m, y^m)
% iRadius              -- the radius of sum of absolute difference window   
% 
% Output:
% mMatch               -- the corresponding point of each main camera pixel in the auxiliary camera
function[ mMatch ] = Func_PhaseMatchSAD( mPhaseLeft, mPhaseRight, mDispMin, mDispMax, iRadius )
[iCameraHeight, iCameraWidth] = size(mPhaseLeft);
mMatch = nan( iCameraHeight, iCameraWidth );
iDiameter = 2*iRadius + 1;
mPhaseLeft(isnan(mPhaseLeft)) = 0;
mPhaseRight(isnan(mPhaseRight)) = 0;
iNum = iDiameter^2;
parfor i = 1 + iRadius:iCameraHeight - iRadius
    vMatch = nan( 1, iCameraWidth );
    for j = 1 + iRadius:iCameraWidth - iRadius    
        if( mPhaseLeft( i, j ) == 0 )
            continue;
        end
        mPhiTemplate = mPhaseLeft( i - iRadius:i + iRadius, j - iRadius: j + iRadius);
        iDispRange = mDispMax(i, j) - mDispMin(i, j);
        vCost = nan(iDispRange, 1);
        xFrom = j + mDispMin(i, j);
        xTo = j + mDispMax(i, j);
        iCount = 0;
        vXCandidate  = xFrom:xTo;
        
        vXBlock = - iRadius: iRadius;
        vYBlock = - iRadius: iRadius;
        for idx = xFrom: xTo
            iCount = iCount + 1;
            if( xFrom - iRadius < 1 || xTo + iRadius >= iCameraWidth )
                continue;
            end
            mPhiTarget = mPhaseRight( i - iRadius:i + iRadius, idx - iRadius: idx + iRadius);
            mCost = min( abs( mPhiTarget - mPhiTemplate ), 2.*pi - abs( mPhiTarget - mPhiTemplate ) );
            vCost(iCount) = sum(mCost(:)) ./ iNum;
        end
        iValid = sum(~isnan(vCost));
        if(iValid < 0.2*iDispRange)
            continue;
        end
        [ ~, iIdxMin ] = min( vCost );
        if( iIdxMin == 1 || iIdxMin == iDispRange + 1 )
            vMatch(j) = xFrom + iIdxMin - 1;
        else
            vVariable = vXCandidate( iIdxMin - 1:iIdxMin + 1 )';
            Temp = vVariable(1);
            vVariable = vVariable - Temp;
            vValue = vCost( iIdxMin - 1:iIdxMin + 1 );
            mA = [ vVariable.^2, vVariable, ones( 3, 1 ) ];
            vABC = mA \ vValue;
            dXMatch = - vABC( 2 ) / vABC( 1 ) / 2 + Temp;
            vMatch(j) = dXMatch;
        end
    end
    mMatch(i,:) = vMatch;
end
end

% This function conducts stereo matching based on the wrapped phase maps.
% This is a single pixel version with 1x1 stereo matching windows
%
% Input parameters:
% mPhaseLeft           -- wrapped phase map of left camera
% mPhaseRight          -- wrapped phase map of right camera
% mDispMin, mDispMax   -- disparity range in auxiliary camera corresponding to each main camera pixel (x^m, y^m)
% 
% Output:
% mMatch               -- the corresponding point of each main camera pixel in the auxiliary camera
function[ mMatch ] = Func_PhaseMatchPoint( mPhaseLeft, mPhaseRight, mDispMin, mDispMax )
[iCameraHeight, iCameraWidth] = size(mPhaseLeft);
mMatch = nan( iCameraHeight, iCameraWidth );
parfor i = 1:iCameraHeight
    vMatch = nan( 1, iCameraWidth );
    for j = 1:iCameraWidth
        dPhiTemplate = mPhaseLeft( i, j );
        if( isnan( dPhiTemplate ) )
            continue;
        end
        vPhiCandidate = mPhaseRight( i, : );
        xFrom = max( 1, j + mDispMin(i, j) );
        xTo = min( iCameraWidth, j + mDispMax(i, j) );
        iDispRange = xTo - xFrom;
        vCost = nan(iDispRange, 1);
        vXCandidate  = xFrom:xTo;
        vCost = min( abs( vPhiCandidate(xFrom:xTo) - dPhiTemplate ), 2.*pi - abs( vPhiCandidate(xFrom:xTo) - dPhiTemplate ) )';
        if( isempty(vCost) )
            continue;
        end

        iValid = sum(~isnan(vCost));
        if(iValid < 0.05*iDispRange)
            continue;
        end
        vCost(isnan(vCost)) = 99999;
        [ dCost, iIdxMin ] = min( vCost );      
        if( iIdxMin == 1 || iIdxMin == iDispRange + 1 )
            vMatch(j) = xFrom + iIdxMin - 1;
        else
            vVariable = vXCandidate( iIdxMin - 1:iIdxMin + 1 )';
            Temp = vVariable(1);
            vVariable = vVariable - Temp;
            vValue = vCost( iIdxMin - 1:iIdxMin + 1 );
            mA = [ vVariable.^2, vVariable, ones( 3, 1 ) ];
            vABC = Func_Inverse3(mA) * vValue;
            dXMatch = - vABC( 2 ) / vABC( 1 ) / 2 + Temp;
            vMatch(j) = dXMatch;
        end
    end
    mMatch(i,:) = vMatch;
end
end