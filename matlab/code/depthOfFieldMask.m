
function [dofMask] = depthOfFieldMask(I)
%%%%% for a 2k image, to help deal with
%%%%% shadows/foreground-background/DepthOfField%
%%%%% made with our Database colors in mind, although
%%%%% highlights and shadows are filtered otherwise
%%%%% by rescaling the image and removing local outliers

channel1Min = 0.000;
channel1Max = 255.000;
channel2Min = 0.000;
channel2Max = 118.000;
channel3Min = 1.000;
channel3Max = 131.000;
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;
seLarge = strel('square', 4);
dofMask = imopen( BW, seLarge );

minNumPix = 200;
conn = 4;
dofMask = bwareaopen(dofMask, minNumPix, conn);
dofMask = ~dofMask;
dofMask = bwareaopen(dofMask, minNumPix, conn);
dofMask = ~dofMask;


end

