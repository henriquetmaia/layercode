function [labeled, nred, nblue] = labelRBPix( filename, withFilter )

[RGB,~,AlphaChannel] = imread(filename);

%% HSV picker
[r, b] = rbSegmentRGB( RGB );

%% Label reliable pixels

% shrink borders to avoid flood leaks
mask = imerode(AlphaChannel,strel('square',7));
mask(mask~= 255) = 0;

% Remove shaded pixels
if withFilter
    dofMask = depthOfFieldMask(RGB);
else
    dofMask = zeros(size(AlphaChannel));
end
% fim(dofMask);


% Remove tiny regions
minNumPix = floor(size(AlphaChannel,1) / 20);
conn = 4;
cleanedR = bwareaopen(r, minNumPix, conn);
cleanedB = bwareaopen(b, minNumPix, conn);

cleanedR = cleanedR & ~cleanedB & mask & ~dofMask;
cleanedB = cleanedB & ~cleanedR & mask & ~dofMask;

cleanedR = bwareaopen(cleanedR, minNumPix, conn);
cleanedB = bwareaopen(cleanedB, minNumPix, conn);

% Label colors with unique indices
[rlabel, nred] = bwlabel( cleanedR, conn );
[blabel, nblue] = bwlabel( cleanedB, conn );
tmp = (blabel == 0);
blabel = blabel + nred;
blabel(tmp) = 0;
labeled = blabel + rlabel;

end

