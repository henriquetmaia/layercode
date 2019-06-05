function [r,b] = rbSegmentRGB(RGB)

orighsv = rgb2hsv(RGB);
rHmin = 0.946;
rHmax = 0.049;
rSmin = 0.311;
rSmax = 0.539;
rVmin = 0.629;
rVmax = 0.858;
r = ( (orighsv(:,:,1) >= rHmin) | (orighsv(:,:,1) <= rHmax) ) & ...
    (orighsv(:,:,2) >= rSmin ) & (orighsv(:,:,2) <= rSmax) & ...
    (orighsv(:,:,3) >= rVmin ) & (orighsv(:,:,3) <= rVmax);
bHmin = 0.564;
bHmax = 0.888;
bSmin = 0.000;
bSmax = 0.311;
bVmin = 0.393;
bVmax = 0.629;
b = ( (orighsv(:,:,1) >= bHmin) | (orighsv(:,:,1) <= bHmax) ) & ...
    (orighsv(:,:,2) >= bSmin ) & (orighsv(:,:,2) <= bSmax) & ...
    (orighsv(:,:,3) >= bVmin ) & (orighsv(:,:,3) <= bVmax);

end

