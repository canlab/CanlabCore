%++
%|
%|  Robert C. Welsh
%|  University of Michigan
%|  Department of Radiology 
%|
%|  A routine to mask an image given the input range.
%|
%|  Copyright 2000/01
%|
%|  Version 0.8 - October 19, 2000
%|
%|     function [maskedImage maskingImage]] = ...
%|             maskImg(inputImage,maskLow,maskHi)
%|
%--
 
function [maskedImage, maskingImage] = maskImg(inputImage,maskLow,maskHi)
 
maskedLowImage = max(0,inputImage - maskLow);
 
maskedHiImage  = max(0,inputImage - maskHi);
 
maskingImage = xor(maskedLowImage,maskedHiImage);
 
maskedImage = maskingImage.*inputImage;
 
%
% End
%
return
