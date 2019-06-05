
function maskedImage = maskImage( originalImage, mask )
    maskedImage = bsxfun(@times, originalImage, cast(mask,class(originalImage)));