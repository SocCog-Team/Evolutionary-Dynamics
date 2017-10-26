function [meanX, stdX] = weighted_mean_std(x, w)
  sumW = sum(w);
  meanX = x*w/sumW;
  dlt = bsxfun(@minus, x, meanX);
  stdX = sqrt((dlt.^2)*w/sumW);
end