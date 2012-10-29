function scale = estimate_scale(M)

%%=====================================================================
%% Project:   Point Set Registration using Gaussian Mixture Model
%% Module:    $RCSfile: estimate_scale.m,v $
%% Language:  MATLAB
%% Author:    $Author: bing.jian $
%% Date:      $Date: 2008-11-14 08:34:29 +1100 (Fri, 14 Nov 2008) $
%% Version:   $Revision: 109 $
%%=====================================================================

[m,d] = size(M);
scale = det(M'*M/m);
for i=1:d
    scale = sqrt(scale);
end