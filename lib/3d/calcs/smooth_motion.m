function mdls = smooth_motion(mdls)
% This function smooths the motion of all the models in MDLS
%

addpath(genpath('tendon_loa_wip'))


% Get the proxy for joint angle:
disp('Calculating a variable required for smoothing...')
proxy = joint_dt_proxy(mdls);

%proxy = [1:numel(proxy)]';
disp('Smoothing motion')
mdls = smooth_models(mdls,proxy);


% ------------------------------------------------------------------------
function mdls = smooth_models(mdls,theta)

% See also tendon_loa.m for dependency on filter_width
global filter_width
filter_width = 0.99;
width = filter_width;

% if matlabpool('size') == 0
%     matlabpool open 4
% end
n = numel(mdls);
for j = 1:n     % ==> upgrate to parfor
    if ~isempty(mdls(j).q) && ~all(isnan(mdls(j).q))
        qj = mdls(j).qraw.unwrap;   %\_ raw data
        %qj = mdls(j).qraw;
        xj = mdls(j).xraw;          %/
        
        sfun = @(y)smooth(theta,y,width,'rloess');
        
        qj = qj.smooth(sfun,1);
        for c = 1:3
            xj(:,c) = sfun(xj(:,c));
        end
        
        mdls(j).q = qj;
        mdls(j).x = xj;
    end
end


