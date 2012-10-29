function [R,T,xyzout,err] = register(target,source,solver)
%REGISTER Register the hi-res source point cloud with low-res target
%
% Inputs
%   TARGET  3-by-M matrix of points which are stationary
%   SOURCE  3-by-N matrix of points to move onto TARGET
%   SOLVER  A string from any one of the options in the switch statement
%
% Outputs
%       R   3-by-3 rotation matrix
%       T   1-by-3 translation vector
%  xyzout   The points of SOURCE superimposed optimally on TARGET
%     err   An error measure (if calculted by the solver specified)
%

t = tic;

% Default solver:
if nargin < 3
    solver = 'libicp';
end

% Call selected solver:
switch solver
    case 'aditeration'
        assert(2==exist('aditeration','file'),'ADITERATION.M not on the path')
        
        % The homogeneous transformation matrix produced by ADITERATION is
        % of a different form, so it's easier to just fall back on the
        % procrustes call down the bottom to recover R & T.
        [xyzout, ~, err] = aditeration(target,source);
        
    case 'Kjer & Wilm ICP'
        assert(2==exist('icp','file'),'ICP.M by Kjer & Wilm not on the path')
        
        TOL = 1e-7;
        err = [inf;0];
        interim = source;
        while abs(diff(err)) > TOL
            [TR,TT,err] = icp(target',interim',2);%   'WorstRejection',0.2);
            interim = ( TR*interim' + repmat(TT,1,size(interim',2)) )';
        end
        xyzout = interim;
        err = err(end);
        
        
    case 'Kroon Finite ICP'
        assert(2==exist('ICP_finite','file'),'ICP_finite.m by D.Kroon not on the path');
        
        % Options.Registration: 'Rigid', Translation and Rotation (default)
        %         'Size', Rigid + Resize
        %         'Affine', Translation, Rotation, Resize
        %                       and Shear.
        % Options.TolX: Registration Position Tollerance, default is the
        %       largest side of a volume containing the points divided by 1000
        % Options.TolP: Allowed tollerance on distance error default
        %       0.001 (Range [0 1])
        % Options.Optimizer : optimizer used, 'fminlbfgs' (default)
        %       ,'fminsearch' and 'lsqnonlin'.
        % Options.Verbose : if true display registration information (default)
        opts.Optimizer = 'lsqnonlin';        
        [xyzout,H] = ICP_finite(target,source,opts);
        R = H(1:3,1:3); 
        T = H(1:3,4)';
        err = NaN;
        
    case 'libicp'
        assert(3==exist('icpMex','file'),'icpMex matlab mex file not found on the path')
        assert(3==exist('sparsifyMex','file'),'sparsifyMex matlab mex file not found on the path')
        
        % LIBICP uses a fixed inlier threshold distance to cull outliers.
        %   We could try making a version which set the inlier criteria to
        %   2.5*sigma according to Masuda 1996.
        indist = 10;                % -1 == no outlier rejection; >0 == inlier distance  
        Hinit = eye(4);             % 4-by-4, initial transformation matrix
        method = 'point_to_point';  % 'point_to_point' | 'point_to_plane'
        H = icpMex(target',source',Hinit,indist,method);
        R = H(1:3,1:3);
        T = H(1:3,4)';
        % Create these only if called for:
        if nargout > 2
            xyzout  = (R*source' + T*ones(1,size(source',2)) )';
            err = NaN;
        end
        
       
    case 'Gaussian Mixtures'
        % We have two options here
        %   - compiled binary (Windows only at the moment)
        %   - matlab mex implementation
        
        xyzout = gmmreg_wrapper(target,source);
                
    otherwise
        error('Solver option not recognised: %s',solver)

end 

% Now if the solver did not give us R and T, we calculate them directly:
if exist('R','var')~=1  ||  exist('T','var')~=1
    [~,~,tform] = procrustes( xyzout, source, 'scaling', false, 'reflection', false);
    R = tform.T';
    T = tform.c(1,:);
end

toc(t)


