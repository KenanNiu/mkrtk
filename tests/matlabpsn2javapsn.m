function jpsn = matlabpsn2javapsn(mpsn)
% Convert matlab screen positions to java screen positions
%
% Matlab: origin is lower left
% Java:   origin is upper left

ss = get(0,'ScreenSize');
jpsn = zeros(size(mpsn));               % Initialise
jpsn(:,1) = mpsn(:,1);              % Zero-base
jpsn(:,2) = ss(4) - (mpsn(:,2));    % Zero-base & invert

end %matlabpsn2javapsn()