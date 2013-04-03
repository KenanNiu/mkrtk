function marm = calcMomentArms(marm,hax,loa)
% CALCMOMENTARMS Calculate moment arms from helical axis and line of action
%
%   Called by:
%       Registration > analysisUpdater
%       Registration > MI_TendonLoA_Callback
%       

for j = 1:numel(marm)
    % Get the appropriate line of action and helical axis for calculating
    % moment arm (j): 
    [line1,line2] = find_structs_with_tag(marm.Item1,marm.Item2,hax,loa);
    [p1,p2] = points_from_struct(line1);
    [p3,p4] = points_from_struct(line2);
    
    % Calculate the moment arm for each phase:
    for k = size(p1,1):-1:1
        [d(k,1),P1(k,:),P2(k,:)] = momentArm(p1(k,:),p2(k,:),p3(k,:),p4(k,:));
    end
    
    % Store calcs:
    marm(j).Point1 = P1;
    marm(j).Point2 = P2;
    marm(j).Length = d;
    
end



% ---------------------------------------
function [p1,p2] = points_from_struct(s)
% Return two points from the structure
%
% Input structure S could either have the fields:
%   -Point
%   -Vector
% or
%   -Point
%   -Axis
%
% We use the point and the "vector" or "axis" to return two points

assert(numel(s) == 1,'POINTS_FROM_STRUCT does not hande arrays of structs')

p1 = s.Point;

if isfield(s,'Vector')
    v = s.Vector;
else
    v = s.Axis;
end

p2 = p1 + v;
