function p = planeThroughTrace(xyz)
% XYZ is a N-by-3 matrix of points defining a traces
[pt1, pt2, pt3] = selectPoints();
p = createPlane(pt1,pt2,pt3);
if any(isnan(p))
    keyboard
end

    function [p1,p2,p3] = selectPoints()
        % Create some random point sets, then select the point set that
        % defines the triangle with the most included area
        N = 20; % number of combinations;
        n = size(xyz,1);
        randsel = @(np) round(rand(np,1)*(n-1))+1;  %Randomly select ids
        idxlist = [randsel(N) randsel(N) randsel(N)]; %Random combinations
        A = NaN(N,1);
        for j = 1:N
            p1 = xyz(idxlist(j,1),:);
            p2 = xyz(idxlist(j,2),:);
            p3 = xyz(idxlist(j,3),:);
            a = sqrt( sum((p1-p2).^2) );     % length side a
            b = sqrt( sum((p2-p3).^2) );     % length side b
            c = sqrt( sum((p1-p3).^2) );     % length side c
            s = (a+b+c)/2;              % semi-perimeter
            A(j) = sqrt( s*(s-a)*(s-b)*(s-c) );
        end
        [~,id] = max(A);
        p1 = xyz(idxlist(id,1),:);
        p2 = xyz(idxlist(id,2),:);
        p3 = xyz(idxlist(id,3),:);
    end

end %planeThroughTrace()

