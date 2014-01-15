function theta = wrap_to_range(theta,tmin,tmax)
%WRAP_TO_RANGE Wrap THETA to specified range [TMIN TMAX]
dt = tmax-tmin;
theta = mod( theta+tmin, dt ) + tmin;
theta(theta==tmin) = tmax;

