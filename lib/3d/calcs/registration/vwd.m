function [transformation_matrix, transformation_error] = vwd(initial_position, final_position)
% This function is intended to replace the VWD function from UBC written by
% Dave Wilson, Nick Hill, and Jeff Vavasour, and included with Agnes' knee
% program.
%
% The old VWD function would occasionally produce imaginary numbers which
% caused the program to break down.  
%
% This function instead uses MATLAB's PROCRUSTES function (part of the
% Statistics Toolbox), which works quite nicely.  We then have to format
% the output to be equivalent to the form of the old VWD function.
%
% Joshua Martin, 2012

[d,Z,tform] = procrustes(final_position,initial_position,'scaling',false,'reflection',false);

% Now convert to the form required by Agnes' code:
rotation_matrix = tform.T';
translation_vector = tform.c(1,:);



% Now use code basically ripped from old VWD.M with minor improvements:


% Finding the Transformation Error:
% ---------------------------------

for i = size(initial_position,1) : -1 : 1
    transformed_initial_matrix(i,:) = [translation_vector + (rotation_matrix*(initial_position(i,:))')'];
end

point_error = transformed_initial_matrix - final_position;

total_error = 0;

max_error = 0.0000000001;

for j = 1:size(initial_position,1)
    total_error = total_error + norm(point_error(j,:));
    if norm(point_error(j,:)) > max_error
        max_error = norm(point_error(j,:));
    end
end

transformation_error = total_error/size(initial_position,1);


% Creating the Transformation Matrix:
% -----------------------------------

transformation_matrix(1,:) = [1 0 0 0];
transformation_matrix(2:4,1) = translation_vector';
transformation_matrix(2:4,2:4) = rotation_matrix;
