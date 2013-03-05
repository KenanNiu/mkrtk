function test_suite = testall %#ok<STOUT> 
%TESTALL Main testing function
%
% Requires the xUnit framework to perform tests

initTestSuite;



% Functions:
test_ColourSpec;

% Classes:
roi_class_spec;
cloud_class_spec;
quaternion_classspec;

selectiveBrush_class_spec;
figlocker_test;

%bone_class_spec;  %% In progress


% Behaviours / user input
roi_user_interaction_test;

