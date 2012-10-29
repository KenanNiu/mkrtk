function test_suite = testall %#ok<STOUT> 
%TESTALL Main testing function
%
% Requires the xUnit framework to perform tests

initTestSuite;



% Functions:
test_ColourSpec;

% Classes:
roi_class_spec;

% Behaviours / user input
roi_user_interaction_test;

