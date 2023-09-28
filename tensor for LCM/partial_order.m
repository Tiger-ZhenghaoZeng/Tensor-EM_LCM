function [gamma] = partial_order(theta, tole)
% Compute the partial order matrix based on item parameters theta
% Inputs: theta the estimated item parameters, tole the deviation allowed
% for a class to be considered as mastering the skills needed for an item
% Output: a matrix gamma indicating the partial order among different
% classes
max_mat = max(theta, [], 2);
gamma = abs(theta-max_mat) < tole;
end

