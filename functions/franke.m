function [shape] = franke(N,d)
% This function generates the Franke function and its derivatives
% for a given N and d.
%
% Inputs:
%   N: Number of centers
%   d: Diameter of the circle enclosing centeres.
%
shape= 0.8*N.^0.25/d;

end