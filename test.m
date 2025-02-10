addpath('domain/')
clear all;
clear;

%First define the domain
% Unit square domain
square_domain = Square([0, 1, 0, 1]);
square_domain = square_domain.generateGrid(20);
%square_domain = square_domain.generateBoundaryPoints(20);
square_domain.scatterPlot("Square Domain", true, false);
