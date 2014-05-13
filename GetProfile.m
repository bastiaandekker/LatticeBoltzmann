function [ data, fit, y] = GetProfile( fileName, tau, deltaV)
%UNTITLED Summary of thisdata function goes here
%   Detailed explanation goes here
data=importdata(fileName); 
lY=size(data,1);

%% Calculate fitting parabola and Plot data
visc = (2*tau - 1)/6;
c = 1;
curvature = c*deltaV/visc;
y = -(lY-1)/2:(lY-1)/2;
fact = 0.5;
fit = -fact*curvature*(y).^2 + fact*curvature*(y(2)-1/2)^2; 

end

