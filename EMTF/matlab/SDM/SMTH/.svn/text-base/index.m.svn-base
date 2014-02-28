function [ii,jj,button] = index(X,Y)
%  for a plot of an n x m array, return the index of the
%  element clicked on with the mouse + the button clicked
%  X is an n+1 vector giving the limits in the x direction of
%   the plotted cells, Y is an m+1 vector giving the limits in the
%  y direction of the cells
n1 = length(X);m1 = length(Y);
waitforbuttonpress
point = get(gca,'CurrentPoint');
x = point(1,1);
y = point(1,2);
button = get(gcf,'SelectionType');

ii = sum(X < x);
jj = sum(Y < y);
