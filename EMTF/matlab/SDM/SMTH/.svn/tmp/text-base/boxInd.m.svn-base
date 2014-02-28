%  box_ind  using rubber-band box, find enclosed range of grid indices
%  for array of size neig_plt,nb
k = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');
finalRect = rbbox;
point2 = get(gca,'CurrentPoint');
ii = [point1(1,1) point2(1,1)];
jj = [point1(1,2) point2(1,2)];
ii = sort(ii);
jj = sort(jj);
button = get(gcf,'SelectionType');

i1 = floor(ii(1));
i1 = max(1,i1);
i2 = floor(ii(2));
i2 = min(nbt,i2);

j1 = floor(jj(1));
j1 = max(1,j1);
j2 = floor(jj(2));
j2 = min(neig_plt,j2);


