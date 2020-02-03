function [NODE,PANEL,BDRY] = ConfigEggbox(x_divs,y_divs,theta,a,b,gmma)
theta=pi*theta/180;
gmma = gmma*pi/180;
if gmma>theta, error('This configuration does not exist!'); end
numx = 2*x_divs; numy = 2*y_divs;
cosb = cos(theta)/cos(gmma); 
h1 = a*cos(gmma); h2 = b*cosb;
h = h1+h2;
s = b*sqrt(1-cosb^2);
l = a*sin(gmma);
[X,Y] = meshgrid(linspace(0,s*numx,numx+1),linspace(0,l*numy,numy+1));
Z = 0*X;
Z(2:2:end,1:2:end) = h1; Z(1:2:end,2:2:end) = h2; 
Z(2:2:end,2:2:end) = h;
NODE = [reshape(X,[],1),reshape(Y,[],1),reshape(Z,[],1)];

BDRY = zeros(2*numx+2*numy,2);
Node_idx = reshape(1:size(NODE,1),size(X));
Lbdry = Node_idx(:,1);
Rbdry = Node_idx(:,end);
Bbdry = Node_idx(1,:);
Tbdry = Node_idx(end,:);

% Specify boundary
BDRY(1:numel(Lbdry)-1,:) = [Lbdry(1:end-1) Lbdry(2:end)];
count = numel(Lbdry)-1;
BDRY(count+1:count+numel(Bbdry)-1,:) = [Bbdry(1:end-1)' Bbdry(2:end)'];
count = count+numel(Bbdry)-1;
BDRY(count+1:count+numel(Rbdry)-1,:) = [Rbdry(1:end-1) Rbdry(2:end)];
count = count+numel(Rbdry)-1;
BDRY(count+1:count+numel(Tbdry)-1,:) = [Tbdry(1:end-1)' Tbdry(2:end)'];
count = count+numel(Tbdry)-1;

k = 0; PANEL = cell(numx*numy,1);
for j=1:numy, for i=1:numx
        k = k+1;
        n1 = (i-1)*(numy+1)+j; n2 = i*(numy+1)+j;
        PANEL{k} = [n1 n2 n2+1 n1+1];
end, end

