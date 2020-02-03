function [] = Write2OBJ(filename, Node, Trigl, Bars, bend)
sizebend = size(bend,1);
vp = [max(Node(:,1))/2,max(Node(:,2))/2,10*(max(Node(:,3))+1)]; % specify a view point to unify vertex winding
%% write data
fileID = fopen([filename,'.obj'],'w');
fprintf(fileID,'# This is an output from MERLIN2\n');

% write vertex data
fprintf(fileID,'v %.10g %.10g %.10g\n',Node');
fprintf(fileID,'# %d vertices\n', size(Node,1));

% write polygonal panel data
for i=1:size(Trigl,1)
    a = Node(Trigl(i,2),:) - Node(Trigl(i,1),:);
    b = Node(Trigl(i,3),:) - Node(Trigl(i,2),:);
    np = sum(Node(Trigl(i,:),:),1)/3-vp;
    % Unify vertex winding order
    if cross(a,b)*np'<0
        Trigli = Trigl(i,:);
    else
        Trigli = Trigl(i,[1 3 2]);
    end
    fprintf(fileID,'f %d %d %d\n', Trigli);
end
% The vertex-rewinding algorithm only works for simple surface meshes
fprintf(fileID,'# %d triangles\n', size(Trigl,1));

if nargin>3
    % write edge data
    fprintf(fileID,'#e %d %d 1\n', Bars(1:sizebend,:)');
    fprintf(fileID,'# %d bending lines\n', sizebend);
    fprintf(fileID,'#e %d %d 0\n', Bars((sizebend+1):end,:)');
    fprintf(fileID,'# %d generic (folding/boundary) edges\n', size(Bars,1)-sizebend);
end

fclose(fileID);