function [Node, Panel] = ReadOBJ(file)
if ~exist(file,'file')
        error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
               ', be sure to specify the full path to the file.'], file);
end

fid = fopen(file,'r');   
Node = nan(1,3); 
Panel = cell(1);
cntv = 0; cntf = 0;

% parse .obj file 
while 1    
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end  % exit at end of file 
    ln = sscanf(tline,'%s',1); % line type 
     %disp(ln)
    switch ln
        case 'v'   % mesh vertexs
            cntv = cntv+1;
            Node(cntv,:) = sscanf(tline(2:end),'%f')';
        case 'f'   % face grouping
            cntf = cntf+1;
            str = textscan(tline(2:end),'%s'); 
            nn = length(strfind(str{1},'/'));
            % number of vertices
            stri = str{1};
            nf = length(strfind(stri{1},'/'))+1;
            % number of fields with this face vertices
            allind = textscan(tline(2:end),'%f',nf*nn,'Delimiter','/'); allind = allind{1};
            fv = allind(1:nf:end);
            Panel{cntf} = fv';
    end
end
Panel = Panel';

fclose(fid);