function [] = PlotOri(Node,Panel,Trigl,varargin)
if length(varargin)<1, varargin = {'FaceVertexColor',[]}; end
options = POParseOptions(varargin{:});
if strcmpi(options.ShowNumber,'on')
    options.PanelColor = 'none';
    options.BendEdgeStyle = ':';
end
if ~isempty(options.FaceVertexColor)
    if size(options.FaceVertexColor,1)==size(Node,1)
        options.PanelColor = 'interp';
    else
        options.PanelColor = 'flat';
    end
end
if isempty(options.EdgeColor)
    if ~isempty(Trigl)
        patch('faces', Trigl, 'vertices', Node, 'facecolor', options.PanelColor, ...
              'linestyle', options.BendEdgeStyle, 'facelighting', 'flat',...
              'edgecolor', (1-options.EdgeShade)*[1 1 1], 'FaceVertexCData',options.FaceVertexColor);
    end
    hold on;
    panelsize = cellfun(@numel,Panel);
    panels = nan(length(Panel),max(panelsize));
    for i = 1:length(Panel), panels(i,1:panelsize(i)) = Panel{i}; end

    if ~isempty(Trigl)
        patch('faces', panels, 'vertices', Node, 'facecolor', 'none', 'facelighting', 'flat',...
              'linestyle', options.FoldEdgeStyle, 'linewidth', 1, 'edgecolor', (1-options.EdgeShade)*[1 1 1]);
    else
        patch('faces', panels, 'vertices', Node, 'facecolor', options.PanelColor, 'facelighting', 'flat',...
              'linestyle', options.FoldEdgeStyle, 'linewidth', 1, 'edgecolor', (1-options.EdgeShade)*[1 1 1],...
              'FaceVertexCData',options.FaceVertexColor);
    end
    
    if strcmpi(options.ShowNumber,'on')
        for i=1:size(Node,1) 
            text(Node(i,1)+0.1,Node(i,2)-0.1,Node(i,3),num2str(i),'Fontsize',14);
        end
    end
    
elseif ~isempty(options.EdgeColor)
    if isempty(Trigl)
        error('Edge Coloring Mode requires triangulation information!')
    else
        patch('faces', Trigl, 'vertices', Node, 'facecolor', 0.85*[1 1 1], ...
              'facelighting', 'flat','edgecolor', 'none');
    end
    hold on
    for i=1:options.NumBendHinge 
        XYZ = [Node(options.Bars(i,1),:) Node(options.Bars(i,2),:)];
        plot3(XYZ([1,4])',XYZ([2,5])',XYZ([3,6])',':','LineWidth',1.5,'Color',options.EdgeColor(i,:))
    end
    for j=(options.NumBendHinge+1):size(options.Bars,1) 
        XYZ = [Node(options.Bars(j,1),:) Node(options.Bars(j,2),:)];
        plot3(XYZ([1,4])',XYZ([2,5])',XYZ([3,6])','-','LineWidth',2,'Color',options.EdgeColor(j,:))
    end
end
end

function options = POParseOptions(varargin)
IP = inputParser;
IP.addParameter('FoldEdgeStyle', '-', @ischar)
IP.addParameter('EdgeShade', 1, @isnumeric);
IP.addParameter('PanelColor', 'g');
IP.addParameter('BendEdgeStyle', 'none', @ischar)
IP.addParameter('ShowNumber', 'off', @ischar)
IP.addParameter('FaceVertexColor', [], @isnumeric)
IP.addParameter('EdgeColor', [], @isnumeric)
IP.addParameter('Bars', [], @isnumeric)
IP.addParameter('NumBendHinge', 0, @isnumeric)
IP.parse(varargin{:});
options = IP.Results;
end