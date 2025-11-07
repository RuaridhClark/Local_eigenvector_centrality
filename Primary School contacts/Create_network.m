%% Create weighted adjacency matrix from GEXF-like XML
filename = 'sp_data_school_day_2_g.gexf';  % <-- change to your actual file name

% Read XML file
xml = xmlread(filename);

% --- Extract Nodes ---
nodeList = xml.getElementsByTagName('node');
numNodes = nodeList.getLength;
nodeIDs = zeros(numNodes,1);

for k = 0:numNodes-1
    node = nodeList.item(k);
    nodeIDs(k+1) = str2double(char(node.getAttribute('id')));
end

% Create a mapping from original IDs to consecutive indices
[~, ~, nodeIndexMap] = unique(nodeIDs);
numNodes = numel(nodeIDs);

% --- Extract Edges ---
edgeList = xml.getElementsByTagName('edge');
numEdges = edgeList.getLength;

% Initialize sparse adjacency matrix
A = sparse(numNodes, numNodes);

for e = 0:numEdges-1
    edge = edgeList.item(e);
    sourceID = str2double(char(edge.getAttribute('source')));
    targetID = str2double(char(edge.getAttribute('target')));
    
    % Handle optional weight attribute (duration or count)
    weight = 1; % default weight
    if edge.hasAttribute('weight')
        weight = str2double(char(edge.getAttribute('weight')));
    else
        % Check if duration attribute exists
        attvalues = edge.getElementsByTagName('attvalue');
        for a = 0:attvalues.getLength-1
            att = attvalues.item(a);
            if strcmp(char(att.getAttribute('for')), '2') % duration
                weight = str2double(char(att.getAttribute('value')));
                break;
            end
        end
    end
    
    % Map original IDs to matrix indices
    i = nodeIndexMap(nodeIDs == sourceID);
    j = nodeIndexMap(nodeIDs == targetID);
    
    % Fill adjacency matrix (undirected)
    A(i,j) = weight;
    A(j,i) = weight;
end

% Convert sparse to full (optional)
A = full(A);

%% Display results
fprintf('Nodes: %d\n', numNodes);
fprintf('Edges: %d\n', numEdges);
disp('Sample adjacency submatrix:');
disp(A(1:min(5,numNodes), 1:min(5,numNodes)));

%% Save adjacency matrix (optional)
% writematrix(A, 'adjacency_from_gexf.csv');
