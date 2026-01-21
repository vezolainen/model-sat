%% Vertex Coloring fo Sphere Delaunay Triangulation
%% such that no two adjacent vertices share the same color

clear; close all; clc;

% Initialize constellation
fprintf('Initializing constellation...\n');
%constellation = Constellation('Starlink');
%constellation = Constellation('singlePolarGlobal');
%constellation = Constellation('singlePolar');
%constellation = Constellation('biWalkerGlobal');
constellation = Constellation('singleDelta');

constellation.getInitialState();
epochs = 0;
constellation.propagateJ2(epochs);

N = size(constellation.state.eci, 1);
r = constellation.earthRadius;

points = zeros(N, 3);
epochIdx = 1;
for satIdx = 1:N
	points(satIdx,:) = constellation.state.eci(satIdx,:,epochIdx);
end
points = points ./sqrt(sum(points.^2, 2));

% Delaunay triangulation via convex hull
fprintf('Computing Delaunay triangulation via convex hull...\n');
tri = convhull(points(:, 1), points(:, 2), points(:, 3));
fprintf('Found %d triangles\n', size(tri,1));

% Adjacency matrix for the vertices
fprintf('Building adjacency matrix...\n');
adjacency = sparse(N, N);
for iTriangle = 1:size(tri, 1)
	v1 = tri(iTriangle, 1);
	v2 = tri(iTriangle, 2);
	v3 = tri(iTriangle, 3);
	
	adjacency(v1, v2) = 1;
	adjacency(v2, v1) = 1;
	
	adjacency(v1, v3) = 1;
	adjacency(v3, v1) = 1;
	
	adjacency(v2, v3) = 1;
	adjacency(v3, v2) = 1;
end
adjacency = adjacency - diag(diag(adjacency));

degrees = full(sum(adjacency, 2));
fprintf('Maximum degree: %d\n', max(degrees));
fprintf('Average degree: %.2f\n', mean(degrees));

function [colors, numColors] = dsaturColoring(adjacency, N)
	% DSATUR (Degree of Saturation) coloring algorithm
	colors = zeros(N, 1);
	
	% Initial degrees, saturation, uncolored vertices
	degrees = full(sum(adjacency, 2));
	saturation = zeros(N, 1);
	uncolored = true(N, 1);

	while any(uncolored)
		% Find uncolored vertex with maximum degree
		candidates = find(uncolored);
		max_sat = max(saturation(candidates));
		max_sat_vertices = candidates(saturation(candidates) == max_sat);
		
		if length(max_sat_vertices) > 1
			cand_degrees = degrees(max_sat_vertices);
			[~, idx] = max(cand_degrees);
			v = max_sat_vertices(idx);
		else
			v = max_sat_vertices(1);
		end
		
		% Find available colors for v1
		neighbours = find(adjacency(v, :));
		usedColors = unique(colors(neighbours(colors(neighbours) > 0)));
		
		% Assign smallest available color
		color = 1;
		while ismember(color, usedColors)
			color = color + 1;
		end
		colors(v) = color;
		uncolored(v) = false;
		
		%Update saturation of neighbours
		for u = neighbours
			if uncolored(u)
				% Check if the neighbour sees a new color
				nbrsU = find(adjacency(u, :));
				colorsAroundU = unique(colors(nbrsU(colors(nbrsU) > 0)));
				saturation(u) = length(colorsAroundU);
			end
		end
	end
	
	numColors = max(colors);
end

% Run coloring algorithm
fprintf('Using DSATUR coloring...\n');
[colors, numColors] = dsaturColoring(adjacency, N);

% Verify the coloring
fprintf('Verifying the coloring...\n');
valid = true;
for i = 1:N
	neighbours = find(adjacency(i, :));
	if any(colors(neighbours) == colors(i))
		fprintf('ERROR: Vertex %d shares color with neighbour!\n', i);
		valid = false;
		break;
	end
end

if valid
	fprintf('Coloring is valid. No adjacent vertices share the same color\n');
end

% Print results
for colorIdx = 1:numColors
	count = sum(colors == colorIdx);
	fprintf('  Color %2d: %4d vertices (%.1f%%)\n',...
	colorIdx, count, 100 * count / N);
end

% Visualization with random colormap and some statistics
colorMap = hsv(numColors);
shuffleOrder = randperm(numColors);
colorMap = colorMap(shuffleOrder, :);

vertexColors = colorMap(colors, :);

figure('Position', [100, 100, 1400, 600]);
subplot(1, 3, 1);
scatter3(points(:, 1), points(:, 2), points(:, 3), 20, vertexColors, 'filled');
axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title(sprintf('Vertex Coloring on Sphere\n%d Colors Used', numColors));
light('Position', [1, 1, 1]);
lighting gouraud;

subplot(1, 3, 2);
patch('Faces', tri, 'Vertices', points, ...
	'FaceVertexCData', vertexColors, ...
	'FaceColor', 'interp', ...
	'EdgeColor', 'k', ...
	'EdgeAlpha', 0.2);
axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Interpolated coloring\n');
light('Position', [1, 1, 1]);
lighting gouraud;

subplot(1, 3, 3);
colorCounts = histc(colors, 1:numColors);
bar(1:numColors, colorCounts);
xlabel('Color Index');
ylabel('Number of Vertices');
title(sprintf('Color Distribution\nMean: %.1f, Std: %.1f', ...
	mean(colorCounts), std(colorCounts)));
grid on;





