function [N, vertices, meshes, Gamma, trans, G2m, e] = init_meshes(h, R)
model = createpde(1);
geometryFromEdges(model, decsg(R));
mesh = generateMesh(model, 'Hmax', h, 'GeometricOrder', 'linear');
N = size(mesh.Nodes, 2);
mesh = generateMesh(model, 'Hmax', h);
[p, e, t] = meshToPet(mesh);
%pdemesh(p, e, t, 'NodeLabels', 'on', 'ElementLabels', 'on');
vertices = p';
meshes = t([1: 3, 5, 6, 4], :)';
idx = find(e(5, :) == 1);
Gamma = [e(1, idx(1): 2: idx(end)-1)', e(2, idx(1)+1: 2: idx(end))', e(1, idx(1)+1: 2: idx(end))'];
trans = zeros(size(vertices, 1), N);
for i = 1: size(meshes, 1)
    trans(meshes(i, :), meshes(i, 1: 3)) = [eye(3); ones(3)/2 - eye(3)/2];
end
trans = sparse(trans);
G2m = zeros(size(Gamma, 1), 1);
for i = 1: size(meshes, 1)
    if vertices(meshes(i, 4), 2) == 1/2
        G2m(Gamma(:, 3) == meshes(i, 4)) = i;
    end
end
assert(min(G2m) >= 1);
end