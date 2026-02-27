function [Avu, Avxux, Avxuy, Avyuy] = init_stiff(vertices, meshes, T, P2P2, P2xP2x, P2xP2y, P2yP2y)
Nv = size(vertices, 1);
Avu = zeros(Nv);
Avxux = zeros(Nv);
Avxuy = zeros(Nv);
Avyuy = zeros(Nv);
for i = 1: size(meshes, 1)
    v1 = vertices(meshes(i, 1), :); %1-by-2
    v2 = vertices(meshes(i, 2), :);
    v3 = vertices(meshes(i, 3), :); 
    Jacobi = T(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    Avu(meshes(i, :), meshes(i, :)) = Avu(meshes(i, :), meshes(i, :)) + Jacobi*P2P2();
    Avxux(meshes(i, :), meshes(i, :)) = Avxux(meshes(i, :), meshes(i, :)) + Jacobi*P2xP2x(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    Avxuy(meshes(i, :), meshes(i, :)) = Avxuy(meshes(i, :), meshes(i, :)) + Jacobi*P2xP2y(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
    Avyuy(meshes(i, :), meshes(i, :)) = Avyuy(meshes(i, :), meshes(i, :)) + Jacobi*P2yP2y(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2));
end
Avu = sparse(Avu);
Avxux = sparse(Avxux);
Avxuy = sparse(Avxuy);
Avyuy = sparse(Avyuy);
end