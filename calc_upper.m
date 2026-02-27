function [AvSuS, AvLuL, AvSuL, AvsusG, AvuxcnG, AvuycnG, Avn1S, Avn2S, Avn1L, Avn2L, bS1, bS2] = calc_upper(shuttle, epsilon, verticesS, meshesS, GammaS, G2mS, verticesL, meshesL, GammaL, G2mL, c, gamma, P2P2xP2nG, P2P2yP2nG)
n = size(GammaL, 1);
NvS = size(verticesS, 1);
NvL = size(verticesL, 1);
AvSuS= zeros(NvS);
AvLuL = zeros(NvL);
AvSuL = zeros(NvS, NvL);
AvsusG = zeros(NvS);
AvuxcnG = zeros(NvS, NvL);
AvuycnG = zeros(NvS, NvL);
Avn1S = zeros(NvS, n);
Avn2S = zeros(NvS, n);
Avn1L = zeros(NvL, n);
Avn2L = zeros(NvL, n);
bS1 = zeros(NvS, 1);
bS2 = zeros(NvS, 1);
for i = 1: n
    v1 = verticesL(GammaL(i, 1), :);
    v2 = verticesL(GammaL(i, 2), :);
    len = norm(v2 - v1, 2);
    n1 = (v2(2) - v1(2))/len;
    n2 = (v1(1) - v2(1))/len;    
    switch shuttle
        case 1
            lambda = 1;
        case 2
            lambda = 2*len*n;
        case 3
            lambda = 2*len*n - 1 - log(len*n);
        otherwise
            assert(false, int2str(shuttle));
    end
    AvSuS(GammaS(i, :), GammaS(i, :)) = AvSuS(GammaS(i, :), GammaS(i, :)) + len*[4, -1, 2; -1, 4, 2; 2, 2, 16]/30;
    AvLuL(GammaL(i, :), GammaL(i, :)) = AvLuL(GammaL(i, :), GammaL(i, :)) + len*[4, -1, 2; -1, 4, 2; 2, 2, 16]/30;
    AvSuL(GammaS(i, :), GammaL(i, :)) = AvSuL(GammaS(i, :), GammaL(i, :)) + len*[4, -1, 2; -1, 4, 2; 2, 2, 16]/30;
    AvsusG(GammaS(i, :), GammaS(i, :)) = AvsusG(GammaS(i, :), GammaS(i, :)) + lambda*reshape([37, -3, 36; 7, 7, -4; -44, -4, -32; 7, 7, -4; -3, 37, 36; -4, -44, -32; -44, -4, -32; -4, -44, -32; 48, 48, 64]*gamma(i, :)', [3, 3])/30/len;
    Avn1S(GammaS(i, :), i) = len*n1*[1; 1; 4]/6;
    Avn2S(GammaS(i, :), i) = len*n2*[1; 1; 4]/6;
    Avn1L(GammaL(i, :), i) = len*n1*[1; 1; 4]/6;
    Avn2L(GammaL(i, :), i) = len*n2*[1; 1; 4]/6;
    bS1(GammaS(i, :)) = bS1(GammaS(i, :)) - n1*epsilon*[37, 7, -44, 7, -3, -4, -44, -4, 48; -3, 7, -4, 7, 37, -44, -4, -44, 48; 36, -4, -32, -4, 36, -32, -32, -32, 64]*reshape(c(GammaL(i, :))*c(GammaL(i, :))', [9, 1])/30/len;
    bS2(GammaS(i, :)) = bS2(GammaS(i, :)) - n2*epsilon*[37, 7, -44, 7, -3, -4, -44, -4, 48; -3, 7, -4, 7, 37, -44, -4, -44, 48; 36, -4, -32, -4, 36, -32, -32, -32, 64]*reshape(c(GammaL(i, :))*c(GammaL(i, :))', [9, 1])/30/len;
    v1 = verticesL(meshesL(G2mL(i), 1), :); %1-by-2
    v2 = verticesL(meshesL(G2mL(i), 2), :);
    v3 = verticesL(meshesL(G2mL(i), 3), :); 
    AvuxcnG(meshesS(G2mS(i), :), meshesL(G2mL(i), :)) = AvuxcnG(meshesS(G2mS(i), :), meshesL(G2mL(i), :)) + len*reshape(P2P2xP2nG(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2))*c(meshesL(G2mL(i), :)), [6, 6]);
    AvuycnG(meshesS(G2mS(i), :), meshesL(G2mL(i), :)) = AvuycnG(meshesS(G2mS(i), :), meshesL(G2mL(i), :)) + len*reshape(P2P2yP2nG(v1(1), v2(1), v3(1), v1(2), v2(2), v3(2))*c(meshesL(G2mL(i), :)), [6, 6]);

end
AvSuS = sparse(AvSuS);
AvLuL = sparse(AvLuL);
AvSuL = sparse(AvSuL);
AvsusG = sparse(AvsusG);
AvuxcnG = sparse(AvuxcnG);
AvuycnG = sparse(AvuycnG);
Avn1S = sparse(Avn1S);
Avn2S = sparse(Avn2S);
Avn1L = sparse(Avn1L);
Avn2L = sparse(Avn2L);
end