function FEM = StiffnessMatrixPFEM(FEM)
% FEM = StiffnessMatrixPFEM(FEM)
%--------------------------------------------------------------------------
%{
Copyright (C) <2016>  <Son Nguyen-Hoang, Khai Chau-Nguyen, Hung Nguyen-Xuan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

Dmatrix = FEM.E / (1 - FEM.nu*FEM.nu)*[1, FEM.nu, 0; FEM.nu, 1, 0; 0, 0, (1 - FEM.nu) / 2];
FEM.K = sparse(FEM.Mesh.Connect.SysDof, FEM.Mesh.Connect.SysDof);
for iel = 1 : FEM.Mesh.Connect.NumElem
    nod = FEM.Mesh.Elements{iel};
    Co = FEM.Mesh.Nodes(nod, :);
    nnel = length(nod);
    edof = FEM.Mesh.Connect.Dof*nnel;
    Ke = sparse(edof, edof);
    [CoISO, ConnectISO] = PolyTrnglt(nnel, [0, 0]);
    for isub = 1 : nnel
        CoSubTri = CoISO(ConnectISO(isub, :), :);
        [W, P] = TriangularQuad(CoSubTri, 2);
        for igaus = 1 : length(W)
            [Ng, dNdxi] = Wachspress(CoISO(1 : nnel, :), P(igaus, :));
            % Jacobian matrix on the Gauss point
            Jac = dNdxi' * Co;
            % Shape functions' derivatives in global coordinates
            dNdxy = Jac \ dNdxi';
            Be = Bmat2D(nnel, dNdxy(1, :), dNdxy(2, :));        % element strain-disp matrix
            Ke = Ke + Be'*Dmatrix*Be*W(igaus)*det(Jac); % element stiffness matrix
        end
    end
    % Assemble element stiffness matrices and force vectors
    index = assembly(nod, nnel, FEM.Mesh.Connect.Dof);
    FEM.K(index, index) = FEM.K(index, index) + Ke;
end
end