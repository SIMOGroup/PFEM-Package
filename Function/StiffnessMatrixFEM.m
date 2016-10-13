function FEM = StiffnessMatrixFEM(FEM)
% FEM = StiffnessMatrixFEM(FEM)
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
nnel = FEM.TypeElement.NumNodeElem;
FEM.K = sparse(FEM.Mesh.Connect.SysDof, FEM.Mesh.Connect.SysDof); 
[P, W] = Quadrature(FEM.TypeElement.Type, 2);
[Ng, dNdxi, dNdeta] = ShapeFunc(FEM.TypeElement.Type, nnel, P);
for iel = 1 : FEM.Mesh.Connect.NumElem
    nod = FEM.Mesh.Elements{iel};
    Co = FEM.Mesh.Nodes(nod, :);
    edof = FEM.Mesh.Connect.Dof*nnel;
    Ke = sparse(edof, edof);
    for igaus = 1 : length(W)
        % Jacobian matrix on the Gauss point
        Jac = [dNdxi(igaus, :); dNdeta(igaus, :)]*Co;
        % Shape functions' derivatives in global coordinates
        dNdxy = Jac \ [dNdxi(igaus, :); dNdeta(igaus, :)];
        %dNdxy = Jac\[dNdxi_igaus;dNdeta_igaus];
        Be = Bmat2D(nnel, dNdxy(1, :), dNdxy(2, :)); % element strain-disp matrix  
        Ke = Ke + Be'*Dmatrix*Be*W(igaus)*det(Jac);% element stiffness matrix
    end
    % Assemble element stiffness matrices and force vectors
    index = assembly(nod, nnel, FEM.Mesh.Connect.Dof);
    FEM.K(index, index) = FEM.K(index, index) + Ke;
end
end