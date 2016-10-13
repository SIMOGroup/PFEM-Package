function FEM = ErrorNormsFEM(FEM)
% FEM = ErrorNormsFEM(FEM)
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
dnorm = 0;
enorm=0;
dnorm1 = 0;
enorm1 = 0;
nnel = FEM.TypeElement.NumNodeElem;
[P, W] = Quadrature(FEM.TypeElement.Type, 2);
[N, dNdxi, dNdeta] = ShapeFunc(FEM.TypeElement.Type, nnel, P);
for iel = 1 : FEM.Mesh.Connect.NumElem
    nod = FEM.Mesh.Elements{iel};
    Co = FEM.Mesh.Nodes(nod, :);
    index = assembly(nod, nnel, FEM.Mesh.Connect.Dof);
    disElem = [FEM.U(index(1 : 2 : 2*nnel), 1), FEM.U(index(2 : 2 : 2*nnel), 1)];
    edisp = FEM.U(index);
    for igaus = 1 : length(W)
        % Shape functions on Gauss point igaus
        Ng = N(igaus, :)';
        Jac = [dNdxi(igaus, :); dNdeta(igaus, :)]*Co;
        % Shape functions' derivatives in global coordinates
        dNdxy = Jac \ [dNdxi(igaus, :); dNdeta(igaus, :)];
        Be = Bmat2D(nnel, dNdxy(1, :), dNdxy(2, :));        % element strain-disp matrix
        % Exact values of displacements & stresses at Gauss integration point
        X0 = Co'*Ng; 
        [u0, s0, e0] = AnalyticalSolution(X0(1), X0(2), FEM.Force, FEM.E, FEM.nu, FEM.ly, FEM.lx);
        % Displacements & stresses at Gauss integration point
        uh = disElem'*Ng; 
        eh = Be*edisp;
        % Error norms
        dnorm = dnorm + (uh - u0)'*(uh - u0)*W(igaus)*det(Jac);
        enorm = enorm + (eh - e0)'*(Dmatrix)*(eh - e0)*W(igaus)*det(Jac);
        
        dnorm1 = dnorm1 + (u0)'*(u0)*W(igaus)*det(Jac);
        enorm1 = enorm1 + (e0)'*(Dmatrix)*(e0)*W(igaus)*det(Jac);
    end
    clear edisp
end
FEM.ErrorNorms.Displacement = sqrt(dnorm) / sqrt(dnorm1);
FEM.ErrorNorms.Energy = sqrt(0.5*enorm) / sqrt(0.5*enorm1);
end