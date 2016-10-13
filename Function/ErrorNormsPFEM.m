function FEM = ErrorNormsPFEM(FEM)
% FEM = ErrorNormsPFEM(FEM)
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
enorm = 0;
dnorm1 = 0;
enorm1 = 0;
for iel = 1 : FEM.Mesh.Connect.NumElem
    nod = FEM.Mesh.Elements{iel};  % extract nodes for (iel)-th element
    Co = FEM.Mesh.Nodes(nod, :);
    nnel = length(nod);
    index = assembly(nod, nnel, FEM.Mesh.Connect.Dof);
    disElem = [FEM.U(index(1 : 2 : 2*nnel), 1), FEM.U(index(2 : 2 : 2*nnel), 1)];
    edisp = FEM.U(index);
    [CoISO, ConnectISO] = PolyTrnglt(nnel, [0, 0]);
    for isub = 1 : nnel
        CoSubTri = CoISO(ConnectISO(isub, :), :);
        [W, P] = TriangularQuad(CoSubTri, 2);
        for igaus = 1 : size(W, 1)
            [Ng, dNdxi] = Wachspress(CoISO(1:nnel, :), P(igaus, :));
            % Jacobian matrix on the Gauss point
            Jac = dNdxi'*Co;
            dNdxy = Jac \ dNdxi';
            Be = Bmat2D(nnel, dNdxy(1, :), dNdxy(2, :));
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
    end
    clear edisp
end
FEM.ErrorNorms.Displacement = sqrt(dnorm) / sqrt(dnorm1);
FEM.ErrorNorms.Energy = sqrt(0.5*enorm) / sqrt(0.5*enorm1);
end
