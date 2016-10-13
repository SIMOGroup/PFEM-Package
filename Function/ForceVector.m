function FEM = ForceVector(FEM)
% FEM = ForceVector(FEM)
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

FEM.F = zeros(FEM.Mesh.Connect.SysDof, 1);
I = FEM.ly ^ 3 / 12;
% Define quadrature rule
W = [1; 1];
Q = [1 / sqrt(3); -1 / sqrt(3)];

for e = 1 : size(FEM.BC.NeumannBC, 1)                          % loop over the elements in the right edge
    sctr = FEM.BC.NeumannBC(e, :);                           % scatter vector for the element
    sctry = FEM.Mesh.Connect.Dof*sctr;                      % y scatter vector
    for q = 1 : size(W, 1)                                     % quadrature loop
        xi = Q(q, :);                                        % quadrature point
        wt = W(q);                                          % quadrature weight
        N = ([1 - xi, 1 + xi] / 2)';                               % element shape functions
        dNdxi = [-1; 1] / 2;
        J0 = dNdxi'*FEM.Mesh.Nodes(sctr, :);                 % Jacobian
        detJ0 = norm(J0);                                   % determiniat of Jacobian
        yPt = N'*FEM.Mesh.Nodes(sctr, 2);                    % y coordinate at quadrature point
        ff = -FEM.Force*(FEM.ly^2 / 4 - yPt^2) / (2*I);           % y traction at quadrature point
        FEM.F(sctry) = FEM.F(sctry) + N*ff*detJ0*wt;          % scatter force into global force vector
    end
end
end
