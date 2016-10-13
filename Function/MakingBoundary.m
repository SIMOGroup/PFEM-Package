function FEM = MakingBoundary(FEM)
% FEM = MakingBoundary(FEM)
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

tol = 1e-3;
FEM.BC.BCDof = [];
FEM.BC.BCValue = [];
RigthBC = find(abs(FEM.Mesh.Nodes(:, 1) - FEM.lx) < tol);
for i = 1 : length(RigthBC) - 1
    for j = i + 1 : length(RigthBC)
        if FEM.Mesh.Nodes(RigthBC(i), 2) > FEM.Mesh.Nodes(RigthBC(j), 2)
            temp = RigthBC(i);
            RigthBC(i) = RigthBC(j);
            RigthBC(j) = temp;
        end
    end
end
FEM.BC.NeumannBC = [RigthBC(1 : end-1), RigthBC(2 : end)];
DirichletBC = find(abs(FEM.Mesh.Nodes(:, 1) - 0) <tol );
for i = 1 : size(DirichletBC)
    FEM.BC.BCDof = [FEM.BC.BCDof (DirichletBC(i)*2-1) (DirichletBC(i)*2)];
    X0 = FEM.Mesh.Nodes(DirichletBC(i),:);
    [u0, s0, e0] = AnalyticalSolution(X0(1), X0(2), FEM.Force, FEM.E, FEM.nu, FEM.ly, FEM.lx);
    FEM.BC.BCValue=[FEM.BC.BCValue u0(1) u0(2)];
end
end


