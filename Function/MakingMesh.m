function FEM = MakingMesh(FEM, iter)
% FEM = MakingMesh(FEM, iter)
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

%=============================  Mesh ======================================
FEM = BeamDomain(FEM, iter);
%==========================================================================
% Input data for nodal connectivity of each element                       
%==========================================================================
FEM.Mesh.Connect.Dof = 2;
if isnumeric(FEM.Mesh.Elements)
    FEM.Mesh.Connect.NumElem = size(FEM.Mesh.Elements,1);                     % number of elements
    FEM.Mesh.Connect.NumNode = size(FEM.Mesh.Nodes, 1);                        % total number of nodes in system
    FEM.Mesh.Connect.SysDof = FEM.Mesh.Connect.NumNode*FEM.Mesh.Connect.Dof;  % total system dofs
else
    FEM.Mesh.Connect.NumElem = length(FEM.Mesh.Elements);  
    FEM.Mesh.Connect.NumNode = length(FEM.Mesh.Nodes); 
    FEM.Mesh.Connect.SysDof = FEM.Mesh.Connect.NumNode*FEM.Mesh.Connect.Dof; 
end
end




