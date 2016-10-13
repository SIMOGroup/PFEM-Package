function FEM = BeamDomain(FEM, iter)
% FEM = BeamDomain(FEM, iter)
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

FEM.lx = 48;
FEM.ly = 12;
if FEM.TypeElement.Type == 2
    if iter == 1
        load dataMeshBeam1 
    elseif iter == 2
        load dataMeshBeam2
    elseif iter == 3
        load dataMeshBeam3 
    end
    FEM.Mesh.Nodes = gcoord; 
    FEM.Mesh.Elements = ele_nods;
else
    if iter == 1
        nx = 16;
        ny = 4;
    elseif iter == 2
        nx = 32;
        ny = 8;
    elseif iter == 3
        nx = 65;
        ny = 16;
    end
    [FEM.Mesh.Nodes, FEM.Mesh.Elements] = rectangleMesh(FEM.TypeElement.Type, FEM.TypeElement.NumNodeElem, 0, 48, -6, 6, nx, ny);
    Temp = cell(size(FEM.Mesh.Elements, 1), 1);
    for i = 1 : size(FEM.Mesh.Elements, 1)
        Temp{i} = FEM.Mesh.Elements(i, :);
    end
    FEM.Mesh.Elements = Temp;
end
end


