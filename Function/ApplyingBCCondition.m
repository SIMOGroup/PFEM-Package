function FEM = ApplyingBCCondition(FEM)
% FEM = ApplyingBCCondition(FEM)
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

n = size(FEM.BC.BCDof, 2);
sdof = size(FEM.K);
for i = 1 : n
    c = FEM.BC.BCDof(i);
    for j = 1 : sdof
        FEM.K(c, j) = 0;
    end
    FEM.K(c, c) = 1;
    FEM.F(c, 1) = FEM.BC.BCValue(i);
end
end

