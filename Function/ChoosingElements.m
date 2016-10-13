function FEM = ChoosingElements(FEM, option)
% FEM = ChoosingElements(FEM, option)
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

%------------------------------------------%
%          Choose type of element          % 
%------------------------------------------%
if option == 1
    FEM.TypeElement.Name = 'T3';                
    FEM.TypeElement.Type = 1;        % type of element: [1]: three nodes   [0]: four nodes
    FEM.TypeElement.NumNodeElem = 3; % number of element nodes
elseif option == 2
    FEM.TypeElement.Name = 'Q4';                 
    FEM.TypeElement.Type = 0;
    FEM.TypeElement.NumNodeElem = 4;   
elseif option == 3
    FEM.TypeElement.Name = 'Polygon_W';               
    FEM.TypeElement.Type = 2; 
    FEM.TypeElement.NumNodeElem = 'NoLimit';      
end
end


