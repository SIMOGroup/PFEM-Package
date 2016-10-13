function PlotPolyMesh(Node, Element, flag)
% PlotPolyMesh(Node, Element, flag)
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

NElem = length(Element);
if isnumeric(Element)
    for iel = 1 : NElem
        Element1{iel} = Element(iel, :);
    end
    Element = Element1;
end
Element = Element(1 : NElem)'; % Only plot the first block
MaxNVer = max(cellfun(@numel, Element)); % Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1, MaxNVer - numel(E))]; % Pad cells with NaN
ElemMat = cellfun(PadWNaN, Element, 'UniformOutput', false);
ElemMat = vertcat(ElemMat{:}); % Create padded element matrix
figure
axis equal
hh = patch('Faces', ElemMat, 'Vertices', Node, 'FaceColor','none'); 
set(hh, 'edgecolor', 'b', 'linewidth', 1);

if flag == 1
    % Plot node numbers
    for i = 1 : length(Node)
        h1 = text(Node(i, 1), Node(i, 2), num2str(i));
        set(h1, 'fontsize', 12, 'color', 'k');
        hold on
    end
    % Plot element numbers
    for e = 1 : length(Element)
        xc = mean(Node(Element{e}, :));
        h2 = text(xc(1), xc(2), sprintf('%d', e));
        set(h2, 'fontsize', 12, 'color', 'r');
        hold on
    end
end