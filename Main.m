%==========================================================================
% A simple educational matlab code for Polygonal Finite Element Method (PFEM)

%----------
% References:
% [1]  H. Nguyen-Xuan, Son Nguyen-Hoang, T. Rabczuk, K. Hackl
%      A polytree-based adaptive approach to limit analysis of cracked structures
%      Computer Methods in Applied Mechanics and Engineering, in press, 2016

% [2]  Son Nguyen-Hoang, H. Nguyen-Xuan
%      A polytree-base adaptive polygonal finite element method for topology optimization
%      International Journal for Numerical Methods in Engineering, in press, 2016

%----------
% Coded by:
% [1]  Son Nguyen-Hoang
%      Emails: mrnguyenhoangson@gmail.com

% [2]: Khai Chau-Nguyen
%      Emails: chaunguyenkhai@gmail.com

% [3]: Hung Nguyen-Xuan
%      Emails: h.nguyenxuan@gmail.com

%==========================================================================
% Problem: A rectangular cantilever beam (Plane stress, Poisson = 0.3)
% Ref: Example 4th in "Int. J. Numer. Meth. Engng 2008; 74 :175–208"

clc
clear
close all
format long g
addpath Function
% Choose elements
FEM.OptionElement = [1, 2, 3];
% [1]:T3
% [2]:Q4
% [3]:Polygon with "Wachspress shape fucntion"

FEM.Results = zeros(3, 3, length(FEM.OptionElement));
for ielement = 1 : length(FEM.OptionElement) % loop elements
    %------------------------
    FEM = ChoosingElements(FEM, FEM.OptionElement(ielement));
    % Materials
    FEM.E = 2.1e4; FEM.nu = 1 / 3; FEM.Force = -1000;
    for iter = 1 : 3 % loop meshes
        FEM = MakingMesh(FEM, iter);
        FEM = MakingBoundary(FEM);
        % Making force vector
        FEM = ForceVector(FEM);
        % Making stiffness matrix
        switch FEM.TypeElement.Name
            case {'T3','Q4'}
                FEM = StiffnessMatrixFEM(FEM);
            case {'Polygon_W'}
                FEM = StiffnessMatrixPFEM(FEM);
        end
        % Solving problem
        FEM = ApplyingBCCondition(FEM);
        FEM.U = FEM.K \ FEM.F;
        % Error norms
        switch FEM.TypeElement.Name
            case {'T3','Q4'}
                FEM = ErrorNormsFEM(FEM);
            case {'Polygon_W'}
                FEM = ErrorNormsPFEM(FEM);
        end
        FEM.Results(iter, 1, ielement) = FEM.ErrorNorms.Displacement;
        FEM.Results(iter, 2, ielement) = FEM.ErrorNorms.Energy;
        FEM.Results(iter, 3, ielement) = 1 / sqrt(FEM.Mesh.Connect.SysDof);
    end
    PlotPolyMesh(FEM.Mesh.Nodes, FEM.Mesh.Elements, 0)
end

figure
% Plot results
line_type = ['b-<','g-o','r-d'];

for i = 1 : length(FEM.OptionElement)
    plot(log10(FEM.Results(:, 3, i)), log10(FEM.Results(:, 1, i)), line_type(1, (i*3 - 2) : i*3));
    hold on
end
ylabel('Relative error in displacement norm')
xlabel('1/sqrt (d.o.f)');
legend('T3','Q4','Poly')

figure
for i = 1 : length(FEM.OptionElement)
    plot(log10(FEM.Results(:, 3, i)), log10(FEM.Results(:, 2, i)), line_type(1, (i*3 - 2) : i*3));
    hold on
end
ylabel('Relative error in energy norm')
xlabel('1/sqrt (d.o.f)');
legend('T3','Q4','Poly')


