clear
clc

%%
NbSpec = 120;
Mu = randn(1,NbSpec);
STD = abs(randn(1,NbSpec));
Sigma_noise = abs(randn(NbSpec));
for iSpec=1:NbSpec
    Sigma(iSpec,iSpec) = STD(iSpec);
end
Data = mvnrnd(Mu, Sigma, 20);

%%
Var = {};
Var(1).level = 1:5;
Var(1).name = {'A', 'B', 'C', 'D', 'E'};
Var(2).level = 1:4;
Var(2).name = {'A', 'B', 'C', 'D'};
Var(3).level = 1:3;
Var(3).name = {'A', 'B', 'C'};
Var(4).level = 1:2;
Var(4).name = {'A', 'B'};

Sets = {Var(:).level};
[a, b, c, d] = ndgrid(Sets{:});
SpecsDef = [a(:), b(:), c(:), d(:)];

SpecMatrix = [];
SpecMatLabel = {};
for iVar = 1:numel(Var)
    for iLevel = 1:numel(Var(iVar).level)
        SpecMatrix(end+1,:) = ones(1,NbSpec); %#ok<SAGROW>
        Idx = find(SpecsDef(:,iVar)==Var(iVar).level(iLevel));
        SpecMatrix(end,Idx) = 0;
        
        SpecMatLabel{end+1} = [num2str(iVar) ' - ' Var(iVar).name{iLevel}];
    end
end



%%
Meam_emp = mean(Data);
STD_emp = std(Data);

H = ttest(Data,0);

[Meam_emp,I] = sort(Meam_emp);
STD_emp = STD_emp(I);
H = logical(H(I));
SpecMatrix = SpecMatrix(:,I);

SpecMatrix = repmat(SpecMatrix, [1,1,3]);
SpecMatrix(:,H,1) = 1;

%%
close all
figure('name', 'specification curve')

subplot(2,1,1)
hold on
for iSpec=1:NbSpec
    g(iSpec) = errorbar(iSpec,Meam_emp(iSpec),STD_emp(iSpec), ' .k', ...
        'linewidth', 1, 'markersize', 10);
end
plot([1 NbSpec], [0 0], '--k')
xlabel('specification number')

set(g(H), 'color' , 'r')


subplot(2,1,2)
% colormap gray
image(SpecMatrix)
xlabel('specification number')
axis([0 size(SpecMatrix,2)+1 0 size(SpecMatrix,1)+1])
set(gca, 'ytick', 1:size(SpecMatrix,1),'yticklabel', SpecMatLabel)
