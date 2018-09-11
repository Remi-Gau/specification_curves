clear
clc

%% Generate data
NbSpec = 120;
Mu = randn(1,NbSpec);
STD = abs(randn(1,NbSpec));
Sigma_noise = abs(randn(NbSpec));
for iSpec=1:NbSpec
    Sigma(iSpec,iSpec) = STD(iSpec);
end
Data = mvnrnd(Mu, Sigma, 20);

%% Get all the possible specifications 
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
clear sets a b c d

%% Creates a matrix to vizualize which level of each variable was used for each specification
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

%% Defines label for multiverse analysis
CombinedVar = {[1 4];[2 3]}; % how the variables are combined to form a 2D matrix of values

X = prod(cellfun(@max,{Var(CombinedVar{1}).level}));
Sets = {Var(CombinedVar{1}).level};
[a, b] = ndgrid(Sets{:}); 
X_label_list = [a(:), b(:)];
clear sets a b

for iX = 1:X
    ToPrint = [];
    for iVar = 1:size(X_label_list,2)
        ToPrint = [ToPrint, CombinedVar{1}(iVar), X_label_list(iX,iVar)];  %#ok<*AGROW>
    end
    X_labels{iX,1} = sprintf('Var%i-Lvl%i  ', ToPrint);
end
X_labels


Y = prod(cellfun(@max,{Var(CombinedVar{2}).level}));
Sets = {Var(CombinedVar{2}).level};
[a, b] = ndgrid(Sets{:}); 
Y_label_list = [a(:), b(:)];
clear sets a b

for iY = 1:Y
    ToPrint = [];
    for iVar = 1:size(Y_label_list,2)
        ToPrint = [ToPrint, CombinedVar{1}(iVar), Y_label_list(iY,iVar)];
    end
    Y_labels{iY,1} = sprintf('Var%i-Lvl%i  ', ToPrint);
end
Y_labels

MultiverseMatrix= [];
for iX=1:X
    CurrentSpec = [];
    CurrentSpec(CombinedVar{1}) = X_label_list(iX,:); %#ok<*SAGROW>
    for iY=1:Y
        CurrentSpec(CombinedVar{2}) = Y_label_list(iY,:);
        MultiverseMatrix(iX,iY,:) = CurrentSpec;
    end
end


%% Compute values to plot
Meam_emp = mean(Data);
STD_emp = std(Data);
[H, p] = ttest(Data,0);

PValMat = nan(X,Y);
for iSpec = 1:NbSpec
    CurrentSpec = reshape(SpecsDef(iSpec,:),[1,1,size(SpecsDef,2)]);
    CurrentSpec = repmat(CurrentSpec, [size(MultiverseMatrix,1),size(MultiverseMatrix,2),1]);
    [I,J] = find(all(MultiverseMatrix==CurrentSpec,3));
    PValMat(I,J) = p(iSpec);
    clear I J
end

% sort the speficications
[Meam_emp,I] = sort(Meam_emp);
STD_emp = STD_emp(I);
H = logical(H(I));
SpecMatrix = SpecMatrix(:,I);

% Make red the significant specifications (leave the other ones black)
SpecMatrix = repmat(SpecMatrix, [1,1,3]);
SpecMatrix(:,H,1) = 1;


%% print multiverse analysis  
close all

figure('name', 'multiverse analysis')

subplot(2,1,1)
hold on
hist(p, 100)
ax = axis;
axis([0 1 0 (ax(4))])

subplot(2,1,2)
hold on
imagesc(PValMat)
for iX=1:X
    for iY=1:Y
        t = text(iY-.2, iX, sprintf('%.3f',PValMat(iX,iY)) );
        set(t, 'color', 'b')
    end
end
colormap gray
axis tight
set(gca, ...
    'xtick', 1:Y, 'xticklabel', Y_labels, ...
    'ytick', 1:X, 'yticklabel', X_labels)
colorbar

%% print specification curve
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




