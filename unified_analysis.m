% new script; combine all analysis--- reference other scripts (functions), organize
% correctly and simply
% create all figures based on the same loaded mutation data


cd('C:\Users\rmathis\Dropbox (MIT)\mutations and evolution');
warning('off','MATLAB:table:ModifiedVarnames');

%paths
path_cur = 'C:\Users\rmathis\Dropbox (MIT)\mutations and evolution\';
path_mut='C:\Users\rmathis\Dropbox (MIT)\mutations and evolution\maf_files\';
path_rna = 'C:\Users\rmathis\Dropbox (MIT)\mutations and evolution\rna_files\';
path_NC = 'C:\Users\rmathis\Dropbox (MIT)\mutations and evolution\multivariate analyses\replication timing\';
path_msigdb = 'C:\Users\rmathis\Dropbox (MIT)\mutations and evolution\molsigdb\';

%figure properties
FigProps.InvertHardcopy = 'off';
FigProps.Units = 'inches';
FigProps.PaperPositionMode = 'Auto';
FigProps.Renderer='painters';
FigProps.Color = 'none';

%axis properties.
AxProps.Color = 'none';
AxProps.FontName = 'Arial';
AxProps.FontSize = 10;
AxProps.TickDir = 'out';
AxProps.TickLength = [0.04,0.03];
AxProps.LineWidth = 1;
AxProps.Box = 'off';
AxProps.TitleFontSizeMultiplier=1;
AxProps.LabelFontSizeMultiplier = 1; 

% which tumor types to load in? 
ttypes = {'SKCM','LUAD','COADREAD','UCEC','BRCA','HNSC','PRAD','STAD','GBMLGG','LIHC','BLCA'};
%melanomas, lung adenocarcinomas, colorectal adenocarcinomas, liver hepatocellular carcinomas, gliomas, and breast invasive carcinomas
%ttypes = {'SKCM','LUAD','COADREAD','LIHC','GBMLGG','BRCA'};


%% Loading up mutation data
% we want the missense mutations in each tumor type and the missense
% mutations across all tumor types, dropping recurrent mutations

mut_temp = cell(size(ttypes));
mut_single_tumors = cell(size(ttypes));
tumor_number = zeros(size(ttypes));
test=datastore(horzcat(path_mut,'*.txt'),'delimiter','\t','whitespace',' \b');

for i = 1:length(test.Files)
    %pull out tumor type ID
    nm = strsplit(test.Files{i},'\');
    fname = nm{end};
    nm = strsplit(nm{end},'.');
    nm = strsplit(nm{1},'-');
    nm = nm{1};
    if ~any(strcmpi(ttypes,nm)),
        continue
    end
    t = find(strcmpi(ttypes,nm));
    mut = GetMutationData(path_mut,fname);%read missense and tabulate
    tumor_number(t) = size(mut,2);
    mut_single_tumors{t} = mut;
    mut_temp{t} = table(mut.Properties.RowNames,sum(mut{:,:},2),'VariableNames',{'genenames',horzcat('mis_',nm)},'RowNames',mut.Properties.RowNames);
   
end

%combine data from all tumors
for i=1:length(mut_temp)
    switch i
        case 1
            mut = mut_temp{i};
        otherwise
            mut = outerjoin(mut,mut_temp{i},'key','genenames','MergeKeys',1);
    end
end
mut.Properties.RowNames = mut.genenames;
mut.genenames = [];

%replace NaNs with 0s
for i=1:width(mut)
    mut{isnan(mut{:,i}),i} = 0;
end
clear nm t fname  test i 

%% Get RNA data
expr_threshold = 3;
stringency = 0.95;    
test = datastore(horzcat(path_rna,'*.txt'),'delimiter','\t','whitespace',' \b');
is_expr = false(size(mut));
is_nexpr = false(size(mut));
expr_lvl = zeros(size(mut));
expr_single_tumors = cell(size(mut_single_tumors));
for i = 1:length(test.Files)
    %pull out tumor type ID
    nm = strsplit(test.Files{i},'\');
    fname = nm{end};
    nm = strsplit(nm{end},'.');
    nm = strsplit(nm{1},'-');
    nm = nm{1};
    if ~any(strcmpi(ttypes,nm)),
        continue
    end
    t = find(strcmpi(ttypes,nm));
    expr_data = GetRNAData(path_rna,fname,mut.Properties.RowNames);
    expr_single_tumors{i} =  GetRNAData(path_rna,fname,mut.Properties.RowNames,mut_single_tumors{i}.Properties.VariableNames);
%    expr_data1 = expr_data{:,~all(isnan(expr_data{:,:}),1)};
    is_expr(:,t) = sum(expr_data{:,:} >expr_threshold,2) > stringency*size(expr_data,2);
    is_nexpr(:,t) = sum(expr_data{:,:} >expr_threshold,2) < (1-stringency)*size(expr_data,2);
    expr_lvl(:,t) = mean(expr_data{:,:},2);
end
is_expr = array2table(is_expr,'RowNames',mut.Properties.RowNames,'VariableNames',ttypes);
is_nexpr = array2table(is_nexpr,'RowNames',mut.Properties.RowNames,'VariableNames',ttypes);
expr_lvl = array2table(expr_lvl,'RowNames',mut.Properties.RowNames,'VariableNames',ttypes);
clear expr_threshold stringency test nm fname t i expr_data

%% Get other data
NoncodingMutRate = GetNoncodingMutRate(mut.Properties.RowNames, path_NC);
length_AAs = GetGenelengths(mut.Properties.RowNames,path_cur);

sets = {'c5.all.v5.1.symbols.gmt',...
    'c2.cp.reactome.v5.1.symbols.gmt',...
    'c2.cp.kegg.v5.1.symbols.gmt'...
    'c2.cp.biocarta.v5.1.symbols.gmt'...
    'c2.cp.v5.1.symbols.gmt','h.all.v5.1.symbols.gmt'};

for i =1:length(sets);
    [SetNames_t,MolSigDB_t,NumSets_t] = GrabSets(path_msigdb,sets{i}); %grab molsigdb data
    switch i 
        case 1
            SetNames = SetNames_t;
            MolSigDB = MolSigDB_t;
            NumSets = NumSets_t;
        otherwise
            SetNames = vertcat(SetNames,SetNames_t);
            MolSigDB = vertcat(MolSigDB,MolSigDB_t);
            NumSets = vertcat(NumSets,NumSets_t);
    end
end
[~,ia,~] = unique(SetNames); %drop repeat sets
SetNames = SetNames(ia);
MolSigDB = MolSigDB(ia);
clear sets SetNames_t MolSigDB_t NumSets_t NumSets

%filter sets for genes we have data for;
%and turn sets into a logical;
gnms = mut.Properties.RowNames;
in_set_mat = false(length(SetNames),height(mut));
for i = 1:length(SetNames) % for each set
    gs = MolSigDB{i}; %grab set
    in_set_mat(i,:)  = ismember(gnms,gs);
end
clear gnms gs i 

%load in published datasets from pooled CRISPR screens;
%data from Hart et al (DOI: 10.1016/j.cell.2015.11.015)
hart_et_al = readtable(horzcat(path_cur,'essentiality_from_other_papers\hart_et_al_T1.xlsx'),'ReadRowNames',1);

%data from tzelepis et al (DOI: 10.1016/j.celrep.2016.09.079)
tzelepis_et_al = readtable(horzcat(path_cur,'essentiality_from_other_papers\tzelepis_et_al.xlsx'),'Sheet','Dropout at FDR10%','ReadRowNames',1);
tz_genes = table2array(readtable(horzcat(path_cur,'essentiality_from_other_papers\tzelepis_et_al.xlsx'),'Sheet','CRISPR dropout summary','ReadRowNames',0,'Range','A2:A18011'));
tz_data1 = tzelepis_et_al(:,'x_All');
ne = tz_genes(~ismember(tz_genes,tz_data1.Properties.RowNames)); %not essential genes
tz_data = vertcat(tz_data1,table(zeros(size(ne)),'RowNames',ne,'VariableNames',tz_data1.Properties.VariableNames)); %append nonessential genes examined
clear tzelepis_et_al tz_genes tz_data1 ne

%data from Wang et al
wang_et_al_i = readtable(horzcat(path_cur,'aac7041_SM_Table_S3.txt'),'delimiter','\t','ReadRowNames',true);
wang_et_al = wang_et_al_i(:,{'KBM7AdjustedP_value','K562AdjustedP_value','JiyoyeAdjustedP_value','RajiAdjustedP_value'});
w_data = array2table(sum(wang_et_al{:,:}<0.01,2),'RowNames',wang_et_al.Properties.RowNames,'VariableNames',{'x_All'});
w_data.avgP = geomean(wang_et_al{:,:},2);
clear wang_et_al_i wang_et_al



%%  Identifying genes under purifying selection in multiple tumor types 
ttypes_used = {'SKCM','LUAD','COADREAD','LIHC','GBMLGG','BRCA'}; 
keep = ~isnan(length_AAs{:,:}) & ~isnan(NoncodingMutRate{:,:}) & all(is_expr{:,ismember(ttypes,ttypes_used)},2); % only keep mutations from genes which we have all the other data for, and are expressed;
tmp = table();
tmp.genenms = mut.Properties.RowNames(keep); % gene names
tmp.observed_muts = sum(mut{keep,ismember(ttypes,ttypes_used)},2); % sum mutations in these tumor types
RelativeNCMutRate = NoncodingMutRate{keep,:} ./ mean(NoncodingMutRate{keep,:});
ExpectedNCMutRate = RelativeNCMutRate * sum(tmp.observed_muts) / sum(3*length_AAs{keep,:}); 
tmp.expected_muts = ExpectedNCMutRate .* (3 * length_AAs{keep,:}); % calculate expected muts
tmp.Pval = cdf('Poisson',tmp.observed_muts,tmp.expected_muts); % determine p value
tmp.obs_exp = tmp.observed_muts ./ tmp.expected_muts; % determine fold change
tmp.Q = mafdr(tmp.Pval,'BHFDR',1);
tmp.SigDep = tmp.Pval < 0.01 & tmp.obs_exp <.5; % significant depletion
tmp.SigDep2 = tmp.Q < 0.01;
tmp = sortrows(tmp,'Q');
if exist('NewFigureGeneration\TableS2_genes_alltumors.txt','file')>0
    delete('NewFigureGeneration\TableS2_genes_alltumors.txt');
end
writetable(tmp,'C:\Users\rmathis\Dropbox (MIT)\mutations and evolution\NewFigureGeneration\TableS2_genes_alltumors.txt','delimiter','\t');
clear ttypes_used keep RelativeNCMutRate ExpectedNCMutRate 

% Make figure 3A
never  = tmp.expected_muts < 10;
ranged = tmp.obs_exp;
ranged(ranged > 2^3) = 2^3;
ranged(ranged < 2^-3) = 2^-3;
[~,edges] = histcounts(log2(ranged(~never)),'BinWidth',0.2);
edges =2.^edges;
fig = figure();
set(fig,FigProps);
fig.Position(3:4) = [3.5 1.5];
hold on
h=histogram(ranged(~never),edges);
h.EdgeColor=[1 1 1];
h.FaceAlpha = 1;
h2=histogram(ranged(~never & tmp.SigDep),edges);
h2.EdgeColor=[1 1 1];
h2.FaceAlpha = 1;
ax=gca;
set(ax,AxProps);
ax.XScale='log';
ax.XTick = 2.^[-3 -2 -1 0 1 2 3];
Label = ax.XTick;
Label(ax.XTick<1) = 1./ax.XTick(ax.XTick<1);
ax.XTickLabel=Label;
xlabel('Fold Depletion | Fold Enrichment','FontName','Arial','FontSize',10)
ylabel('No. Genes','FontName','Arial','FontSize',10)
title('All tumors','FontName','Arial','FontSize',10,'FontWeight','Normal');
L = legend({'All genes','UnderSelection'},'Location','NE');
L.Box = 'off';
L.FontName = 'Arial';
L.FontSize = 10;

print('NewFigureGeneration\Fig3A_histogram_alltumors_depl_muts.svg','-dsvg');

clear ranged never edges tmp L ax fig h h2 Label 

%% Identifying gene sets under purifying selection in multiple tumor types
%keep sets with size >= 10 and <= 100
%also keep if have >50% expressed genes in each tumor type;
ttypes_used = {'SKCM','LUAD','COADREAD'};
keepG = ~isnan(length_AAs{:,:}) & ~isnan(NoncodingMutRate{:,:}); % only keep mutations from genes which we have all the other data for

L = sum(in_set_mat(:,keepG),2);
expr_t = double(in_set_mat(:,keepG)) * double(is_expr{keepG,ismember(ttypes,ttypes_used)}); % # of genes expressed in each set, in each tumor type
keep = L >=10 & L <= 400 ...
    & all(expr_t > .5 .* repmat(L,1,size(ttypes_used,2)),2);
in_set_mat1 = in_set_mat(keep,keepG);

observed_muts = sum(mut{keepG,ismember(ttypes,ttypes_used)},2); % sum mutations in these tumor types
RelativeNCMutRate = NoncodingMutRate{keepG,:} ./ nanmean(NoncodingMutRate{keepG,:});
ExpectedNCMutRate = RelativeNCMutRate * (sum(observed_muts) / sum(3*length_AAs{keepG,:})); 
expected_muts = ExpectedNCMutRate .* (3 * length_AAs{keepG,:}); %calculate expected mutations

tmp = table();
tmp.observed_muts = double(in_set_mat1) * observed_muts; %sum up observed mutations in these sets
tmp.expected_muts = double(in_set_mat1) * expected_muts; %sum up observed mutations in these sets
tmp.Diff = cdf('Poisson',tmp.observed_muts,tmp.expected_muts); %compute a pvalue
tmp.obs_exp = tmp.observed_muts ./ tmp.expected_muts; %also obs / expected
tmp.SetSizes = sum(in_set_mat1,2);
tmp.SetNames = SetNames(keep);

%Now find p values
rep = 1*10^4; %replicates
expr = double(is_expr{keepG & any(in_set_mat(keep,:),1)',ismember(ttypes,ttypes_used)}); %precompute expressed 
observed_m = observed_muts(any(in_set_mat1,1)); %prefilter observed muts
expected_m = expected_muts(any(in_set_mat1,1)); %prefilter expected muts
%Pull random gene sets of different lengths
setsizesLim = [min(sum(in_set_mat1,2)),max(sum(in_set_mat1,2))];
gns = 1:nnz(any(in_set_mat1,1)); %gene labels; all gns in sets after filtering
SetSizes = tmp.SetSizes;
[tmp.emp_pval,null_diff] = GeneSetRandomization(rep,observed_m,expected_m,setsizesLim,gns,tmp.Diff,expr,SetSizes);

%for very unlikely things, fit null_diff to estimate pvalues
ssizes = setsizesLim(1):setsizesLim(2);
tmp.est_pval = ones(size(tmp.emp_pval));
for i = 1:size(null_diff,2)
    nd = null_diff(:,i);
    p = 1/length(nd) : 1/length(nd) : 1 ;
    Q = quantile(nd,p);
    x = -log10(Q);
    y = -log10(p);
    mdl1 = fitlm(x,y,'linear','exclude',Q>0.01 | y>3.5);
    s = tmp.SetSizes == ssizes(i);
    est = predict(mdl1,-log10(tmp.Diff(s)));
    tmp.est_pval(s) = 10.^-est;
end

clear i nd p Q x y mdl1 s est

tmp.Q = mafdr(tmp.emp_pval,'BHFDR',1);
tmp = sortrows(tmp,'Q');
if exist('NewFigureGeneration\TableS3_genesets_alltumors.txt','file')>0
    delete('NewFigureGeneration\TableS3_genesets_alltumors.txt');
end
writetable(tmp,'C:\Users\rmathis\Dropbox (MIT)\mutations and evolution\NewFigureGeneration\TableS3_genesets_alltumors.txt'...
    ,'delimiter','\t','WriteRowNames',0);

% Create figure 3B
tmp_fig = tmp(tmp.Q<0.01,{'SetNames','est_pval','obs_exp'});
tmp_fig.DepletionOfMutations = 1-tmp_fig.obs_exp;
tmp_fig.obs_exp = [];
tmp_fig.Properties.VariableNames = {'SetName','Pvalue','DepletionOfMutations'};

%how many of the genes in each one of these sets have fewer mutations than expected?
alltumor_gene_input = readtable(horzcat(path_cur,'\NewFigureGeneration\TableS2_genes_alltumors.txt'),'ReadRowNames',1,'Delimiter','\t');
alltumor_input = tmp;
genesets = alltumor_input.SetNames;
%for each set that's under selection, pull out those genes
FracDepleted = zeros(size(genesets));
for s = 1:length(genesets)
    gns = MolSigDB{strcmp(SetNames,genesets{s})};
    gns = gns(ismember(gns,alltumor_gene_input.Properties.RowNames));%filter for genes we have obs_exp data for
    obs_exp = alltumor_gene_input{gns,'obs_exp'};
    FracDepleted(s) = nnz(obs_exp < 1) / length(obs_exp);
end
FracDepleted = array2table(FracDepleted,'RowNames',genesets,'VariableNames',{'FracDepleted'});
FracDepleted.SigDep = alltumor_input.Q<0.05 & alltumor_input.obs_exp<0.8; % genesets we found to be under purifying selection

tmp_fig.FracGenesWithReducedMutations = FracDepleted{tmp_fig.SetName,'FracDepleted'}; %now record the fraction of genes with reduced mutations
tmp_fig = sortrows(tmp_fig,'DepletionOfMutations','descend'); %sort by depletion of mutations
tmp_fig2 = tmp_fig(1:20,:);
%save as excel document
%if this file already exists delete it
if exist('NewFigureGeneration\Fig3B.xlsx','file')>0
    delete('NewFigureGeneration\Fig3B.xlsx');
end
writetable(tmp_fig2,'C:\Users\rmathis\Dropbox (MIT)\mutations and evolution\NewFigureGeneration\Fig3B.xlsx');

% Create Figure 3C
fig = figure();
set(fig,FigProps);
fig.Position(3:4) = [3 2];
hold on
[y,edgs] = histcounts(FracDepleted.FracDepleted(~FracDepleted.SigDep),'Normalization','probability','BinWidth',.07);
x = edgs(2:end)-diff(edgs)/2;
plot(x,y,'-o','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',4);
[y,edgs] = histcounts(FracDepleted.FracDepleted(FracDepleted.SigDep),'Normalization','probability','BinWidth',.07);
x = edgs(2:end)-diff(edgs)/2;
plot(x,y,'-o','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',4);
p = ranksum(FracDepleted.FracDepleted(FracDepleted.SigDep),FracDepleted.FracDepleted(~FracDepleted.SigDep));
if p < 1*10^-20;
    txt = 'p < 1E-20';
else
    txt = horzcat('p = ',num2str(p));
end
ax=gca;
set(ax,AxProps);
ax.XLim = [0 1];
ylabel('Gene set frequency','FontSize',10,'FontName','Arial')
xlabel({'Genes depleted for mutations'},'FontSize',10,'FontName','Arial');
L=legend({'No Selection','Under Selection'},'Box','off','Location','NW','FontSize',10,'FontName','Arial','Units','Normalized');
L.Position(1:2) = [0.15 0.75];
text(.05,.68,txt,'Units','Normalized','FontSize',10,'FontName','Arial');
print('NewFigureGeneration\Fig3C_histogram_set_genes_depl.svg','-dsvg');

clear gnms in_set_mat1 i gs ttypes_used keep RelativeNCMutRate ExpectedNCMutRate...
    tmp keepG expected_muts observed_muts expr gns gnsinset ssizes setsizes exprnum obs_r r exp_r...
    s s2 indi indi2 rnds rn fails works gns n observed_m expected_m null_diff G rep L ia expr_t SetSizes setsizesLim...
    tmp_fig2 tmp_fig gns obs_exp genesets alltumor_gene_input alltumor_input fig x y edgs p ax L txt FracDepleted



%% Essentiality analysis of genes under purifying selection across all tumors
%Load in the analysed data--- gene sets under purifying selection across
%tumors
alltumor_input = readtable(horzcat(path_cur,'\NewFigureGeneration\TableS3_genesets_alltumors.txt'),'ReadRowNames',0,'Delimiter','\t');
genesets_sel = alltumor_input.SetNames(alltumor_input.Q<0.05 & alltumor_input.obs_exp<0.8); % genesets we found to be under purifying selection
temp2 = MolSigDB(ismember(SetNames,genesets_sel));
genes_sel = unique(vertcat(temp2{:}));
setgenes = unique(vertcat(MolSigDB{:}));
genenms = mut.Properties.RowNames(~isnan(length_AAs{:,:}) & ~isnan(NoncodingMutRate{:,:}));
genes_sel = genes_sel(ismember(genes_sel,genenms)); %filter for genes analyzed
setgenes = setgenes(ismember(setgenes,genenms)); %filter for genes analyzed
alltumor_data = array2table(ismember(setgenes,genes_sel),'RowNames',setgenes,'VariableNames',{'PurifyingSelection'});
clear temp2 genenms


% Make figure 3D
genes = intersect(w_data.Properties.RowNames,intersect(tz_data.Properties.RowNames,intersect(alltumor_data.Properties.RowNames, hart_et_al.Properties.RowNames)));
t = table();
t.genes = genes;
t.Pur_alltumors = alltumor_data{genes,:};
t.Hart_et_al = hart_et_al{genes,'numTKOHits'};   %for each gene, number of cell lines in which it is deemed essential
t.Tzelepis_et_al = tz_data{genes,'x_All'};
t.Wang_et_al = w_data{genes,'x_All'};
t.Properties.RowNames = t.genes;
t.genes = [];
t.summed = sum(t{:,{'Hart_et_al','Wang_et_al','Tzelepis_et_al'}},2);
pur = t.Pur_alltumors;
ess_cls = t.summed ./ max(t.summed);
ess_cls_thr = ess_cls;
tr = .5;
ess_cls_thr(ess_cls_thr >tr) = tr;
fig = figure();
set(fig,FigProps);
fig.Position(3:4) = [3 2];
hold on
[y,edgs] = histcounts(ess_cls_thr(~pur),'Normalization','probability','BinWidth',0.1,'BinLimits',[-0.05,0.55]);
x = edgs(2:end)-diff(edgs)/2;
plot(x,y,'-o','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',4);
[y,edgs] = histcounts(ess_cls_thr(pur),'Normalization','probability','BinWidth',0.1,'BinLimits',[-0.05,0.55]);
x = edgs(2:end)-diff(edgs)/2;
plot(x,y,'-o','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',4);
ax=gca;
set(ax,AxProps);
ax.XTickLabels(strcmp(ax.XTickLabels,num2str(tr))) = {horzcat('>=',num2str(tr))};
ylabel('Gene frequency','FontSize',10,'FontName','Arial')
xlabel({'Essentiality score (CRISPR)'},'FontSize',10,'FontName','Arial');
p = ranksum(ess_cls(pur),ess_cls(~pur));
if p < 1*10^-20;
    txt = 'p < 1E-20';
else
    txt = horzcat('p = ',num2str(p));
end
L=legend({'No Selection','Under Selection'},'Box','off','Location','NW','FontSize',10,'FontName','Arial','Units','Normalized');
L.Position(1:2) = [0.45 0.75];
text(.55,.68,txt,'Units','Normalized','FontSize',10,'FontName','Arial');
print('NewFigureGeneration\Fig3D_histogram_set_genes_ess.svg','-dsvg');

clear X Y L label fig ax t i s AUC e linespec pred_scores 
clear genes t pur ess_cls_thr tr ess_cls x y edgs L ax fig txt...
    all_tumor_input genesets_sel setgenes genes_sel genenms alltumor_data genes t pur ess_cls_thr tr ess_cls x y edgs...
    ttypes_used g i2 ax1 ax2 b

%% Plot ROC curve for prediction of essentiality from purifying selection
% (Figure 3E)
alltumor_input = readtable(horzcat(path_cur,'\NewFigureGeneration\TableS3_genesets_alltumors.txt'),'ReadRowNames',0,'Delimiter','\t');
genesets_sel = alltumor_input.SetNames(alltumor_input.Q<0.05 & alltumor_input.obs_exp<0.8); % genesets we found to be under purifying selection
temp2 = MolSigDB(ismember(SetNames,genesets_sel));
genes_sel = unique(vertcat(temp2{:}));
setgenes = unique(vertcat(MolSigDB{:}));
genenms = mut.Properties.RowNames(~isnan(length_AAs{:,:}) & ~isnan(NoncodingMutRate{:,:}));
genes_sel = genes_sel(ismember(genes_sel,genenms)); %filter for genes analyzed
setgenes = setgenes(ismember(setgenes,genenms)); %filter for genes analyzed
alltumor_data = array2table(ismember(setgenes,genes_sel),'RowNames',setgenes,'VariableNames',{'PurifyingSelection'});
clear temp2 genenms

ttypes_used = {'SKCM','LUAD','COADREAD'}; % tumor types used in this analysis

genes = intersect(tz_data.Properties.RowNames,alltumor_data.Properties.RowNames);
t = table();
t.genes = genes;
t.Pur_alltumors = alltumor_data{genes,:};
t.Tzelepis_et_al = tz_data{genes,'x_All'} >= 5;
t.Properties.RowNames = t.genes;
t.Expression = mean(expr_lvl{genes,ismember(ttypes,ttypes_used)},2);

%for each gene grab the minimum P value of sets that gene is in.
tmp = sortrows(alltumor_input,'est_Pval','descend');
t.est_pval = ones(size(t.genes));
i2 = ismember(mut.Properties.RowNames,genes)'; 
for s = 1:height(tmp);
    i = strcmp(SetNames,tmp.SetNames(s)); %find genes in set
    g = mut.Properties.RowNames(in_set_mat(i,:) & i2); % genes in set
    t{g,'est_pval'} = tmp{s,'est_Pval'};
end
truelabels = t.Tzelepis_et_al;
fig = figure();
label = cell(2,1);
set(fig,FigProps);
fig.Position(3:4) = [3 2];
ax1 = axes('Position',[0.20 0.2520 0.7530 0.7230]);
ax1.ColorOrderIndex = 2.0;
ax2 = axes('Position',[0.7 0.3 0.28 0.28]);
AUC = zeros(3,1);
hold(ax1,'on');
hold(ax2,'on');
for e = 1:2
    switch e
        case 1
            pred_scores = 1-t.est_pval;
            L = 'Purifying selection';
            linespec = '-';
            fcolor = [0.8500 0.3250 0.0980];
            bline = '-';
        case 3
            pred_scores = 1-t.Wang_et_al;
            L = 'CRISPR-w';
            linespec = '-k';
            fcolor = 'k';
            bline = '-';    
        case 2  %seems less important; definitionally going to be x=y
            pred_scores = 1-t.est_pval(randperm(height(t)));
            L = 'Randomized scores';
            linespec = ':k';
            fcolor = 'none';
            bline = ':';
        case 4
            pred_scores = t.Expression;
            L = 'Expression';
            linespec = '-b';
            fcolor = 'b';
            bline = '-';
    end
    
    [X,Y,~,AUC(e)] = perfcurve(truelabels,pred_scores,1);
    plot(ax1,X,Y,linespec,'LineWidth',3,'MarkerSize',1);
    bar(ax2,e,AUC(e),'FaceColor',fcolor,'LineStyle',bline,'LineWidth',1.5);

    label{e} = L;
    %label{e} = horzcat(L,'; AUC = ',num2str(round(AUC(e),2)));
end
set(ax1,AxProps);
ylabel(ax1,'Specificity','FontName','Arial','FontSize',10);
xlabel(ax1,'1 - Sensitivity','FontName','Arial','FontSize',10);
legend(ax1,label,'Box','off','Location','NE','FontSize',10,'FontName','Arial','Units','Normalized');
%line(0:1,0:1,'LineStyle',':','Color','k','Parent',ax1)
fig.Color = [1 1 1];
ax2.YLim = [0.5 1];
ylabel(ax2,'AUC');
set(ax2,AxProps);
ax2.FontSize = 8;
ax2.XTick = [];

print('NewFigureGeneration\Fig3E_ROC.svg','-dsvg');

clear X Y L label fig ax t i s AUC e linespec pred_scores 
clear genes t pur ess_cls_thr tr ess_cls x y edgs L ax fig txt...
    all_tumor_input genesets_sel setgenes genes_sel genenms alltumor_data genes t pur ess_cls_thr tr ess_cls x y edgs...
    ttypes_used g i2 ax1 ax2 b bline truelabels fcolor

%% Identifying genes under tumor type-specific purifying selection
% build expected value based on other tumors
ttypes_used = ttypes;  %tumors to build expected?
ttypes_analyzed = {'SKCM','LUAD'};

% drop genes with > 2.5x more mutations than expected across tumors
keep = ~isnan(length_AAs{:,:});
RelativeNCMutRate = NoncodingMutRate{keep,:} ./ nanmean(NoncodingMutRate{keep,:});
ExpectedNCMutRate = RelativeNCMutRate * (sum(mut{keep,ismember(ttypes,ttypes_used)}(:)) / sum(3*length_AAs{keep,:})); 
exp_muts = ExpectedNCMutRate .* (3 * length_AAs{keep,:});
obs_muts = sum(mut{keep,ismember(ttypes,ttypes_used)},2);
obs_exp = obs_muts./exp_muts;
out_gr_exp = obs_exp > 2.5; 
gnms = mut.Properties.RowNames(keep);
out_gr_exp = ismember(mut.Properties.RowNames,gnms(out_gr_exp));
clear RelativeNCMutRate ExpectedNCMutRate exp_muts obs_muts obs_exp gnms keep

% drop the top 10 genes with the highest mutations / length in each tumortype
outliers = {};
keep = ~isnan(length_AAs{:,:});
mutrates = sum(mut{keep,ismember(ttypes,ttypes_used)}) ./ sum(length_AAs{keep,:});
exp_muts = length_AAs{keep,:}* mutrates;
obs_exp = array2table(mut{keep,ismember(ttypes,ttypes_used)} ./ exp_muts,'RowNames',mut.Properties.RowNames(keep));
for t = 1:size(obs_exp,2)
    temp = sortrows(obs_exp,t,'descend');
    outliers = union(outliers,temp.Properties.RowNames(1:10));
end
out_top10 = ismember(mut.Properties.RowNames,outliers);
clear mutrates keep exp_muts t obs_exp temp outliers

keep = ~isnan(length_AAs{:,:}) & ~all(is_nexpr{:,ismember(ttypes,ttypes_used)},2) & ~out_top10 & ~out_gr_exp;
    % only keep mutations from genes which we have all the other data for, and expressed in at least some tumors

warning('off','MATLAB:xlswrite:AddSheet');
tum_sum = sum(mut{keep,:},1);
if exist('NewFigureGeneration\genes_ttype.xlsx','file')>0
    delete('NewFigureGeneration\genes_ttype.xlsx','NewFigureGeneration\genes_PID_ttype.xlsx');
end
   
for t = 1:width(mut); %for each tumor type being looked at
    mut_out = mut(keep,t);
    mut_out.LengthAAs =length_AAs{keep,:};
    mut_out.Properties.RowNames = mut.Properties.RowNames(keep);
    mut_out.SumOthers = sum(mut{keep,:},2) - mut{keep,t};
    TotalOtherMuts = sum(tum_sum) - tum_sum(t);
    %calculate pval using Poisson; expected is (sum(all others) / sum(comparitor)
    mut_out.exp = sum(mut{keep,t}) .* (mut_out.SumOthers ./ TotalOtherMuts);
    exp = mut_out.exp;
    mut_out.obs_exp = log2(mut{keep,t} ./ exp);
    mut_out.Pval = cdf('Poisson',mut{keep,t},exp);
    mut_out.sigdep = mut_out.Pval < 0.02 & mut_out.obs_exp < log2(.5);
    mut_out.never = mut_out.exp < 10;
    mut_out.Q = mafdr(mut_out.Pval,'BHFDR',true);
    hits = mut_out.Properties.RowNames(mut_out.sigdep); 
    tmp_out = OverlapCalculatorUniverse3(hits); %for these genes, check overlap with PID
    
    %record output
    writetable(sortrows(mut_out,'Q'),'NewFigureGeneration\genes_ttype.xlsx','Sheet',ttypes{t},'WriteRowNames',1);
    writetable(tmp_out,'NewFigureGeneration\genes_PID_ttype.xlsx','Sheet',ttypes{t},'WriteRowNames',1);
    
    if any(strcmp(ttypes_analyzed,ttypes{t}))
        switch ttypes{t}
            case 'SKCM'
                nm1 = 'TableS6_genes_melanoma_increased_pur';
                nm2 = 'TableS7_genesets_enriched_melanoma_increased_pur';
                tumortype = 'Melanoma';
            case 'LUAD'
                nm1 = 'TableS4_genes_lungadeno_increased_pur';
                nm2 = 'TableS5_genesets_enriched_lungadeno_increased_pur';
                tumortype = 'Lung adenocarcinoma';
        end
        %if these files already exist delete them
        if exist(horzcat('NewFigureGeneration\',nm1,'.txt'),'file')>0
            delete(horzcat('NewFigureGeneration\',nm1,'.txt'),horzcat('NewFigureGeneration\',nm2,'.txt'));
        end
        writetable(mut_out,horzcat('NewFigureGeneration\',nm1,'.txt'),'delimiter','\t','WriteRowNames',1);
        writetable(tmp_out,horzcat('NewFigureGeneration\',nm2,'.txt'),'delimiter','\t','WriteRowNames',1);
        %Make a histogram of the percent depletion of mutations of all
        %genes and uniquely negatively selected genes
        %Figure 4A
        never  = mut_out.never;
        ranged = 2.^mut_out.obs_exp;
        ranged(ranged > 3) = 2^3;
        ranged(ranged < 2^-3) = 2^-3;
        [~,edges] = histcounts(log2(ranged(~never)));
        edges =2.^edges;
        fig = figure();
        set(fig,FigProps);
        fig.Position(3:4) = [2.25 1.5];
        hold on
        h=histogram(ranged(~never),edges);
        h.EdgeColor=[1 1 1];
        h.FaceAlpha = 1;
        h2=histogram(ranged(~never & mut_out.sigdep),edges);
        h2.EdgeColor=[1 1 1];
        h2.FaceAlpha = 1;
        ax=gca;
        set(ax,AxProps);
        ax.XScale='log';
        ax.XTick = 2.^[-3 -2 -1 0 1 2 3];
        Label = ax.XTick;
        Label(ax.XTick<1) = 1./ax.XTick(ax.XTick<1);
        ax.XTickLabel=Label;
        xlabel('Fold Depletion | Fold Enrichment','FontName','Arial','FontSize',10)
        ylabel('No. Genes','FontName','Arial','FontSize',10)
        title(tumortype,'FontName','Arial','FontSize',10,'FontWeight','Normal');
        L = legend({'All genes','UnderSelection'},'Location','NE');
        L.Box = 'off';
        L.FontName = 'Arial';
        L.FontSize = 10;
        print(horzcat('NewFigureGeneration\Fig4A_histogram_',tumortype,'_depl_muts.svg'),'-dsvg');
        clear ranged never edges tmp L ax fig h h2 Label 
    end
end
clear hits exp TotalOtherMuts t i tum_sum nm1 nm2 keep tmp_out out_gr_exp out_top10 X2_out ttypes_analyzed tumortype


% Make figure 4C (melanoma DNA repair pathways)
if exist('NewFigureGeneration\Fig4C.txt','file')>0
    delete('NewFigureGeneration\Fig4C.txt');
end
fig_genes = {'FANCA','FANCC','FANCG','FANCM','FANCB','FANCF','C17orf70','FANCL','ATR','ATM',...
    'MRE11A','RAD50','NBN','TOPBP1','CHEK1','CHEK2','SMC1A','SMC3','STAG2','RAD21','BRCA1',...
    'FANCD2','FANCI','BRCA2','BACH1','TP53','BARD1','MLH1','MSH6','MSH2','BLM','BRIP1','DDB2','XPC','ATF1',...
    'PALB2','RAD51','PCNA','RAD54L','PRKDC','POLM','ERCC5','POLE','RFC1','ERCC2','RFC4','APEX2'};
mut_out = readtable('NewFigureGeneration\TableS6_genes_melanoma_increased_pur.txt','delimiter','\t');
%which genes are under purifying selection?
is_pur = ismember(fig_genes,mut_out.Row(logical(mut_out.sigdep)));
present = ismember(fig_genes,mut_out.Row);
t = table(is_pur',present','RowNames',fig_genes,'VariableNames',{'is_pur','present'});
writetable(t,'NewFigureGeneration\Fig4C.txt','delimiter','\t','WriteRowNames',1);
clear t is_pur mut_out fig_genes present

% Make supplementary figure 3
if exist('NewFigureGeneration\FigS3.txt','file')>0
    delete('NewFigureGeneration\FigS3.txt');
end
fig_genes = {'EGFR','ERBB2','ERBB3','ERBB4','PLCG1','ITPR1','IL6ST', 'NRAS','PIK3R1'};
mut_out = readtable('NewFigureGeneration\TableS4_genes_lungadeno_increased_pur.txt','delimiter','\t');
%which genes are under purifying selection?
is_pur = ismember(fig_genes,mut_out.Row(logical(mut_out.sigdep)));
present = ismember(fig_genes,mut_out.Row);
t = table(is_pur',present','RowNames',fig_genes,'VariableNames',{'is_pur','present'});
writetable(t,'NewFigureGeneration\FigS3.txt','delimiter','\t','WriteRowNames',1);

clear is_pur present t mut_out fig_genes

%% Identifying gene sets under purifying selection in specific tumor types
ttypes_used = {'LUAD','SKCM','COADREAD','LIHC','GBMLGG','BRCA'};  %tumors to build expected?
ttypes_analyzed = {'SKCM','LUAD'};
ttypes_analyzed = {'LUAD'};

ttypes_used = {'COADREAD','LIHC','GBMLGG','BRCA','UCEC'};
ttypes_analyzed = {'COADREAD'};

% delete output if they already exist
if exist('NewFigureGeneration\TableS8_genesets_melanoma_increased_pur.txt','file')>0 ...
        || exist('NewFigureGeneration\TableS9_genesets_lungadeno_increased_pur.txt','file')>0
    delete('NewFigureGeneration\TableS8_genesets_melanoma_increased_pur.txt', 'NewFigureGeneration\TableS9_genesets_lungadeno_increased_pur.txt');
end

% drop genes with > 2.5x more mutations than expected across tumors
keep = ~isnan(length_AAs{:,:});
RelativeNCMutRate = NoncodingMutRate{keep,:} ./ nanmean(NoncodingMutRate{keep,:});
ExpectedNCMutRate = RelativeNCMutRate * (sum(mut{keep,ismember(ttypes,ttypes_used)}(:)) / sum(3*length_AAs{keep,:})); 
exp_muts = ExpectedNCMutRate .* (3 * length_AAs{keep,:});
obs_muts = sum(mut{keep,ismember(ttypes,ttypes_used)},2);
obs_exp = obs_muts./exp_muts;
out_gr_exp = obs_exp > 2.5; 
gnms = mut.Properties.RowNames(keep);
out_gr_exp = ismember(mut.Properties.RowNames,gnms(out_gr_exp));
clear RelativeNCMutRate ExpectedNCMutRate exp_muts obs_muts obs_exp gnms keep

% drop the top 10 genes with the highest mutations / length in each tumortype
outliers = {};
keep = ~isnan(length_AAs{:,:});
mutrates = sum(mut{keep,ismember(ttypes,ttypes_used)}) ./ sum(length_AAs{keep,:});
exp_muts = length_AAs{keep,:}* mutrates;
obs_exp = array2table(mut{keep,ismember(ttypes,ttypes_used)} ./ exp_muts,'RowNames',mut.Properties.RowNames(keep));
for t = 1:size(obs_exp,2)
    temp = sortrows(obs_exp,t,'descend');
    outliers = union(outliers,temp.Properties.RowNames(1:10));
end
out_top10 = ismember(mut.Properties.RowNames,outliers);
clear mutrates keep exp_muts t obs_exp temp outliers

keepG = ~isnan(length_AAs{:,:}) & ~out_top10 & ~out_gr_exp; % only keep mutations from genes which we have all the other data for and aren't outliers


L = sum(in_set_mat(:,keepG),2);
expr_t = double(in_set_mat(:,keepG)) * double(is_expr{keepG,ismember(ttypes,ttypes_used)}); % # of genes expressed in each set, in each tumor type
keep = L >=10 & L <= 400 ...
    & all(expr_t > .5 .* repmat(L,1,size(ttypes_used,2)),2);
in_set_mat1 = in_set_mat(keep,keepG);
tum_sum = sum(mut{keepG,ismember(ttypes,ttypes_used)},1);
for t = 1:length(ttypes_analyzed)
    i = strcmp(ttypes,ttypes_analyzed{t});
    observed_muts = mut{keepG,i}; % mutations in tumor type
    SumOthers = sum(mut{keepG,ismember(ttypes,ttypes_used)},2) - mut{keepG,i};
    TotalOtherMuts = sum(tum_sum) - tum_sum(i);
    % expected is (sum(all others) / sum(comparitor)
    expected_muts = sum(mut{keepG,i}) .* (SumOthers ./ TotalOtherMuts);    
    tmp = table();
    tmp.observed_muts = double(in_set_mat1) * observed_muts; %sum up observed mutations in these sets
    tmp.expected_muts = double(in_set_mat1) * expected_muts; %sum up observed mutations in these sets
    tmp.Diff = cdf('Poisson',tmp.observed_muts,tmp.expected_muts); %compute a significance
    tmp.obs_exp = tmp.observed_muts ./ tmp.expected_muts; %also obs / expected
    tmp.SetSizes = sum(in_set_mat1,2);
    tmp.SetNames = SetNames(keep);
    %Now find p values
    rep = 1*10^4; %replicates
    expr = double(is_expr{keepG & any(in_set_mat(keep,:),1)',ismember(ttypes,ttypes_used)}); %precompute expressed 
    observed_m = observed_muts(any(in_set_mat1,1)); %prefilter observed muts
    expected_m = expected_muts(any(in_set_mat1,1)); %prefilter expected muts
    %Pull random gene sets of different lengths
    setsizesLim = [min(sum(in_set_mat1,2)),max(sum(in_set_mat1,2))];
    gns = 1:nnz(keepG & any(in_set_mat(keep,:),1)'); %gene labels
    SetSizes = tmp.SetSizes;
    [tmp.emp_pval,null_diff] = GeneSetRandomization(rep,observed_m,expected_m,setsizesLim,gns,tmp.Diff,expr,SetSizes);
    tmp.Q = mafdr(tmp.emp_pval,'BHFDR',1);
    %for very unlikely things, fit null_diff to estimate pvalues
    ssizes = setsizesLim(1):setsizesLim(2);
    tmp.est_pval = ones(size(tmp.emp_pval));
    for i = 1:size(null_diff,2)
        nd = null_diff(:,i);
        p = 1/length(nd) : 1/length(nd) : 1 ;
        Q = quantile(nd,p);
        x = -log10(Q);
        y = -log10(p);
        mdl1 = fitlm(x,y,'linear','exclude',Q>0.01 | y>3.5);
        s = tmp.SetSizes == ssizes(i);
        est = predict(mdl1,-log10(tmp.Diff(s)));
        tmp.est_pval(s) = 10.^-est;
    end

    clear i nd p Q x y mdl1 s est

    tmp = sortrows(tmp,'Q');
    tmp2 = tmp(:,{'SetNames', 'SetSizes', 'observed_muts', 'expected_muts', 'obs_exp', 'emp_pval', 'est_pval', 'Q'});
    i = strcmp(ttypes,ttypes_analyzed{t});
    switch ttypes{i}
           case 'SKCM'
                nm1 = 'TableS8_genesets_melanoma_increased_pur';
                
           case 'LUAD'
                nm1 = 'TableS9_genesets_lungadeno_increased_pur';
           case 'COADREAD'
               nm1 = 'Tables___genesets_coadread_increased_pur';
                
    end
    
    writetable(tmp2,horzcat('C:\Users\rmathis\Dropbox (MIT)\mutations and evolution\NewFigureGeneration\',nm1,'.txt')...
        ,'delimiter','\t','WriteRowNames',0);
     
end
clear tmp tmp2 nm1 SetSizes gns setsizesLim expected_m observed_m expr rep t i observed_muts...
    expected_muts keepG keep expr_t L  ttypes_used ttypes_analyzed SumOthers in_set_mat1 ia null_diff

% Make figure 4B
% load gene set mutation data
tumor_input = readtable(horzcat(path_cur,'\NewFigureGeneration\TableS8_genesets_melanoma_increased_pur.txt'),'ReadRowNames',0,'Delimiter','\t');
genesets_sel = tumor_input.SetNames(tumor_input.Q<0.055 & tumor_input.obs_exp<0.6); % genesets we found to be under purifying selection
anno = readtable(horzcat(path_cur,'annotated_setnames.txt'),'delimiter','\t','ReadRowNames',0,'ReadVariableNames',0);%load in annotations
anno.Properties.VariableNames ={'setnames','anno'};
tmp = array2table(genesets_sel,'VariableNames',{'setnames'});
tmp.anno = cell(size(tmp.setnames));
for s = 1:height(tmp)
     tmp.anno(s) = anno{strcmp(anno.setnames,tmp.setnames{s}),'anno'}; %for each set append annotation
end
tmp.anno = categorical(tmp.anno);
anno = categories(tmp.anno); %list of annotations
% for each annotation, which gene sets are annotated with it? How many
% genes are in those sets? What is the percentage depletion of mutations?
% first calculate observed and expected genes
i = strcmp(ttypes,'SKCM');
keepG = ~isnan(length_AAs{:,:}) & ~out_top10 & ~out_gr_exp; % only keep mutations from genes which we have all the other data for and aren't outliers
observed_muts = mut{keepG,i}; % mutations in tumor type
SumOthers = sum(mut{keepG,ismember(ttypes,ttypes_used)},2) - mut{keepG,i};
tum_sum = sum(mut{keepG,ismember(ttypes,ttypes_used)},1);
TotalOtherMuts = sum(tum_sum) - tum_sum(i);
% expected is (sum(all others) / sum(comparitor)
expected_muts = sum(mut{keepG,i}) .* (SumOthers ./ TotalOtherMuts);    
in_set_mat1 = in_set_mat(:,keepG);
anno_t = array2table(zeros(length(anno),2),'RowNames',anno,'VariableNames',{'obs','exp'});
anno_t.genes = cell(size(anno_t.obs));
for a = 1:length(anno);
    sts = tmp{tmp.anno==anno{a},'setnames'}; %sets annotated with that
    i = any(in_set_mat1(ismember(SetNames,sts),:),1); %genes in sets
    anno_t{a,'obs'} = sum(observed_muts(i));
    anno_t{a,'exp'} = sum(expected_muts(i));
    anno_t{a,'genes'} = {i'};
end
anno_t.obs_exp = anno_t.obs ./ anno_t.exp;
anno_t.depl = 1 - anno_t.obs_exp;
anno_t.growth = ismember(anno_t.Properties.RowNames,{'SHH','WNT','NFkappaB','PI3K','MAPK','INFgamma'});
anno_t.survival = ismember(anno_t.Properties.RowNames,{'SHH','DNA repair','NFkappaB','Proteasome'});
anno_t.immune_supression = ismember(anno_t.Properties.RowNames,{'TNFalpha'});
anno_t = sortrows(anno_t,'depl','descend');
anno_t2 = anno_t(any(anno_t{:,{'growth','survival','immune_supression'}},2),:);

%Of the pathways in this figure 4B: how many genes, out of the total genes
%under selection (in gene sets under purifying selection) are reflected?
total_genes = nnz(any(in_set_mat1(ismember(SetNames,genesets_sel),:),1));
tmp = false(size(in_set_mat1,2),1);
for a = 1:height(anno_t2)
    i = anno_t2{a,'genes'};
    tmp = tmp | i{1};
end
is_anno = nnz(tmp);

clear genesets_sel tumor_input anno tmp out_top10 out_gr_exp i observed_muts...
    expected_muts keepG in_set_mat1 SumOthers tum_sum a anno_t sts is_anno total_genes...
    tmp

%% Determining the fraction of essential genes under purifying selection 

% get the genes in sets
setgenes = unique(vertcat(MolSigDB{:}));
genenms = mut.Properties.RowNames(~isnan(length_AAs{:,:}) & ~isnan(NoncodingMutRate{:,:}));
setgenes = setgenes(ismember(setgenes,genenms)); %filter for genes analyzed

% grab the genes under general purifying selection across tumors
alltumor_input = readtable(horzcat(path_cur,'\NewFigureGeneration\TableS3_genesets_alltumors.txt'),'ReadRowNames',0,'Delimiter','\t');
genesets_sel = alltumor_input.SetNames(alltumor_input.Q<0.05 & alltumor_input.obs_exp<0.8); % genesets we found to be under purifying selection
temp2 = MolSigDB(ismember(SetNames,genesets_sel));
genes_sel = unique(vertcat(temp2{:}));
pur_data = array2table(ismember(setgenes,genes_sel),'RowNames',setgenes,'VariableNames',{'Gp'});
clear genes_sel genesets_sel alltumor_input

% grab the genes under increased purifying selection in melanomas
tumor_input = readtable(horzcat(path_cur,'\NewFigureGeneration\TableS8_genesets_melanoma_increased_pur.txt'),'ReadRowNames',0,'Delimiter','\t');
genesets_sel = tumor_input.SetNames(tumor_input.Q<0.055 & tumor_input.obs_exp<0.6); % genesets we found to be under purifying selection
temp2 = MolSigDB(ismember(SetNames,genesets_sel));
genes_sel = unique(vertcat(temp2{:}));
pur_data.GpM = ismember(setgenes,genes_sel);
clear temp2 genenms genes_sel setgenes genesets_sel tumor_input
tmp = array2table(zeros(1,2),'VariableNames',{'overlap','n_essential'});
tmp.essential_datasource = cell(1,1);
tmp.purifying_datasource = cell(1,1);
n=1;
for i = 1:2
    for s =1:3
        switch s
            case 1 %use Hart et al data                    
                ess_cls_t = hart_et_al(:,'numTKOHits');  %for each gene, number of cell lines in which it is deemed essential
                L1 = 'Hart et al';
                cln = 5;
            case 2 %use Tzelepis et al data         
                ess_cls_t = tz_data(:,'x_All');  %for each gene, number of cell lines in which it is deemed essential
                L1 = 'Tzelepis et al';
                cln = 7;
            case 3 % use Wang et al data
                ess_cls_t = w_data(:,'x_All');
                L1 = 'Wang et al';
                cln = 4;
        end
        genes = intersect(pur_data.Properties.RowNames,ess_cls_t.Properties.RowNames);
        ess_cls = ess_cls_t{genes,:};
        switch i
            case 1
                L2 = 'All Tumor Types';
                Pur = pur_data{genes,'Gp'};
            case 2 
                L2 = 'Melanoma';
                Pur = pur_data{genes,'GpM'};
        end
        Ess = ess_cls >= cln - 2;
        x = nnz(Ess & Pur); % overlap
        %M = length(Pur); % num genes in universe
        K = nnz(Ess); % num essential genes
        %N = nnz(Pur); % num genes under purifying selection
        %p = hygecdf(x,M,K,N,'upper');
        tmp{n,{'overlap','n_essential'}} = [x,K];
        tmp{n,'essential_datasource'} = {L1};
        tmp{n,'purifying_datasource'} = {L2};
        n=n+1;
        clear x M K N Ess Pur
    end
end

% Make figure 4D
t = horzcat( tmp.n_essential - tmp.overlap,tmp.overlap);
fig = figure();
set(fig,FigProps)
fig.Position(3:4) = [3 2.5];
pieprops.FontName = 'Arial';
pieprops.FontSize = 8;
pieprops.Color = [1 1 1];
for i = 1:size(t,1);
    subplot((size(t,1))/3,3,i)
    P1 = pie(t(i,:));
    ax = gca;
    set(ax,AxProps);
    Ps = findobj(P1,'Type','Text');
    set(Ps,pieprops);
    P1(1).FaceColor = [0    0.4470    0.7410];
    P1(3).FaceColor = [0.8500    0.3250    0.0980];
    P1(2).Position(1:2) = [-0.56 .5];
    P1(4).Position(1:2) =[0.56 .5];
    if i <=3;
        title(tmp.essential_datasource{i},'FontName','Arial','FontSize',10,'FontWeight','Normal');
    end
end
print('NewFigureGeneration\Fig4D_pie.svg','-dsvg');
clear L cln ess_cls_t s i pur_data L1 L2 ess_cls genes tmp n t fig ax Ps pieprops P1 


%what happens if we reduce the number of mutations in the data?  eg, if you
%droped mutations from the data, how many fewer hits would be returned?
%for each tumor type

%%  There exists an estimate of the amount of problematicness for observed missense mutations.
% If purifying selection dropped important mutations, we should see the
% most damaging mutations gone, right?  Non-expressed genes as a control.

%% Is there any evidence of increased purifying genes in chromosomes that are dropped to one copy?  
% Split chromosomes into dropped more frequently, not dropped more frequently.  
% Determine purifying selection in genes on these two groups of chromosomes and
% compare the distributions.  

%genes under purifying selection?
alltumor_input = readtable(horzcat(path_cur,'\NewFigureGeneration\TableS3_genesets_alltumors.txt'),'ReadRowNames',0,'Delimiter','\t');
genesets_sel = alltumor_input.SetNames(alltumor_input.Q<0.05 & alltumor_input.obs_exp<0.8); % genesets we found to be under purifying selection
temp2 = MolSigDB(ismember(SetNames,genesets_sel));
genes_sel = unique(vertcat(temp2{:}));
genenms = mut.Properties.RowNames(~isnan(length_AAs{:,:}) & ~isnan(NoncodingMutRate{:,:}));
genes_sel = genes_sel(ismember(genes_sel,genenms)); %filter for genes analyzed
is_pur = ismember(mut.Properties.RowNames,genes_sel);
clear temp2 genenms clear setgenes genenms genes_sel 
  
Chr_tg = table();
Chr_tg(:,{'Pos','Chr'}) = GetChrData(path_cur,mut.Properties.RowNames); %chromosome and position


ttypes_analyzed = {'SKCM','LUAD'};
chr = horzcat({'X'},arrayfun(@num2str,1:22,'UniformOutput',false));
for t = 1:length(ttypes_analyzed)
    %observed and expected mutations
    keep = ~isnan(length_AAs{:,:}) & ~isnan(NoncodingMutRate{:,:}) & is_expr{:,ismember(ttypes,ttypes_analyzed{t})} & ~cellfun(@isempty,Chr_tg.Chr); % only keep mutations from genes which we have all the other data for, and are expressed;
    gene_t = table();
    gene_t.genenms = mut.Properties.RowNames(keep); % gene names
    gene_t.observed_muts = mut{keep,ismember(ttypes,ttypes_analyzed{t})}; % sum mutations in these tumor types
    RelativeNCMutRate = NoncodingMutRate{keep,:} ./ mean(NoncodingMutRate{keep,:});
    ExpectedNCMutRate = RelativeNCMutRate * sum(gene_t.observed_muts) / sum(3*length_AAs{keep,:}); 
    gene_t.expected_muts = ExpectedNCMutRate .* (3 * length_AAs{keep,:}); % calculate expected muts
    gene_t.Properties.RowNames = gene_t.genenms;
    gene_t.is_expr =  is_expr{keep,ismember(ttypes,ttypes_analyzed{t})};
    gene_t.Chr = Chr_tg{keep,'Chr'};
    gene_t.is_pur = is_pur(keep);
    
    % Which chromosomes are depleted or amplified in this tumor type?
    chr_i = readtable(horzcat(path_cur,'GISTIC\',ttypes_analyzed{t},'_GISTIC2_arm_level.txt'),'delimiter','\t');
    chr_i.chr = cellfun(@(x) x(1:end-1),chr_i.Arm,'UniformOutput',false);
    
    % chromosomes depleted, amplified?
    % for each chromosome, determine the observed and expected mutations
    chr_t = array2table(false(length(chr),2),'RowNames',chr,'VariableNames',{'is_dep','is_amp'});
    chr_t.observed = zeros(size(chr_t.is_dep));
    chr_t.expected = zeros(size(chr_t.is_dep));
    
    for c = 1:height(chr_t)
        chr_t{c,'is_dep'} = all(chr_i{strcmpi(chr_i.chr,chr{c}),'DelQ_value'} < 0.05);
        chr_t{c,'is_amp'} = all(chr_i{strcmpi(chr_i.chr,chr{c}),'AmpQ_value'} < 0.05);
        chr_t{c,'observed'} = sum(gene_t{strcmpi(gene_t.Chr,chr{c}) & gene_t.is_expr & gene_t.is_pur,'observed_muts'});
        chr_t{c,'expected'} = sum(gene_t{strcmpi(gene_t.Chr,chr{c}) & gene_t.is_expr & gene_t.is_pur,'expected_muts'});
    end
    chr_t.perc_depl = 1-(chr_t.observed./chr_t.expected);
    chr_t.L = cell(size(chr_t.is_dep));
    chr_t.L(chr_t.is_dep) = {'Dep'};
    chr_t.L(chr_t.is_amp) = {'Amp'};
    chr_t.L(~chr_t.is_dep & ~chr_t.is_amp) = {'Normal'};
    

    fig = figure();
    set(fig,FigProps);
    fig.Position(3:4) = [7 1.5];
    X = 1;
    hold on
    label = cell(size(chr_t.is_dep));
    lgnd = cell(3,1);
    for s = 1:3
        switch s
            case 1
                L = 'Depleted';
                i = chr_t.is_dep;
                clr = 'r';
            case 2
                L = 'Normal';
                i = ~chr_t.is_dep & ~chr_t.is_amp;
                clr = 'k';
            case 3
                L = 'Amplified';
                i = chr_t.is_amp;
                clr = 'b';
                
        end
        
        bar(X:X - 1 + nnz(i),chr_t{i,'perc_depl'},clr); %make the bar graph   
        label(X:X - 1 + nnz(i)) = chr_t.Properties.RowNames(i);
        X = X + nnz(i); %iterate up x
        lgnd{s} = L;
    end
    legend(lgnd,'Box','off','FontName','Arial','FontSize',10,'Location','eastoutside')
    ax = gca;
    set(ax,AxProps)
    ax.XTick = 1:length(label);
    ax.XTickLabels = label;
    fig.Color = [1 1 1];
    ax.Color = [1 1 1];
    
end

    
clear chr_t c chr_i t chr ttypes_analyzed chr_pos alltumor_data setgenes genes_sel genesets_sel...
    s gene_t gene_t2 gene_t3 b ax fig L i clr lgnd label X keep is_pur ExpectedNCMutRate RelativeNCMutRate


%% Of the genes we found in to be under increased purifying selection, what fraction of them are druggable?

% get the genes in sets
setgenes = unique(vertcat(MolSigDB{:}));
genenms = mut.Properties.RowNames(~isnan(length_AAs{:,:}) & ~isnan(NoncodingMutRate{:,:}));
setgenes = setgenes(ismember(setgenes,genenms)); %filter for genes analyzed

% grab the genes under increased purifying selection in melanomas
tumor_input = readtable(horzcat(path_cur,'\NewFigureGeneration\TableS8_genesets_melanoma_increased_pur.txt'),'ReadRowNames',0,'Delimiter','\t');
genesets_sel = tumor_input.SetNames(tumor_input.Q<0.055 & tumor_input.obs_exp<0.6); % genesets we found to be under purifying selection
temp2 = MolSigDB(ismember(SetNames,genesets_sel));
genes_sel = unique(vertcat(temp2{:}));
genes_sel = genes_sel(ismember(genes_sel,setgenes));

%list of druggable genes?
drg = readtable(horzcat(path_cur,'aag1166_Table S1.xlsx'),'Sheet','Data');
drg = drg(ismember(drg.hgnc_names,genes_sel),:);
drg.Properties.RowNames = drg.hgnc_names;
%Results as a logical-table
r = array2table(false(length(genes_sel),2),'VariableNames',{'small_mol_druggable','bio_druggable'});
r.small_mol_druggable(ismember(genes_sel,drg.hgnc_names)) = strcmpi(drg{genes_sel(ismember(genes_sel,drg.hgnc_names)),'small_mol_druggable'},'Y');
r.bio_druggable(ismember(genes_sel,drg.hgnc_names)) = strcmpi(drg{genes_sel(ismember(genes_sel,drg.hgnc_names)),'bio_druggable'},'Y');
ng = size(r,1);
nsmd = nnz(r.small_mol_druggable);
nbd = nnz(r.bio_druggable);
nd = nnz(r.small_mol_druggable | r.bio_druggable);

clear genes_sel temp2 genesets_sel tumor_input setgenes genenms setgenes drg r
    

