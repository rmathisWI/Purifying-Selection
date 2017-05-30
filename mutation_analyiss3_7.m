do_mut_context = 1;
expr_threshold = 2^3;
stringency = 0.95;

quant_results = table();
[mutfilename,mutpathname] = uigetfile('*.txt','open mutation file');
[expfilename,exppathname] = uigetfile('*.txt','open expression file');

%% tumor type?
switch mutfilename(1:4)
    case 'LUAD'
        tumortype = 'Lung adenocarcinoma';
    case 'SKCM'
        tumortype = 'Melanoma';
    case 'LUSC'
        tumortype = 'Lung squamous cell carcinoma';
    case 'GBML'
        tumortype = 'Glioma';
    case 'COAD'
        tumortype = 'Colorectal adenocarcinoma';
end
        

%% load muation data
try 
    mut=datastore(horzcat(mutpathname,mutfilename),'delimiter','\t','whitespace',' \b','SelectedVariableNames',{'Hugo_Symbol','Variant_Type','Codon_Change','Tumor_Sample_Barcode','dbSNP_RS','is_coding','is_missense','is_silent','is_nonsense','Reference_Allele','Tumor_Seq_Allele2','Transcript_Strand','ref_context','Chromosome','Start_Position','Variant_Classification','Protein_Change'});
catch
        mut=datastore(horzcat(mutpathname,mutfilename),'delimiter','\t','whitespace',' \b','SelectedVariableNames',{'Hugo_Symbol','Variant_Type','Tumor_Sample_Barcode','dbSNP_RS','is_coding','is_missense','is_silent','is_nonsense','Variant_Classification'});
end
mut.TextscanFormats(1,14) ={'%s'};
mut=readall(mut);
mut.not_SNP = strcmp(table2cell(mut(:,'dbSNP_RS')),'') | strcmp(table2cell(mut(:,'dbSNP_RS')),'novel');

%% load expression data 
expr=readall(datastore(horzcat(exppathname,expfilename),'delimiter','\t','whitespace',' \b','NumHeaderLines',1));
expr_genelist = table2cell(expr(:,1));
for x=1:length(expr_genelist)
    temp = expr_genelist(x,1);
    temp = strsplit(temp{1},'|');
    expr_genelist(x) = temp(1,1);
end
clear temp x
expr(:,1) = expr_genelist;
exprdata = expr(:,2:end);

%expr data thresholded; look up mutant genes for expression
exprdata_by_tumor=exprdata{:,:}>expr_threshold;
expr_by_tumor = sum(exprdata_by_tumor,2);
expr_by_tumor_geomean = geomean(1+exprdata{:,:},2);

%expressed in >95% of tumors?  Expressed in <5% of tumors?
is_expressed = expr_by_tumor> stringency*size(exprdata,2);
is_notexpressed = expr_by_tumor< (1-stringency)*size(exprdata,2);
exprdata_by_mut = horzcat(cell2table(expr{:,1}),array2table(is_expressed),array2table(is_notexpressed),array2table(expr_by_tumor_geomean));
exprdata_by_mut = exprdata_by_mut(~strcmp(exprdata_by_mut{:,1},'?'),:);
[~,ia,~] = unique(exprdata_by_mut{:,1});
exprdata_by_mut = exprdata_by_mut(ia,:);
exprdata_by_mut.Properties.RowNames = exprdata_by_mut{:,1};
exprdata_by_mut = exprdata_by_mut(:,2:end);
clear  exprdata_by_tumor expr expr_by_tumor expr_by_tumor_geomean expr_genelist exprdata is_expressed is_notexpressed ia

%% tabulate missense mutations
%strcmp(mut.Variant_Type,'SNP') & ~strcmp(mut.Protein_Change,'') & ~strcmp(mut.Protein_Change,'.')& cellfun(@isempty,strfind(mut.Protein_Change,'_splice')) & cellfun(@isempty,strfind(mut.Protein_Change,'ins'))
mis1 = tabulate(mut.Hugo_Symbol(mut.is_missense & mut.is_coding & mut.not_SNP & strcmp(mut.Variant_Type,'SNP')));
mis1 = tabulate(mut.Hugo_Symbol(mut.is_missense & mut.is_coding & strcmp(mut.Variant_Type,'SNP')));
mis_by_gene = table(mis1(:,1),cell2mat(mis1(:,2)),'VariableNames',{'genenames','mis_by_gene'});
clear mis1


%% load gene length data
genelengths =readtable('uniprot_length_16_06_27.xlsx','ReadRowNames',true,'ReadVariableNames',1);
%genelengths =readtable('uniprot-len.xlsx','ReadRowNames',true,'ReadVariableNames',1);
genelengths.Properties.VariableNames = {'length_AAs'};

%% for each gene: look up missense mutation rate, length, and expression;
%first, only want genes we have length and expression data for
genes = intersect(genelengths.Properties.RowNames,exprdata_by_mut.Properties.RowNames);
genelengths = genelengths(ismember(genelengths.Properties.RowNames,genes),:);
exprdata_by_mut = exprdata_by_mut(ismember(exprdata_by_mut.Properties.RowNames,genes),:);
result_t = join(genelengths,exprdata_by_mut,'keys','RowNames');
result_t = horzcat(result_t.Properties.RowNames,result_t);

%now look up mutations
mis_by_gene = mis_by_gene(ismember(mis_by_gene.genenames,genes),:);
result_t = outerjoin(result_t,mis_by_gene,'LeftKey',result_t.Properties.VariableNames(1),'RightKey',{'genenames'},'MergeKeys',1);
result_t.mis_by_gene(isnan(result_t.mis_by_gene)) = 0;


%how many tumors were sequenced?
tumor_number = length(unique(cellfun(@(x) {x(1:12)},mut.Tumor_Sample_Barcode,'UniformOutput',true)));

%for each gene calculate the log10 number of mutations/aa/scaled#
%tumors
scale_factor = 575 * 1000; %average gene size * 1000 tumors
result_t.norm_mis_by_gene = log10((1+result_t.mis_by_gene)./result_t.length_AAs /(tumor_number / scale_factor)); 

%% What's the codon count of each gene?
cdncounts =readtable('codon_counttable_1.txt','delimiter','\t');
cdncounts = cdncounts(ismember(cdncounts.Row,genes),:);
cdncount = table(table2array(cdncounts(:,2:end)),'VariableNames',{'cdncount_by_gene'});
cdncount.Row = cdncounts.Row;
result_table = outerjoin(result_t,cdncount,'LeftKey',result_t.Properties.VariableNames(1),'RightKey',{'Row'},'MergeKeys',1);
result_table.Properties.VariableNames{1} = 'mut_genes';
clear cdncount cdncounts consv exprdata_by_mut genelengths genes mis_by_gene result_t syn_by_gene expr


FontSize=12;
FontName='Arial';
linewidth = 1.5;


%% Let's look at genes that, when KOed, drop growth.  Are they less mutated in tumors?
%download data from a CRISPR screen for essentiality
%DOI: 10.1126/science.aac7041

growth_phenotype_input = readtable('aac7041_SM_Table_S3.txt','delimiter','\t','ReadRowNames',true);
growth_phenotype = growth_phenotype_input(:,:);

%take the mean essentiality scores from the different cell lines
growth_phenotype.essentiality_score = mean([growth_phenotype.KBM7CS,growth_phenotype.K562CS,growth_phenotype.JiyoyeCS,growth_phenotype.RajiCS],2);

%use thresholding to define essential, not essential
essentiality_pos_thresh = -1.5;
essentiality_neg_thresh = 0;
growth_phenotype.log_essential = growth_phenotype.essentiality_score < essentiality_pos_thresh;
growth_phenotype.log_not_essential = growth_phenotype.essentiality_score > essentiality_neg_thresh;

%add this to the result table
result_table.is_essential = ismember(result_table.mut_genes,growth_phenotype.Properties.RowNames(growth_phenotype.log_essential));
result_table.is_notessential = ismember(result_table.mut_genes,growth_phenotype.Properties.RowNames(growth_phenotype.log_not_essential));


%merge growth phenotype table with a copy of result_table and sort 
min_mut=0;
result_table_temp = result_table;
result_table_temp.Properties.RowNames = result_table.mut_genes;
result_table_temp = result_table_temp(ismember(result_table_temp.Properties.RowNames,growth_phenotype.Properties.RowNames) & result_table_temp.is_expressed & result_table_temp.mis_by_gene >= min_mut,{'norm_mis_by_gene','expr_by_tumor_geomean','mis_by_gene','length_AAs'});
essentiality_t = join(result_table_temp,growth_phenotype,'Key','RowNames');
essentiality_t = essentiality_t(~isinf(essentiality_t.norm_mis_by_gene) & ~isnan(essentiality_t.norm_mis_by_gene),:);
essentiality_ts = sortrows(essentiality_t,'essentiality_score');
ess = essentiality_ts.essentiality_score;
mis = essentiality_ts.mis_by_gene;
len = essentiality_ts.length_AAs;


%split genes into essential and not essential
min_mut = 0;

%take from the table normalized missense rates for essential, ~essential
yesess=essentiality_t(essentiality_t.log_essential & essentiality_t.mis_by_gene > min_mut,:);
notess=essentiality_t(essentiality_t.log_not_essential & essentiality_t.mis_by_gene > min_mut,:);
p_es = ranksum(yesess.norm_mis_by_gene,notess.norm_mis_by_gene,'tail','both');
quant_results.P_essential = p_es;
quant_results.FC_essential = (sum(yesess.mis_by_gene)/sum(yesess.length_AAs)) / (sum(notess.mis_by_gene) / sum(notess.length_AAs));
%quant_results.FC_essential = 1/(10^(median(yesess)-median(notess)));
if log10(p_es) < -20
    p_es = 1*10^-20;
end

%bin genes into bins of 100 genes by essentiality score and plot against the median of mutation rate
w=100;
%how many windows?
n=floor(height(essentiality_ts)/w);

mean_essentiality = zeros(n,1);
mean_mutrate = zeros(n,1);
for i=1:n;
    mean_essentiality(i)=mean(ess((i-1)*w+1:i*w));
    mean_mutrate(i) = sum(mis((i-1)*w+1:i*w))/sum(len((i-1)*w+1:i*w))/(tumor_number/scale_factor);
end
figure()
fig = gcf;
fig.InvertHardcopy='off'; 
fig.Units='inches'; 
fig.Position(3:4) = [4.5 3]; 
fig.PaperPositionMode='Auto'; 
fig.Renderer='painters';
scatter(mean_essentiality,mean_mutrate,25,'filled');
xlabel({'Gene essentiality score','(CRISPR screen)'})
ylabel({'Gene mutation rate','(patient tumors)'})
title(tumortype);
t=horzcat('p < 10^{',num2str(ceil(log10(p_es))),'}');
text(.85,0.1,t,'FontSize',FontSize,'FontName',FontName,'units','normalized')
ax=gca; 
ax.FontName=FontName; 
ax.FontSize=FontSize; 
ax.LabelFontSizeMultiplier=1; 
ax.TitleFontSizeMultiplier=1; 
ax.TickDir='out'; 
ax.LineWidth=linewidth; 
ax.Box='off'; 
ax.Color='none';
fig.Color='none';
print(horzcat(mutfilename(1:4),'essentiality.svg'),'-dsvg','-painters');
fig.Color=[1,1,1];

% make a pie chart of mutations lost:
% What fraction of mutations are estimated to be lost from selection (essentiality)
figure()
fig = gcf; 
fig.InvertHardcopy='off';
fig.Units='inches';
fig.Position(3:4) = [1 1]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
left = quant_results.FC_essential;
P = pie([left,1-left],{'',horzcat(num2str(round((1-left)*100)),'%')});
P(4).FontSize = FontSize;
P(4).FontName = FontName;
P(3).LineStyle = '--';
P(3).EdgeColor = [.4 .4 .4];
colormap([0.8,0.8,0.8; 1,1,1])
ax=gca; 
ax.Color='none';
fig.Color='none';
print(horzcat(mutfilename(1:4),'_pie_essential.svg'),'-dsvg');
fig.Color = [1 1 1]; 


% Does this still work if we control for expression?
yesess=essentiality_t(essentiality_t.log_essential & essentiality_t.mis_by_gene > min_mut,:);
notess=essentiality_t(essentiality_t.log_not_essential & essentiality_t.mis_by_gene > min_mut,:);
r=1000;
FC = zeros(r,1);
yesess.logexpr = log2(yesess.expr_by_tumor_geomean);
notess.logexpr = log2(notess.expr_by_tumor_geomean);
%pull from notess to match the distribution of expression levels in yesess
yesess_mr = yesess;
notess_sampr =  cell(height(yesess),r); 
for g = 1:size(notess_sampr,1) %for each essential gene
    log_expr = notess.Properties.RowNames(abs(notess.logexpr - yesess.logexpr(g)) < 0.005); %find not essential genes with the same expression level
    if ~isempty(log_expr)
        notess_sampr(g,:) = datasample(log_expr,r,'Replace',true);
        %for i = 1:r
        %    notess_sampr(g,i) = randsample(log_expr,1); %pull a random one and get its mutation rate
        %end
    else
        notess_sampr(g,:) = {''};          
    end
end
%drop empties
logempty = cellfun(@isempty,notess_sampr(:,1));
yesess_mr(logempty,:) = [];
notess_sampr(logempty,:) = [];
essmr = (sum(yesess_mr.mis_by_gene) / sum(yesess_mr.length_AAs));
for i = 1:r
    %FC(i) = (1/10^(nanmedian(yesess_mr) - nanmedian(notess_sampr(:,i))));
    FC(i)= essmr / (sum(notess{notess_sampr(:,i),'mis_by_gene'}) / sum(notess{notess_sampr(:,i),'length_AAs'}));   
end
figure()
fig = gcf;
fig.InvertHardcopy='off';
fig.Units='inches';
fig.Position(3:4) = [4.5 3]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters'; 
h=histogram(FC,40,'Normalization','Probability');
h.EdgeColor=[1 1 1];
h.FaceAlpha = 1;
hold on
scatter(quant_results.FC_essential,0,100,'filled');
ax=gca; 
ax.FontName=FontName; 
ax.FontSize=FontSize;
ax.LabelFontSizeMultiplier=1;
ax.TitleFontSizeMultiplier=1;
ax.TickDir='out';
ax.LineWidth=linewidth;
ax.Box='off';
ax.Color='none';
xlabel({'Fold decrease in essential','gene mutation rate'})
ylabel('Fraction of permutations')
ax.TitleFontSizeMultiplier=1;
title(tumortype);
L=legend({'Expression-normalized','Observed'},'Location','northwest','FontSize',FontSize,'FontName',FontName);
L.Box = 'off';
L.FontSize = FontSize;
L.FontName = FontName;
fig.Color='none';
print(horzcat(mutfilename(1:4),'histogram_expressionnorm_essentiality.svg'),'-dsvg','-painters');
fig.Color=[1,1,1];

%% okay, let's try testing this correlation after normalizing for expression differences
%strategy 1: bin by expression rates, calculate mutation rate and
%essentiality.  Correlate essentiality and mutation rates.
essentiality_t.logexpr =  log2(essentiality_t.expr_by_tumor_geomean); % add expression score;
%bin by expression: 
bs = 0.025; %bin size;
bs = 0.25;
mean_mutrate = zeros(height(essentiality_t),1);
mean_expression = zeros(height(essentiality_t),1);
flag = zeros(size(mean_mutrate));
binsize = zeros(size(mean_mutrate));
for b = 1:height(essentiality_t); %for each gene
    ge = essentiality_t.logexpr(b);
    gb = essentiality_t(essentiality_t.logexpr >= (ge - bs) & essentiality_t.logexpr <= (ge + bs),:); %pull genes close in expr
    mean_mutrate(b) = sum(gb.mis_by_gene) / sum(gb.length_AAs); %compute mutrate
    mean_expression(b) = mean(gb.logexpr); %compute average expression
    binsize(b) = height(gb);
    if height(gb) < 50;
        flag(b) = -1;
    end
end
essentiality_t.expected = mean_mutrate .* essentiality_t.length_AAs;
essentiality_t.obs_exp = essentiality_t.mis_by_gene ./ essentiality_t.expected;
p_es = ranksum(essentiality_t.obs_exp(essentiality_t.log_essential),essentiality_t.obs_exp(essentiality_t.log_not_essential));
quant_results.P_essential_nrm = p_es;
quant_results.FC_essential_nrm = (sum(essentiality_t.mis_by_gene(essentiality_t.log_essential)) / sum(essentiality_t.expected(essentiality_t.log_essential))) / ...
(sum(essentiality_t.mis_by_gene(essentiality_t.log_not_essential)) / sum(essentiality_t.expected(essentiality_t.log_not_essential)));
essentiality_ts = sortrows(essentiality_t,'essentiality_score');
ess = essentiality_ts.essentiality_score;
mis = essentiality_ts.mis_by_gene;
exp = essentiality_ts.expected;
w=100;
%how many windows?
n=floor(height(essentiality_ts)/w);
mean_essentiality = zeros(n,1);
mean_mutrate = zeros(n,1);
for i=1:n;
    mean_essentiality(i)=mean(ess((i-1)*w+1:i*w));
    mean_mutrate(i) = sum(mis((i-1)*w+1:i*w))/sum(exp((i-1)*w+1:i*w));
end
figure()
fig = gcf;
fig.InvertHardcopy='off'; 
fig.Units='inches'; 
fig.Position(3:4) = [4.5 3]; 
fig.PaperPositionMode='Auto'; 
fig.Renderer='painters';
scatter(mean_essentiality,mean_mutrate,25,'filled');
xlabel({'Gene essentiality score','(CRISPR screen)'})
ylabel({'Gene mutation rate','(patient tumors)'})
title(tumortype);
t=horzcat('p < 10^{',num2str(ceil(log10(p_es))),'}');
text(.85,0.1,t,'FontSize',FontSize,'FontName',FontName,'units','normalized')
ax=gca; 
ax.FontName=FontName; 
ax.FontSize=FontSize; 
ax.LabelFontSizeMultiplier=1; 
ax.TitleFontSizeMultiplier=1; 
ax.TickDir='out'; 
ax.LineWidth=linewidth; 
ax.Box='off'; 
ax.Color='none';
fig.Color='none';
print(horzcat(mutfilename(1:4),'essentiality_exprnorm.svg'),'-dsvg','-painters');
fig.Color=[1,1,1];

% What fraction of mutations are estimated to be lost from selection (expression normalized essentiality)
figure()
fig = gcf; 
fig.InvertHardcopy='off';
fig.Units='inches';
fig.Position(3:4) = [1 1]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
left = quant_results.FC_essential_nrm;
P = pie([left,1-left],{'',horzcat(num2str(round((1-left)*100)),'%')});
P(4).FontSize = FontSize;
P(4).FontName = FontName;
P(3).LineStyle = '--';
P(3).EdgeColor = [.4 .4 .4];
colormap([0.8,0.8,0.8; 1,1,1])
ax=gca; 
ax.Color='none';
fig.Color='none';
print(horzcat(mutfilename(1:4),'_pie_exprnorm_essential.svg'),'-dsvg');
fig.Color = [1 1 1]; 


clear ax s fig i n w mean_essentiality mean_mutrate result_table_temp mis growth_phenotype_input ess essentiality_ts p_es len

%% What fraction of mutations are estimated to be lost from selection?
%First grab all of the missense mutations from expr and nonexpr
temp_t = result_table(~isnan(result_table.norm_mis_by_gene) & ~isinf(result_table.norm_mis_by_gene) & ~(~result_table.is_expressed & ~result_table.is_notexpressed),:);
temp_expr_mis = temp_t{temp_t.is_expressed,'norm_mis_by_gene'};
temp_noexpr_mis = temp_t{temp_t.is_notexpressed,'norm_mis_by_gene'};
%make ecdf of expr mut rates 
[q1,m1]=ecdf(temp_expr_mis);
%find nonexpr mut rates of equivilent quantiles
m2=quantile(temp_noexpr_mis,q1);
%convert mut rates of expr to mut rates of equivilent quantiles in nonexpr
expected_t = table(m1,m2);
expected_expr_normmis = zeros(size(temp_expr_mis));
for x=1:length(expected_expr_normmis)
    expected_expr_normmis(x) = min(expected_t{expected_t.m1 == temp_expr_mis(x),'m2'});
end
%turn back into # of missense mutations
expected_expr_mis = (((10.^expected_expr_normmis) * (tumor_number / scale_factor)) .* temp_t{temp_t.is_expressed,'length_AAs'}) -1;
%sum up observed, expected mutationsc
obs_exp = [sum(temp_t.mis_by_gene) sum([temp_t{temp_t.is_notexpressed,'mis_by_gene'}; expected_expr_mis])];
%add to both the observed mutations in unclear expression genes 
obs_exp = obs_exp + sum(result_table{~isnan(result_table.norm_mis_by_gene) & ~isinf(result_table.norm_mis_by_gene) & (~result_table.is_expressed & ~result_table.is_notexpressed),'mis_by_gene'});
lost = (obs_exp(2)-obs_exp(1))/obs_exp(2);
figure()
fig = gcf; 
fig.InvertHardcopy='off';
fig.Units='inches';
fig.Position(3:4) = [1 1]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
quant_results.Avg_muts_lost_to_selection = round((obs_exp(2)-obs_exp(1))/tumor_number);
quant_results.percent_muts_lost_to_selection = (obs_exp(2)/obs_exp(1))*100;
P = pie([1-lost,lost],{num2str(round(obs_exp(1)/tumor_number)),''});
P(2).FontSize = FontSize;
P(2).FontName = FontName;
P(3).LineStyle = '--';
colormap([0.8,0.8,0.8; 1,1,1])
ax=gca; 
ax.Color='none';
fig.Color='none';
print(horzcat(mutfilename(1:4),'_pie_expr.svg'),'-dsvg');
fig.Color = [1 1 1]; 
clear ax fig lost obs_exp expected_expr_mis expected_expr_normmis expected_t m2 m1 q1 temp_noexpr-_mis temp_expr_mis temp_t



%% plot CDFs of mutations in expressed and not-expressed genes
min_mut = 1;
expr_mis_t = result_table(result_table.is_expressed & result_table.mis_by_gene >= min_mut & ~isnan(result_table.norm_mis_by_gene) & ~isinf(result_table.norm_mis_by_gene),:);
temp_expr_mis = expr_mis_t.norm_mis_by_gene;
nexpr_mis_t = result_table(result_table.is_notexpressed & result_table.mis_by_gene >= min_mut & ~isnan(result_table.norm_mis_by_gene) & ~isinf(result_table.norm_mis_by_gene),:);
temp_noexpr_mis = nexpr_mis_t.norm_mis_by_gene;
figure()
fig = gcf;
fig.InvertHardcopy='off';
fig.Units='inches';
fig.Position(3:4) = [4.5 3]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters'; 
[y,x]=ecdf(10.^temp_expr_mis);
semilogx(x,y,'-','LineWidth',3,'Color',[0.2627,0.5765,0.7647])
hold on
[y,x]=ecdf(10.^temp_noexpr_mis);
semilogx(x,y,'-','LineWidth',3,'Color',[0.8392,0.3765,0.3020])
ax=gca; 
ax.FontName=FontName; 
ax.FontSize=FontSize;
ax.LabelFontSizeMultiplier=1;
ax.TitleFontSizeMultiplier=1;
ax.TickDir='out';
ax.LineWidth=linewidth;
ax.Box='off';
ax.Color='none';
p = ranksum(temp_expr_mis,temp_noexpr_mis,'tail','both');
quant_results.P_expr_notexpr = p;
quant_results.FC_expr_notexpr = 1/((sum(expr_mis_t.mis_by_gene)/sum(expr_mis_t.length_AAs))...
    / (sum(nexpr_mis_t.mis_by_gene) / sum(nexpr_mis_t.length_AAs)));
quant_results.percentlost_expr_notexpr = 1- (((sum(expr_mis_t.mis_by_gene)/sum(expr_mis_t.length_AAs)))...
    / ((sum(nexpr_mis_t.mis_by_gene) / sum(nexpr_mis_t.length_AAs))));
quant_results.expr_nexpr_new_muts  = (sum(expr_mis_t.mis_by_gene)/tumor_number) * (1/(1-quant_results.percentlost_expr_notexpr));
if log10(p) <-20
    p=10^-20;
end

%ax.XLim=[min(10.^temp_expr_mis)/2,max(10.^temp_noexpr_mis)*2];
ax.YTick = ax.YLim(1):.5:ax.YLim(2);
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.YGrid='on';
ax.GridLineStyle='--';
xlabel('Missense mutation rate')
ylabel('Cumulative density of genes')
ax.TitleFontSizeMultiplier=1;
title(tumortype);
t=horzcat('p < 10^{',num2str(ceil(log10(p))),'}');
text(.1,0.9,t,'FontSize',FontSize,'FontName',FontName,'units','normalized')
L=legend({'Expressed','Silent'},'Location','southeast','FontSize',FontSize,'FontName',FontName);
L.Box = 'off';
L.FontSize = FontSize;
L.FontName = FontName;
fig.Color='none';
print(horzcat(mutfilename(1:4),'cdf.svg'),'-dsvg','-painters');
fig.Color=[1,1,1];

clear Pos ax binwidth p fc temp_expr_mis temp_noexpr_mis temp_expr_syn temp_noexpr_syn L t x y 

%% negatively-selected genes?
%muts/AA
%calculate the average muts/aa in expr- build expected based on gene length
result_table.exp_muts = result_table.length_AAs .* (sum(result_table.mis_by_gene) / sum(result_table.length_AAs));
result_table.obs_exp = result_table.mis_by_gene ./ result_table.exp_muts;
result_table.Pvaldn = cdf('Poisson',result_table.mis_by_gene,result_table.exp_muts);
result_table.Pvalup = 1-cdf('Poisson',result_table.mis_by_gene,result_table.exp_muts);
%result_table.is_neg_sel = ~result_table.is_notexpressed & result_table.Pvaldn <0.05;
%result_table.is_pos_sel = ~result_table.is_notexpressed & result_table.Pvalup <0.01;
result_table.is_neg_sel = result_table.obs_exp < .5;
result_table.is_pos_sel = result_table.obs_exp > 1.5;
result_table.is_pos_sel = result_table.obs_exp > 2;

%result_table.muts_per_aa = result_table.mis_by_gene ./ result_table.length_AAs;
%threshdn = mean(log10(result_table.muts_per_aa(result_table.is_notexpressed & result_table.mis_by_gene >0)))+log10(0.8);
%result_table.is_neg_sel = result_table.mis_by_gene>0 & log10(result_table.muts_per_aa)< threshdn & result_table.is_expressed;
%threshup = mean(log10(result_table.muts_per_aa(result_table.is_notexpressed & result_table.mis_by_gene >0))) + log10(2);
%result_table.is_pos_sel = result_table.mis_by_gene>0 & log10(result_table.muts_per_aa)> threshup & result_table.is_expressed;

%% evidence of selection in a.a. transitions?
%starting info
code = geneticcode(1);
codons = fieldnames(codoncount('atg'));
AAs={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
AAs = {'R' 'H' 'K' 'D' 'E' 'S' 'T' 'N' 'Q' 'C' 'G' 'P' 'A' 'V' 'I' 'L' 'M' 'F' 'Y' 'W'};

%filter mutations by missense and containing the relevant codon
%information
log_mut = mut.is_missense & mut.is_coding & mut.not_SNP &  strcmp(mut.Variant_Type,'SNP') & ismember(mut.Hugo_Symbol,result_table.mut_genes(sum(result_table.cdncount_by_gene(:,:),2)~=0 &~isnan(sum(result_table.cdncount_by_gene(:,:),2))));
log_mut = mut.is_missense & mut.is_coding &  strcmp(mut.Variant_Type,'SNP') & ismember(mut.Hugo_Symbol,result_table.mut_genes(sum(result_table.cdncount_by_gene(:,:),2)~=0 &~isnan(sum(result_table.cdncount_by_gene(:,:),2))));
mut_temp = mut(log_mut,:);

%let's compare our 3 different things: expression, essentiality, and
%selection
log_comparison = cell(4,1);
log_comparison_genes = cell(4,1);
comparison_labels = {'expr.','ess','sel','sel2'};

%comparing 

%comparing expressed and not expressed
log_comparison{1}  = [ismember(mut_temp.Hugo_Symbol,result_table.mut_genes(result_table.is_expressed & ~result_table.is_pos_sel)) ismember(mut_temp.Hugo_Symbol,result_table.mut_genes(result_table.is_notexpressed))];
log_comparison_genes{1} = [result_table.is_expressed & ~result_table.is_pos_sel & sum(result_table.cdncount_by_gene(:,:),2)~=0 &~isnan(sum(result_table.cdncount_by_gene(:,:),2)) result_table.is_notexpressed & sum(result_table.cdncount_by_gene(:,:),2)~=0 &~isnan(sum(result_table.cdncount_by_gene(:,:),2))];

%comparing essential and not essential
log_comparison{2}  = [ismember(mut_temp.Hugo_Symbol,result_table.mut_genes(result_table.is_essential)) ismember(mut_temp.Hugo_Symbol,result_table.mut_genes(result_table.is_notessential))];
log_comparison_genes{2} = [result_table.is_essential & sum(result_table.cdncount_by_gene(:,:),2)~=0 &~isnan(sum(result_table.cdncount_by_gene(:,:),2)) result_table.is_notessential & sum(result_table.cdncount_by_gene(:,:),2)~=0 &~isnan(sum(result_table.cdncount_by_gene(:,:),2))];

%comparing selected and not selected
log_comparison{3}  = [ismember(mut_temp.Hugo_Symbol,result_table.mut_genes(result_table.is_neg_sel)) ismember(mut_temp.Hugo_Symbol,result_table.mut_genes(~result_table.is_neg_sel & ~result_table.is_pos_sel))];
log_comparison_genes{3} = [result_table.is_neg_sel & sum(result_table.cdncount_by_gene(:,:),2)~=0 &~isnan(sum(result_table.cdncount_by_gene(:,:),2)) ~result_table.is_neg_sel & ~result_table.is_pos_sel & sum(result_table.cdncount_by_gene(:,:),2)~=0 &~isnan(sum(result_table.cdncount_by_gene(:,:),2))];

%comparing selected and not selected (expression-normalized);
%first compute expected:
result_table.logexpr =  log2(result_table.expr_by_tumor_geomean); % add expression score;
%bin by expression: 
bs = 0.25; %bin size;
mean_mutrate = zeros(height(result_table),1);
mean_expression = zeros(height(result_table),1);
flag = zeros(size(mean_mutrate));
binsize = zeros(size(mean_mutrate));
universe = result_table(:,:);
for b = 1:height(result_table); %for each gene
    ge = result_table.logexpr(b);
    gb = universe(universe.logexpr >= (ge - bs) & universe.logexpr <= (ge + bs),:); %pull genes close in expr
    mean_mutrate(b) = sum(gb.mis_by_gene) / sum(gb.length_AAs); %compute mutrate
    mean_expression(b) = mean(gb.logexpr); %compute average expression
    binsize(b) = height(gb);
    if height(gb) < 50;
        flag(b) = -1;
    end
end
result_table.expected = mean_mutrate .* result_table.length_AAs;
result_table.obs_exp = result_table.mis_by_gene ./ result_table.expected;
result_table.P_negsel = cdf('Poisson',result_table.mis_by_gene,result_table.expected);
result_table.is_neg_selected = result_table.obs_exp<=0.7 & result_table.is_expressed;
result_table.is_neg_selected2 = result_table.P_negsel < 0.01;
result_table.is_not_selected = result_table.obs_exp>0.9 & result_table.obs_exp < 1.1 & result_table.is_expressed;
log_comparison{4} =  [ismember(mut_temp.Hugo_Symbol,result_table.mut_genes(result_table.is_neg_selected))...
    ismember(mut_temp.Hugo_Symbol,result_table.mut_genes(result_table.is_not_selected))];
log_comparison_genes{4} = [result_table.is_neg_selected & sum(result_table.cdncount_by_gene(:,:),2)~=0 &~isnan(sum(result_table.cdncount_by_gene(:,:),2))... 
    result_table.is_not_selected & sum(result_table.cdncount_by_gene(:,:),2)~=0 &~isnan(sum(result_table.cdncount_by_gene(:,:),2))];
    
    
%now do this for both comparisons
%store result in change_res
change_res = cell(4,1);
for c=1:4
    
    %collect codon changes
    changes_in = cell(2,1);
    changes_in{1} = cellfun(@upper,cellfun(@(x) [x(end-6:end-4) x(end-2:end)],mut_temp.Codon_Change(log_comparison{c}(:,1)),'UniformOutput',false),'UniformOutput',false);
    changes_in{2} = cellfun(@upper,cellfun(@(x) [x(end-6:end-4) x(end-2:end)],mut_temp.Codon_Change(log_comparison{c}(:,2)),'UniformOutput',false),'UniformOutput',false);

    %tabulate codon changes 
    changes = cell(2,1);
    for n=1:2
        change =tabulate(changes_in{n});
        changes{n} = table(change(:,1),cell2mat(change(:,2)),'VariableNames',{'change','freq'},'RowNames',change(:,1));
    end

    %remove those changes that didn't occur in both fractions and store values in a table
    %change = intersect(changes{1}{changes{1}{:,2}>4,1},changes{2}{changes{2}{:,2}>4,1});
    %change = intersect(changes{1}{:,1},changes{2}{:,1});
    %change_t = join(changes{1}(change,:),changes{2}(change,:),'Keys','change');
    
    %combine changes from both sets into one table
    change_t = outerjoin(changes{1},changes{2},'Keys','change','MergeKeys',1);
    change_t.Properties.VariableNames = {'change','freq_e','freq_not_e'};
    
    %turn NaNs into 0s
    change_t{isnan(change_t.freq_e),'freq_e'} = 0;
    change_t{isnan(change_t.freq_not_e),'freq_not_e'} = 0;
    
    %find start, finish, and translations of changes
    change_t.start = cellfun(@(x) {x(1:3)},change_t.change);
    change_t.finish = cellfun(@(x) {x(4:6)},change_t.change);
    change_t.start_tr = cellfun(@(x) {code.(x)},change_t.start);
    change_t.finish_tr = cellfun(@(x) {code.(x)},change_t.finish);

    %count up codons in genes and add to table
    cdn_count = zeros(length(codons),2); 
    for n=1:2
        cdn_count(:,n) = sum(result_table.cdncount_by_gene(log_comparison_genes{c}(:,n),:),1);
    end
    cdn_count = array2table(cdn_count,'RowNames',codons,'VariableNames',{'cdn_count_e','cdn_count_not_e'});
    change_t.cdn_count_e = cdn_count{change_t.start,'cdn_count_e'};
    change_t.cdn_count_not_e = cdn_count{change_t.start,'cdn_count_not_e'};
    
    %which codon transitions (that move between AAs) are possible?
    pos_t_cdn = zeros(length(codons));
    pos_t_aa = zeros(length(AAs));
    %which codons can go to which amino acid?
    pos_t_cdn_aa = zeros(length(codons),length(AAs));
    nt = {'A','T','G','C'}; 
    for x=1:length(codons)
        cdn = codons{x};
        start = code.(cdn);
        if start == '*'
            continue
        end
        %mutate bp at each position and find new aa
        for i=1:3;
            for t = 1:4
                new_codon = cdn;
                new_codon(i)=nt{t};
                finish = code.(new_codon);
                if finish == '*' || finish == start
                    continue
                end
                pos_t_cdn(strcmpi(codons,cdn), strcmpi(codons,new_codon)) = 1;
                pos_t_aa(strcmpi(AAs,start),strcmpi(AAs,finish)) = 1;
                pos_t_cdn_aa(strcmpi(codons,cdn),strcmpi(AAs,finish)) =1 ;
            end
        end
    end
    pos_t_cdn_aa = array2table(pos_t_cdn_aa,'RowNames',codons,'VariableNames',AAs);
    pos_t_cdn_aa.st_trans = cellfun(@(x) {code.(x)},pos_t_cdn_aa.Properties.RowNames);
    %make a list of possible amino acid transitions
    change = cell(1);
    i=1;
    for x = 1:length(AAs)
        for y = 1:length(AAs)
            if pos_t_aa(x,y) == 1
                change(i,1) = {[AAs{x} AAs{y}]};
                i = i+1;
            end
        end
    end
    
    
    %now for each kind of aa---> aa transition, add up codons that could, and times
    %it did
    change_t.tr = cellfun(@(x,y) {[x y]},change_t.start_tr,change_t.finish_tr);
    %change = unique(change_t.tr);
    change_a = [array2table(cell(length(change),3),'VariableNames',{'start','finish','change'}) array2table(zeros(length(change),4),'VariableNames',{'freq_e','freq_not_e','cdn_count_e','cdn_count_not_e'})];
    for n=1:length(change)
        temp_t = change_t(ismember(change_t.tr,change{n}),:);
        if isequal(height(temp_t),0)
            change_a{n,'change'} = change(n);
            change_a{n,'freq_e'} = 0;
            change_a{n,'freq_not_e'} = 0;
            change_a{n,'start'} = {change{n}(1)};
            change_a{n,'finish'} = {change{n}(2)};
        else
            change_a{n,'change'} = temp_t.tr(1);
            %count up mutations
            change_a{n,'freq_e'} = sum(temp_t.freq_e);
            change_a{n,'freq_not_e'} = sum(temp_t.freq_not_e);
            change_a{n,'start'} = temp_t.start_tr(1);
            change_a{n,'finish'} = temp_t.finish_tr(1);
        end
        %get unique codons in transition and count up
        %temp_t_cdn = unique(temp_t(:,{'start','cdn_count_e','cdn_count_not_e'}),'rows');
        %change_a{n,'cdn_count_e'} = sum(temp_t_cdn.cdn_count_e);
        %change_a{n,'cdn_count_not_e'} = sum(temp_t_cdn.cdn_count_not_e);
        %which codons for st could go to the ending aa? sum them up.
        temp_t_cdn = cdn_count(strcmpi(pos_t_cdn_aa.st_trans,change{n}(1)) & logical(pos_t_cdn_aa.(change{n}(2))),:);
        change_a{n,'cdn_count_e'} = sum(temp_t_cdn.cdn_count_e);
        change_a{n,'cdn_count_not_e'} = sum(temp_t_cdn.cdn_count_not_e);
        
    end
    
   

    %calculate fold change in mutation frequency by aa (selected / not-selected)
    change_a.fc = (change_a.freq_e./change_a.cdn_count_e) ./ (change_a.freq_not_e./change_a.cdn_count_not_e);
    %pval = ones(height(change_a),1);
    %for n=1:height(change_a)
     %   pval(n) = chisquarecont([[change_a.freq_e(n); change_a.cdn_count_e(n) - change_a.freq_e(n)] [change_a.freq_not_e(n); change_a.cdn_count_not_e(n) - change_a.freq_not_e(n)]]);
    %end
    %change_a.pval = pval;
    
    %calculate the expected change--- fold change all mutation frequencies(not-selected / selected)
    R = (sum(change_a.freq_not_e)/sum(change_a.cdn_count_not_e)) /(sum(change_a.freq_e)/sum(change_a.cdn_count_e));
    %what is CDN count for 2 groups?
    cdncnt = zeros(2,1);
    for n=1:2
        for i=1:length(codons)
            if ~strcmp(code.(codons{i}),'*')
                cdncnt(n,1) = cdncnt(n,1) + cdn_count{i,n};
            end
        end
    end
    cdncnt = array2table(cdncnt','VariableNames',{'cdn_count_e','cdn_count_not_e'});
    R = (sum(change_a.freq_not_e)/sum(change_a.cdn_count_not_e)) /(sum(change_a.freq_e)/sum(change_a.cdn_count_e));
    
    
    
    %calculate a X2 value for each transition
    
    X2=zeros(height(change_a),1);
    for i=1:height(change_a)
        A = change_a.freq_e(i);
        B = change_a.freq_not_e(i);
        C = change_a.cdn_count_e(i);
        D = change_a.cdn_count_not_e(i);
        % calculate expected for each change
        exp = zeros(2,1);
        exp(1) = ((A*R + B)/R) * (C / (C + D));
        exp(2) = (A*R + B) * (D /(D + C));
        
        X2(i) = ((A-exp(1))^2)/exp(1) + ((B-exp(2))^2)/exp(2);
    end
    change_a.X2 = X2;
    clear A B C D exp X2;
    
    %calculate P val
    change_a.Pval = 1-chi2cdf(change_a.X2,1);
        
    
            
    
    %calculate a normalized fold change (log2 fold change / mean fold change)
    %change_a.norm_fc = log2(change_a.fc / mean(change_a.fc));
    % change_a.norm_fc = (change_a.fc / mean(change_a.fc));
    %store result in change_res
    change_res{c} = change_a;
    
    %only do the following for expression (more mutations);
    switch c
        case 1
            
            %find all of the bidirectional changes
            change_a.change_rev = cellfun(@(x,y) {[y x]},change_a.start,change_a.finish);
            bidir = change_a(ismember(change_a.change,change_a.change_rev),:);
            %grab pvals and fold changes from reverse 
            bitemp = bidir;
            bitemp.Properties.RowNames = bitemp.change_rev;
            bidir.rev_fc = bitemp{bitemp.change,'fc'};
            bidir.rev_Pval = bitemp{bitemp.change,'Pval'};
            
            %score similar amino acids as ones you can easily go forward and back from
            P_thresh = .251;
            bidir.sym_good = bidir.fc >= (1/R) & bidir.rev_fc >= (1/R) & bidir.Pval <= P_thresh & bidir.rev_Pval <= P_thresh;
            bidir.sym_good = bidir.fc >= (1/R) & bidir.rev_fc >= (1/R) & (bidir.Pval + bidir.rev_Pval)./2 < 0.25;
            

            %record data in adjacency matrix
            pair_sym = array2table(zeros(20,20),'RowNames',AAs,'VariableNames',AAs);
            temp_t = bidir(bidir.sym_good,:);
            for t = 1:height(temp_t)
                pair_sym{temp_t.start(t),temp_t.finish(t)} = 1;
            end
            %Plot similar amino acids connected in network graph
            figure()
            fig = gcf;
            fig.InvertHardcopy='off';
            fig.PaperPositionMode='Auto';
            fig.Renderer='painters';
            G = graph(pair_sym{:,:},AAs);
            P=plot(G,'Layout','force','Iterations',40);
            %P.EdgeCData=G.Edges.Weight;
            P.EdgeAlpha=1;
            P.LineWidth=2;
            P.MarkerSize=5;
            P.NodeLabel = aminolookup(P.NodeLabel);
            P.EdgeColor = 'c';
            highlight(P,{'H' 'K' 'R'},'nodecolor','r')
            highlight(P,{'L' 'V' 'I' 'M' 'A'},'nodecolor','g')
            highlight(P,{'C' 'G' 'P'},'nodecolor','y')
            highlight(P,{'Y' 'F' 'W'},'nodecolor',[ 0.9100 0.4100 0.1700])
            highlight(P,{'N' 'Q' 'T' 'S'},'nodecolor',[0.5 0 0.9])
            highlight(P,{'E' 'D' },'nodecolor','b')
            %P.NodeColor = [.3 .3 .3];
            ax=gca;
            ax.XTick = [];
            ax.YTick = [];
            ax.FontName=FontName;
            ax.FontSize=FontSize;
            ax.TitleFontSizeMultiplier=1;
            title(horzcat(tumortype,' (',comparison_labels{c},')'))
            ax.LineWidth=1;
            ax.Color='none';
            fig.Color='none';
            print(horzcat(mutfilename(1:4),comparison_labels{c},'_aa_selection_netgraph.svg'),'-dsvg','-painters');
            fig.Color=[1,1,1];
            
            %make a similar plot for all transitions
            %bidir_t = bidir;
            %drop if both directions are zero
            %bidir_t(bidir_t.freq_e == 0 & bidir_t.freq_not_e == 0,:)=[];
            %remove repeats (only one line per bidirectional transition)
            %bidir_log = ones(size(bidir_t.change));
            %for n=1:height(bidir_t)
            %    if bidir_log(n)==1
            %        bidir_log(strcmp(bidir_t.change(n),bidir_t.change_rev))=0;
            %    end
            %end
            %bidir_t1 = bidir_t(logical(bidir_log),:);
            
            %make a similar plot for all possible transitions
            %start with a 400 long table of all possible changes

            pos_t = zeros(length(AAs));
            %what aa transitions are even possible by the codon code?
            nt = {'A','T','G','C'}; 
            for x=1:length(codons)
                cdn = codons{x};
                start = code.(cdn);
                if start == '*'
                    continue
                end
                %mutate bp at each position and find new aa
                for i=1:3;
                    for t = 1:4
                        new_codon = cdn;
                        new_codon(i)=nt{t};
                        finish = code.(new_codon);
                        if finish == '*' || finish == start
                            continue
                        end
                        pos_t(strcmpi(AAs,start), strcmpi(AAs,finish)) = 1;
                        
                    end
                end
            end
            %drop synonymous
            for x=1:20
                pos_t(x,x) = 0;
            end
            
            figure()
            fig = gcf;
            fig.InvertHardcopy='off';
            fig.PaperPositionMode='Auto';
            fig.Renderer='painters';
            G_all = graph(pos_t,AAs);
            P=plot(G_all,'Layout','force','Iterations',40);
            highlight(P,{'H' 'K' 'R'},'nodecolor','r')
            highlight(P,{'L' 'V' 'I' 'M' 'A'},'nodecolor','g')
            highlight(P,{'C' 'G' 'P'},'nodecolor','y')
            highlight(P,{'Y' 'F' 'W'},'nodecolor',[ 0.9100 0.4100 0.1700])
            highlight(P,{'N' 'Q' 'T' 'S'},'nodecolor',[0.5 0 0.9])
            highlight(P,{'E' 'D' },'nodecolor','b')
            P.NodeLabel = aminolookup(P.NodeLabel);
            P.EdgeColor = 'c';
            P.MarkerSize=5;
            P.EdgeAlpha=1;
            P.LineWidth=1;
            ax=gca;
            ax.Box = 'off';
            ax.Visible = 'off';
            ax.XTick = [];
            ax.YTick = [];
            ax.LineWidth=1;
            ax.Color='none';
            ax.FontName=FontName;
            ax.FontSize=FontSize;
            fig.Color='none';
            print(horzcat(mutfilename(1:4),comparison_labels{c},'_aa_all_netgraph.svg'),'-dsvg','-painters');
            fig.Color=[1,1,1];
            
            %make a similar plot for Blosum 90 conserved substitutions
            blos = blosum(90,'Order',strjoin(AAs,''));
            
            figure()
            fig = gcf;
            fig.InvertHardcopy='off';
            fig.PaperPositionMode='Auto';
            fig.Renderer='painters';
            G_blos = graph(blos>=0,AAs,'OmitSelfLoops');
            P=plot(G_blos,'Layout','force','Iterations',40);
            highlight(P,{'H' 'K' 'R'},'nodecolor','r')
            highlight(P,{'L' 'V' 'I' 'M' 'A'},'nodecolor','g')
            highlight(P,{'C' 'G' 'P'},'nodecolor','y')
            highlight(P,{'Y' 'F' 'W'},'nodecolor',[ 0.9100 0.4100 0.1700])
            highlight(P,{'N' 'Q' 'T' 'S'},'nodecolor',[0.5 0 0.9])
            highlight(P,{'E' 'D' },'nodecolor','b')
            P.NodeLabel = aminolookup(P.NodeLabel);
            P.EdgeColor = 'c';
            P.MarkerSize=5;
            P.EdgeAlpha=1;
            P.LineWidth=2;
            ax=gca;
            ax.Box = 'off';
            ax.Visible = 'off';
            ax.XTick = [];
            ax.YTick = [];
            ax.LineWidth=1;
            ax.Color='none';
            ax.FontName=FontName;
            ax.FontSize=FontSize;
            fig.Color='none';
            print(horzcat(mutfilename(1:4),comparison_labels{c},'_aa_blosum_netgraph.svg'),'-dsvg','-painters');
            fig.Color=[1,1,1];
            
            %now plot the intersection between all IDed substitutions and
            %blosum-tolerated substitutions
            blos_edge = cellfun(@(y,x) {[y x]},G_blos.Edges{:,1}(:,1),G_blos.Edges{:,1}(:,2));
            blos_pos = rmedge(G_blos,find(~ismember(blos_edge,bidir.change) & ~ismember(blos_edge,bidir.change_rev)));
            figure()
            fig = gcf;
            fig.InvertHardcopy='off';
            fig.PaperPositionMode='Auto';
            fig.Renderer='painters';
            P=plot(blos_pos,'Layout','force','Iterations',40);
            highlight(P,{'H' 'K' 'R'},'nodecolor','r')
            highlight(P,{'L' 'V' 'I' 'M' 'A'},'nodecolor','g')
            highlight(P,{'C' 'G' 'P'},'nodecolor','y')
            highlight(P,{'Y' 'F' 'W'},'nodecolor',[ 0.9100 0.4100 0.1700])
            highlight(P,{'N' 'Q' 'T' 'S'},'nodecolor',[0.5 0 0.9])
            highlight(P,{'E' 'D' },'nodecolor','b')
            P.NodeLabel = aminolookup(P.NodeLabel);
            P.EdgeColor = 'c';
            P.MarkerSize=5;
            P.EdgeAlpha=1;
            P.LineWidth=2;
            ax=gca;
            ax.Box = 'off';
            ax.Visible = 'off';
            ax.XTick = [];
            ax.YTick = [];
            ax.LineWidth=1;
            ax.Color='none';
            ax.FontName=FontName;
            ax.FontSize=FontSize;
            fig.Color='none';
            print(horzcat(mutfilename(1:4),comparison_labels{c},'_aa_blosum_pos_netgraph.svg'),'-dsvg','-painters');
            fig.Color=[1,1,1];
            
            %calculate a pvalue of the overlap here
            blos_edgs = blos_pos.Edges.EndNodes;
            blos_edgs = cellfun(@(x,y) {[y x]},blos_edgs(:,1),blos_edgs(:,2));
            x = nnz(ismember(blos_edgs,bidir.change(bidir.sym_good))); %overlap
            M = height(bidir)/2;  %size of population
            K = length(blos_edgs); %geneset size
            N = nnz(bidir.sym_good)/2; %genes drawn
            quant_results.ConservativeAAsPval = hygecdf(x,M,K,N,'upper')*2; %hypergeometric pval
            
            %Make a scatter-spy showing observed and
            %symmetricly-ok transitions 
                       
            spy_t = change_a(:,{'change'});
            st = zeros(height(spy_t),1);
            fn = zeros(height(spy_t),1);
            for s = 1:height(spy_t)
                st(s)= find(strcmp(AAs,change_a.start{s}));
                fn(s) = find(strcmp(AAs,change_a.finish{s}));
            end
            spy_t.start = st;
            spy_t.finish = fn;
            spy_t.Properties.RowNames = spy_t.change;
            %observed at least n times?
            spy_t.observed = zeros(height(spy_t),1);
            spy_t{change_a.change(change_a.freq_e > 0 & change_a.freq_not_e >= 0),'observed'} = 1;
            
            %symmetric above average (ie less selected than average in both directions?)
            spy_t.sym_good = zeros(height(spy_t),1);
            spy_t{bidir{bidir.sym_good,'change'},'sym_good'} = 1;
                                  
            %drop out not observed enough
            spy_t = spy_t(logical(spy_t.observed),:);
            
            %add in fc score
            change_a.Properties.RowNames = change_a.change;
            %spy_t.norm_fc = ((change_a{spy_t.Properties.RowNames,'norm_fc'})-max(change_a.norm_fc))*-1;
            %spy_t.fc = ((change_a{spy_t.Properties.RowNames,'fc'}).^-1);
            spy_t.fc = -1*(log2(change_a{spy_t.Properties.RowNames,'fc'}));
                        
            figure()
            fig = gcf;
            fig.InvertHardcopy='off';
            fig.PaperPositionMode='Auto';
            fig.Renderer='painters';
            fig.Units = 'inches'; 
            fig.Position(3:4) = [5 4]; 
            scatter(spy_t.start(logical(spy_t.sym_good)),spy_t.finish(logical(spy_t.sym_good)),150,[0 0 0],'filled');
            hold on
            scatter(spy_t.start,spy_t.finish,70,spy_t.fc,'filled')

            plot(1:20,1:20,'--k');
            ax = gca;
            ax.XLim = [0 21];
            ax.YLim = [0 21];
            ax.XTick = 1:20;
            ax.YTick = 1:20;
            ax.LineWidth = linewidth;
            ax.XTickLabel = AAs;
            ax.YTickLabel = AAs;
            ax.Color = 'None';
            ax.FontName = FontName; 
            ax.FontSize = FontSize; 
            ax.LabelFontSizeMultiplier = 1; 
            xlabel('Wild-type a.a.')
            ylabel('Mutant a.a.')
            C = colorbar;
            colormap(cool(7));
            %colormap(flipud(colormap));
            C.FontSize = FontSize;
            C.FontName = FontName;
            fig.Color = 'None';
            print(horzcat(mutfilename(1:4),comparison_labels{c},'_aa_selection_spyplot.svg'),'-dsvg','-painters');
            fig.Color = [1 1 1];
       
            
    end
end




%% Quantifying mutation rates of different base pairs in the different strands
%%and plotting expressed / not expressed mutation rates
nt = {'T','G','A','C'}; 
codons = fieldnames(codoncount('atg'));
switch tumortype
    case 'Melanoma'
        tr = {'TA','TC','TG','CA','CT','CG'};
    case 'Lung adenocarcinoma';
        tr = {'AT','AG','AC','GT','GA','GC'};
        
    case 'Colorectal adenocarcinoma'
        tr = {'AT','AG','AC','GT','GA','GC'};
end
expr_nexpr_muttype = zeros(length(tr),2);
done = zeros(size(expr_nexpr_muttype,1),1);
label = cell(size(expr_nexpr_muttype,1),1);      
log_mut = mut.is_missense & mut.is_coding & strcmp(mut.Variant_Type,'SNP');
mut_temp = mut(log_mut,:);
mut_temp.transcript_strand = strcmp(mut_temp.Transcript_Strand,'+');
i=0;
mut_nums = zeros(length(tr),1);
for t = 1:length(tr)
    st = tr{t}(1);
    fn = tr{t}(2);
    log_transcribed = (strcmp(mut_temp.Reference_Allele,st) & strcmp(mut_temp.Tumor_Seq_Allele2,fn) & ~mut_temp.transcript_strand)...
        | (strcmp(mut_temp.Reference_Allele,seqcomplement(st)) & strcmp(mut_temp.Tumor_Seq_Allele2,seqcomplement(fn)) & mut_temp.transcript_strand);
    log_nottranscribed = (strcmp(mut_temp.Reference_Allele,st) & strcmp(mut_temp.Tumor_Seq_Allele2,fn) & mut_temp.transcript_strand)...
        | (strcmp(mut_temp.Reference_Allele,seqcomplement(st)) & strcmp(mut_temp.Tumor_Seq_Allele2,seqcomplement(fn)) & ~mut_temp.transcript_strand);
    mut_nums(t) = nnz(log_transcribed)+nnz(log_nottranscribed);
    mis1 = tabulate(mut_temp.Hugo_Symbol(log_transcribed));
    transcribed = table(mis1(:,1),cell2mat(mis1(:,2)),'VariableNames',{'genenames','mutations'});
    mis1 = tabulate(mut_temp.Hugo_Symbol(log_nottranscribed));
    nottranscribed = table(mis1(:,1),cell2mat(mis1(:,2)),'VariableNames',{'genenames','mutations'});
    mut_by_gene = outerjoin(transcribed,nottranscribed,'key','genenames','MergeKeys',1);
    mut_by_gene{isnan(mut_by_gene.mutations_transcribed),'mutations_transcribed'} = 0;
    mut_by_gene{isnan(mut_by_gene.mutations_nottranscribed),'mutations_nottranscribed'} = 0;
    %append genes with 0 mutations
    zeromut_g = result_table.mut_genes(~ismember(result_table.mut_genes,mut_by_gene.genenames));
    zeromuts = array2table(zeromut_g,'VariableNames',{'genenames'});
    zeromuts.mutations_transcribed = zeros(height(zeromuts),1);
    zeromuts.mutations_nottranscribed = zeros(height(zeromuts),1);
    mut_by_gene = vertcat(mut_by_gene,zeromuts);
    %count up nucleotides in sequences in coding strand,(len(genes) x 4);
    cdncounts =readtable('codon_counttable_1.txt','delimiter','\t');
    cdncounts = cdncounts(ismember(cdncounts.Row,mut_by_gene.genenames),:);
    cdncount = table(table2array(cdncounts(:,2:end)),'VariableNames',{'cdncount_by_gene'});
    cdncount.Row = cdncounts.Row;
    mut_by_gene = outerjoin(mut_by_gene,cdncount,'LeftKey','genenames','RightKey',{'Row'},'MergeKeys',0);
    cdncount = mut_by_gene.cdncount_by_gene;
    nt_count = zeros(size(mut_by_gene,1),length(nt));
    for n = 1:length(nt) %for each nucleotide
        for c = 1:length(codons) % for each codon
            ntc = sum(codons{c} == nt{n});
            for g = 1:size(nt_count,1) % for each gene
                nt_count(g,n) = nt_count(g,n) + (ntc * cdncount(g,c));           
            end
        end
    end
    mut_by_gene.start_transc = nt_count(:,strcmp(nt,seqcomplement(st)));
    mut_by_gene.start_nottransc = nt_count(:,strcmp(nt,st));
    mut_by_gene.is_expressed = ismember(mut_by_gene.genenames,result_table.mut_genes(result_table.is_expressed));
    mut_by_gene.is_notexpressed = ismember(mut_by_gene.genenames,result_table.mut_genes(result_table.is_notexpressed));
    mut_by_gene_t = mut_by_gene(~isnan(mut_by_gene.start_transc),:);
    MR_expr_transc = sum(mut_by_gene_t.mutations_transcribed(mut_by_gene_t.is_expressed)) / sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_expressed));
    MR_expr_ntransc = sum(mut_by_gene_t.mutations_nottranscribed(mut_by_gene_t.is_expressed)) / sum(mut_by_gene_t.start_nottransc(mut_by_gene_t.is_expressed));
    MR_nexpr_transc = sum(mut_by_gene_t.mutations_transcribed(mut_by_gene_t.is_notexpressed)) / sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_notexpressed));
    MR_nexpr_ntransc = sum(mut_by_gene_t.mutations_nottranscribed(mut_by_gene_t.is_notexpressed)) / sum(mut_by_gene_t.start_nottransc(mut_by_gene_t.is_notexpressed));

    %What is difference in mutation rates?
    expr_nexpr_muttype(t,1) = MR_expr_transc / MR_nexpr_transc;
    expr_nexpr_muttype(t,2)= MR_expr_ntransc / MR_nexpr_ntransc;
    %frac_transc_st = sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_expressed)) / (sum(mut_by_gene_t.start_nottransc(mut_by_gene_t.is_expressed)) + sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_expressed)));

    label{t,1} = horzcat(st,'>',fn);
end

figure()
fig = gcf;
fig.InvertHardcopy='off';
fig.Units='inches';
fig.Position(3:4) = [4.5 3]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters'; 
b=bar(expr_nexpr_muttype,1.2,'EdgeColor','black');
b(1).FaceColor = [1,0,0];
b(2).FaceColor = [0,0,0];
L=legend({'transcribed strand','other strand'});
L.Box = 'off'; L.FontSize = FontSize; L.FontName = FontName;
ylim([0,1]);ylabel({'Expressed / Not Expressed','mutation rate'});
title(tumortype)
ax=gca; ax.XTickLabel = label; ax.FontName=FontName; ax.FontSize=FontSize; ax.LabelFontSizeMultiplier=1; ax.TitleFontSizeMultiplier=1;
ax.TickDir='out'; ax.LineWidth=linewidth; ax.Box='off'; ax.Color='none'; fig.Color='none'; 
print(horzcat(mutfilename(1:4),'bar_strand.svg'),'-dsvg','-painters');
fig.Color = [1 1 1];

%% Quantifying negative selection and transcription coupled repair

codons = fieldnames(codoncount('atg'));
nt = {'A','T','G','C'}; 

%bin genes by expressed & neg selected, not-expressed
log_mut = mut.is_missense & mut.is_coding & strcmp(mut.Variant_Type,'SNP');
mut_temp = mut(log_mut,:);
mut_temp.transcript_strand = strcmp(mut_temp.Transcript_Strand,'+');

%look at different muts in each tumor type; count up mutations in
%transcribed and not transcribed strand
switch tumortype
    case 'Lung adenocarcinoma'
        %look in particular at C--> A / G-->T mutations.  
        st = 'G';
        fn = 'T';
        
    case 'Melanoma'
        %look in particular at C--> T / G-->A mutations.
        st = 'C';
        fn = 'T'; 
    case 'Colorectal adenocarcinoma'
        st = 'A';
        fn = 'C';
end
log_transcribed = (strcmp(mut_temp.Reference_Allele,st) & strcmp(mut_temp.Tumor_Seq_Allele2,fn) & ~mut_temp.transcript_strand) |...
    (strcmp(mut_temp.Reference_Allele,seqcomplement(st)) & strcmp(mut_temp.Tumor_Seq_Allele2,seqcomplement(fn)) & mut_temp.transcript_strand);
log_nottranscribed = (strcmp(mut_temp.Reference_Allele,st) & strcmp(mut_temp.Tumor_Seq_Allele2,fn) & mut_temp.transcript_strand) |...
    (strcmp(mut_temp.Reference_Allele,seqcomplement(st)) & strcmp(mut_temp.Tumor_Seq_Allele2,seqcomplement(fn)) & ~mut_temp.transcript_strand);
        
%for each gene in expressed and not expressed, count up mutations in these two strands
mis1 = tabulate(mut_temp.Hugo_Symbol(log_transcribed));
transcribed = table(mis1(:,1),cell2mat(mis1(:,2)),'VariableNames',{'genenames','mutations'});
mis1 = tabulate(mut_temp.Hugo_Symbol(log_nottranscribed));
nottranscribed = table(mis1(:,1),cell2mat(mis1(:,2)),'VariableNames',{'genenames','mutations'});
mut_by_gene = outerjoin(transcribed,nottranscribed,'key','genenames','MergeKeys',1);
mut_by_gene{isnan(mut_by_gene.mutations_transcribed),'mutations_transcribed'} = 0;
mut_by_gene{isnan(mut_by_gene.mutations_nottranscribed),'mutations_nottranscribed'} = 0;

%load in gene lengths
%genelengths =readtable('uniprot-len.xlsx','ReadRowNames',false,'ReadVariableNames',1);
%genelengths.Properties.VariableNames = {'genenames','length_AAs'};
%mut_by_gene = innerjoin(mut_by_gene,genelengths,'key','genenames');

%append genes with 0 mutations
zeromut_g = result_table.mut_genes(~ismember(result_table.mut_genes,mut_by_gene.genenames));
zeromuts = array2table(zeromut_g,'VariableNames',{'genenames'});
zeromuts.mutations_transcribed = zeros(height(zeromuts),1);
zeromuts.mutations_nottranscribed = zeros(height(zeromuts),1);
mut_by_gene = vertcat(mut_by_gene,zeromuts);


%count up nucleotides in sequences in coding strand,(len(genes) x 4);
cdncounts =readtable('codon_counttable_1.txt','delimiter','\t');
cdncounts = cdncounts(ismember(cdncounts.Row,mut_by_gene.genenames),:);
cdncount = table(table2array(cdncounts(:,2:end)),'VariableNames',{'cdncount_by_gene'});
cdncount.Row = cdncounts.Row;
mut_by_gene = outerjoin(mut_by_gene,cdncount,'LeftKey','genenames','RightKey',{'Row'},'MergeKeys',0);

%cdncount = array2table(result_table.cdncount_by_gene,'RowNames',result_table.mut_genes);
%cdncount = cdncount(mut_by_gene.genenames,:);
cdncount = mut_by_gene.cdncount_by_gene;
nt_count = zeros(size(mut_by_gene,1),length(nt));
for n = 1:length(nt) %for each nucleotide
    for c = 1:length(codons) % for each codon
        ntc = sum(codons{c} == nt{n});
        for g = 1:size(nt_count,1) % for each gene
            nt_count(g,n) = nt_count(g,n) + (ntc * cdncount(g,c));           
        end
    end
end
mut_by_gene.start_transc = nt_count(:,strcmp(nt,seqcomplement(st)));
mut_by_gene.start_nottransc = nt_count(:,strcmp(nt,st));

%calculate normalized mutation rates
mutrate_by_gene = array2table(zeros(height(mut_by_gene),2),'VariableNames',{'mutrate_transcribed','mutrate_nottranscribed'});
mutrate_by_gene.mutrate_transcribed = log10((1 + mut_by_gene.mutations_transcribed) ./ mut_by_gene.start_transc  / (tumor_number / scale_factor));
mutrate_by_gene.mutrate_nottranscribed = log10((1 + mut_by_gene.mutations_nottranscribed) ./ mut_by_gene.start_nottransc  / (tumor_number / scale_factor));
mutrate_by_gene.genenames = mut_by_gene.genenames;
mutrate_by_gene.start_transc = mut_by_gene.start_transc;
mutrate_by_gene.start_nottransc = mut_by_gene.start_transc;
mutrate_by_gene_t = mutrate_by_gene(~isnan(mutrate_by_gene.mutrate_transcribed) & ~isnan(mutrate_by_gene.mutrate_nottranscribed) & mutrate_by_gene.start_transc > 0,:);
mutrate_by_gene_t.is_expressed = ismember(mutrate_by_gene_t.genenames,result_table.mut_genes(result_table.is_expressed));
mutrate_by_gene_t.is_notexpressed = ismember(mutrate_by_gene_t.genenames,result_table.mut_genes(result_table.is_notexpressed));

temp_expr_mis = mutrate_by_gene_t{mutrate_by_gene_t.is_expressed,'mutrate_nottranscribed'};
temp_noexpr_mis = mutrate_by_gene_t{mutrate_by_gene_t.is_notexpressed,'mutrate_nottranscribed'};
temp_expr_mis_transc = mutrate_by_gene_t{mutrate_by_gene_t.is_expressed,'mutrate_transcribed'};
temp_noexpr_mis_transc = mutrate_by_gene_t{mutrate_by_gene_t.is_notexpressed,'mutrate_transcribed'};

figure()
fig = gcf;
fig.InvertHardcopy='off';
fig.Units='inches';
fig.Position(3:4) = [4.5 3]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters'; 

[y,x]=ecdf(10.^temp_expr_mis_transc);
semilogx(x,y,'-','LineWidth',3,'Color',[0,1,0])
hold on

[y,x]=ecdf(10.^temp_expr_mis);
semilogx(x,y,'-','LineWidth',3,'Color',[.3,0.8,.7])

[y,x]=ecdf(10.^temp_noexpr_mis_transc);
semilogx(x,y,'-','LineWidth',3,'Color',[1,0,1])

[y,x]=ecdf(10.^temp_noexpr_mis);
semilogx(x,y,'-','LineWidth',3,'Color',[0.7,0.3,0.7])

ax=gca; 
ax.FontName=FontName; 
ax.FontSize=FontSize;
ax.LabelFontSizeMultiplier=1;
ax.TitleFontSizeMultiplier=1;
ax.TickDir='out';
ax.LineWidth=linewidth;
ax.Box='off';
ax.Color='none';
p = ranksum(temp_expr_mis,temp_noexpr_mis,'tail','both');
quant_results.P_expr_notexpr = p;
quant_results.FC_expr_notexpr = 1/(10^(median(temp_expr_mis) - median(temp_noexpr_mis)));
if log10(p) <-20
    p=10^-20;
end
ax.YTick = ax.YLim(1):.5:ax.YLim(2);
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.YGrid='on';
ax.GridLineStyle='--';
xlabel('Missense mutation rate')
ylabel('Cumulative density of genes')
ax.TitleFontSizeMultiplier=1;
title(tumortype);
t=horzcat('p < 10^{',num2str(ceil(log10(p))),'}');
text(.1,0.9,t,'FontSize',FontSize,'FontName',FontName,'units','normalized')
L=legend({'Expr. genes, transc. strand',...
    'Expr. genes, not transc. strand',...
    'Not expr. genes, transc. strand',...
    'Not expr. genes, not transc. strand'},'Location','southeast','FontSize',FontSize,'FontName',FontName);
L.Box = 'off';
L.FontSize = FontSize;
L.FontName = FontName;
fig.Color='none';
print(horzcat(mutfilename(1:4),'cdf_strand.svg'),'-dsvg','-painters');
fig.Color=[1,1,1];

%Now, what fraction of mutations are missing?
%first, calculate the mutation rates of the different groups of genes
mut_by_gene.is_expressed = ismember(mut_by_gene.genenames,result_table.mut_genes(result_table.is_expressed));
mut_by_gene.is_notexpressed = ismember(mut_by_gene.genenames,result_table.mut_genes(result_table.is_notexpressed));
mut_by_gene_t = mut_by_gene(~isnan(mut_by_gene.start_transc),:);
MR_expr_transc = sum(mut_by_gene_t.mutations_transcribed(mut_by_gene_t.is_expressed)) / sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_expressed));
MR_expr_ntransc = sum(mut_by_gene_t.mutations_nottranscribed(mut_by_gene_t.is_expressed)) / sum(mut_by_gene_t.start_nottransc(mut_by_gene_t.is_expressed));
MR_nexpr_transc = sum(mut_by_gene_t.mutations_transcribed(mut_by_gene_t.is_notexpressed)) / sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_notexpressed));
MR_nexpr_ntransc = sum(mut_by_gene_t.mutations_nottranscribed(mut_by_gene_t.is_notexpressed)) / sum(mut_by_gene_t.start_nottransc(mut_by_gene_t.is_notexpressed));

%What is difference in mutation rates?
MD_transc = MR_expr_transc / MR_nexpr_transc;
MD_ntransc = MR_expr_ntransc / MR_nexpr_ntransc;

frac_transc_st = sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_expressed)) / (sum(mut_by_gene_t.start_nottransc(mut_by_gene_t.is_expressed)) + sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_expressed)));
frac_transc_st_nexpr = sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_notexpressed)) / (sum(mut_by_gene_t.start_nottransc(mut_by_gene_t.is_notexpressed)) + sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_notexpressed)));

% observed is average of expressed strand mutation rates
observed = frac_transc_st*MD_transc + (1-frac_transc_st)*MD_ntransc; % weight by #Gs
observed_MR = frac_transc_st*MR_expr_transc + (1-frac_transc_st)*MR_expr_ntransc; % weight by #Gs
expected_MR = frac_transc_st_nexpr*MR_nexpr_transc + (1-frac_transc_st_nexpr)*MR_nexpr_ntransc; 
noTCR_observed_MR = MR_expr_ntransc;
MD_obs = 1-observed_MR/expected_MR;
MD_noTCR = 1-noTCR_observed_MR/expected_MR;

% TCR is difference in MD_transc and MD_ntransc; only works on transc
% strand
TCR = frac_transc_st*(MD_ntransc - MD_transc);
%TCRd 
TCRd = TCR / observed;
nTCR = observed-TCR;
%%%%%%%%% What about differences in chromatin?
%%%%%%%%% pull out expressed / notexpressed pairs adjacent to each other to
%%%%%%%%% look for differences.
genepos = readtable('gene_positions.txt','delimiter','\t'); %load in data table
genepos.is_expressed = ismember(genepos.HugoGeneSymbol,mut_by_gene_t.genenames(mut_by_gene_t.is_expressed)); 
genepos.is_notexpressed = ismember(genepos.HugoGeneSymbol,mut_by_gene_t.genenames(mut_by_gene_t.is_notexpressed));
genepos = genepos(genepos.is_expressed | genepos.is_notexpressed,:); %keep genes marked expr or notexpr
genepos = unique(genepos,'rows'); %drop repeated genes
chromatin_pairs = zeros(size(mut_by_gene_t,1),1); % set up output
d = 10^5; %distance in bp that genes must be within
for n = 1:24 % iterate through chromosomes
    switch n    
        case n <= 23
            chrom = num2str(n);
        case n == 24
            chrom = 'x';
    end
    genepos_t = genepos(strcmpi(genepos.ChromosomeName,chrom),:); % pull out genes on chromosome
    nexpr = genepos_t{genepos_t.is_notexpressed,'HugoGeneSymbol'}; %pull out notexpressed gene names
    for i = 1:length(nexpr) %for each not expressed gene
        g = strcmp(genepos_t.HugoGeneSymbol,nexpr{i});
        stp = genepos_t{g,'GeneStart_bp_'}; %start position
        fnp = genepos_t{g,'GeneEnd_bp_'}; %end position
        is_close = genepos_t.is_expressed & (abs(genepos_t.GeneStart_bp_ - stp) <= d ...
            | abs(genepos_t.GeneStart_bp_ - fnp) <= d | abs(genepos_t.GeneEnd_bp_ - stp) <= d ...
            | abs(genepos_t.GeneEnd_bp_ - fnp) <= d);%find any expressed genes within "d" distance
        if any(is_close)
            chromatin_pairs(strcmp(genepos.HugoGeneSymbol,nexpr{i})) = 1; %record nexpr gene
            chromatin_pairs(ismember(genepos.HugoGeneSymbol,genepos_t{is_close,'HugoGeneSymbol'})) = 1; %record expr genes
        end
    end
end
mut_by_gene_t.chromatin_pairs = chromatin_pairs; % add chromatin pairs to mutation list
%find mutation rate of expr or nexpr nearby pairs
MR_expr_transc_CP = sum(mut_by_gene_t.mutations_transcribed(mut_by_gene_t.is_expressed & mut_by_gene_t.chromatin_pairs)) / sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_expressed & mut_by_gene_t.chromatin_pairs));
MR_expr_ntransc_CP = sum(mut_by_gene_t.mutations_nottranscribed(mut_by_gene_t.is_expressed & mut_by_gene_t.chromatin_pairs)) / sum(mut_by_gene_t.start_nottransc(mut_by_gene_t.is_expressed & mut_by_gene_t.chromatin_pairs));
MR_nexpr_transc_CP = sum(mut_by_gene_t.mutations_transcribed(mut_by_gene_t.is_notexpressed & mut_by_gene_t.chromatin_pairs)) / sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_notexpressed & mut_by_gene_t.chromatin_pairs));
MR_nexpr_ntransc_CP = sum(mut_by_gene_t.mutations_nottranscribed(mut_by_gene_t.is_notexpressed & mut_by_gene_t.chromatin_pairs)) / sum(mut_by_gene_t.start_nottransc(mut_by_gene_t.is_notexpressed & mut_by_gene_t.chromatin_pairs));

%compute relative mutation rates
MD_transc_CP = MR_expr_transc_CP / MR_nexpr_transc_CP;
MD_ntransc_CP = MR_expr_ntransc_CP / MR_nexpr_ntransc_CP;
%effect of chromatin is difference in mutation rates
frac_transc_st = sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_expressed & mut_by_gene_t.chromatin_pairs))...
    / (sum(mut_by_gene_t.start_nottransc(mut_by_gene_t.is_expressed & mut_by_gene_t.chromatin_pairs))...
    + sum(mut_by_gene_t.start_transc(mut_by_gene_t.is_expressed & mut_by_gene_t.chromatin_pairs)));
chr = frac_transc_st * (MD_transc_CP - MD_transc) + (1-frac_transc_st) * (MD_ntransc_CP - MD_ntransc);
if chr < 0;
    chr = 1E-5;
end

figure()
fig = gcf;
fig.InvertHardcopy='off';
fig.Units='inches';
fig.Position(3:4) = [4.5 3]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters'; 
%b=bar([MD_transc MD_transc_CP; MD_ntransc MD_ntransc_CP],'EdgeColor','black');
b= bar([frac_transc_st * (MD_transc) + (1-frac_transc_st)*MD_ntransc; frac_transc_st * (MD_transc_CP) + (1-frac_transc_st)*MD_ntransc_CP],.5,'EdgeColor','black'); 
label = {'All genes','Close genes'};
%b(2).FaceColor = [1,0,0];
b(1).FaceColor = [1,0,0];
ylim([0,1]);
ylabel({'Expressed / Not Expressed','mutation rate'}); 
title(tumortype);
ax=gca; 
ax.XTickLabel = label;
ax.FontName=FontName; 
ax.FontSize=FontSize; 
ax.LabelFontSizeMultiplier=1; 
ax.TitleFontSizeMultiplier=1;
ax.TickDir='out'; 
ax.LineWidth=linewidth; 
t = horzcat(st,'>',fn);
t=text(.05,.95,t,'Units','Normalized');
t.FontName=FontName;
t.FontSize = FontSize;
ax.Box='off'; 
ax.Color='none'; 
fig.Color='none'; 
print(horzcat(mutfilename(1:4),'bar_chr_strand.svg'),'-dsvg','-painters');
fig.Color = [1 1 1];

%unexplained is rest of difference
unexplained = 1-observed-TCR;
%now create pie chart for differences
figure()
fig = gcf; 
fig.InvertHardcopy='off';
fig.Units='inches';
fig.Position(3:4) = [3 3]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
P = pie([observed,TCR,unexplained],{'obs','TCR','?'});
P(2).FontSize = FontSize;
P(4).FontSize = FontSize;
P(2).FontName = FontName;
P(4).FontName = FontName;
P(7).LineStyle = '--';
P(6).FontName = FontName;
P(6).FontSize = FontSize;
P(8).FontName = FontName;
P(8).FontSize = FontSize;
colormap([.8,.8,.8;    1.0,0.5,0;    1,.35,0;    1 1 1])
title(tumortype);
ax=gca; 
ax.FontName = FontName;
ax.FontSize = FontSize;
ax.Color='none';
fig.Color='none';
print(horzcat(mutfilename(1:4),'_pie_expr_chr_TCR.svg'),'-dsvg');
fig.Color = [1 1 1]; 



%% What fraction of mutations are missing? 
%First grab all of the missense mutations from expr and nonexpr, 
temp_t = result_table(~isnan(result_table.norm_mis_by_gene) & ~isinf(result_table.norm_mis_by_gene) & ~(~result_table.is_expressed & ~result_table.is_notexpressed),:);
temp_expr_mis = temp_t{temp_t.is_expressed,'norm_mis_by_gene'};
temp_noexpr_mis = temp_t{temp_t.is_notexpressed,'norm_mis_by_gene'};
%make ecdf of expr mut rates 
[q1,m1]=ecdf(temp_expr_mis);
%find nonexpr mut rates of equivilent quantiles
m2=quantile(temp_noexpr_mis,q1);
%convert mut rates of expr to mut rates of equivilent quantiles in nonexpr
expected_t = table(m1,m2);
expected_expr_normmis = zeros(size(temp_expr_mis));
for x=1:length(expected_expr_normmis)
    expected_expr_normmis(x) = min(expected_t{expected_t.m1 == temp_expr_mis(x),'m2'});
end
%turn back into # of missense mutations
expected_expr_mis = (((10.^expected_expr_normmis) * (tumor_number / scale_factor)) .* temp_t{temp_t.is_expressed,'length_AAs'}) -1;
%sum up observed, expected mutationsc
obs_exp = [sum(temp_t.mis_by_gene) sum([temp_t{temp_t.is_notexpressed,'mis_by_gene'}; expected_expr_mis])];
%add to both the observed mutations in unclear expression genes 
obs_exp = obs_exp + sum(result_table{~isnan(result_table.norm_mis_by_gene) & ~isinf(result_table.norm_mis_by_gene) & (~result_table.is_expressed & ~result_table.is_notexpressed),'mis_by_gene'});
lost = (obs_exp(2)-obs_exp(1))/obs_exp(2);
figure()
fig = gcf; 
fig.InvertHardcopy='off';
fig.Units='inches';
fig.Position(3:4) = [1 1]; 
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
quant_results.Avg_muts_lost_to_selection = round((obs_exp(2)-obs_exp(1))/tumor_number);
quant_results.percent_muts_lost_to_selection = (obs_exp(2)/obs_exp(1))*100;
P = pie([1-lost,lost],{num2str(round(obs_exp(1)/tumor_number)),''});
P(2).FontSize = FontSize;
P(2).FontName = FontName;
P(3).LineStyle = '--';
colormap([0.8,0.8,0.8; 1,1,1])
ax=gca; 
ax.Color='none';
fig.Color='none';
print(horzcat(mutfilename(1:4),'_pie_expr.svg'),'-dsvg');
fig.Color = [1 1 1]; 
clear ax fig lost obs_exp expected_expr_mis expected_expr_normmis expected_t m2 m1 q1 temp_noexpr-_mis temp_expr_mis temp_t


%% Can we pull out significantly depleted genes?
%controling for expression
result_table.logexpr =  log2(result_table.expr_by_tumor_geomean); % add expression score;
obs_exp_muts = result_table(result_table.is_expressed,:); %only use expressed genes
%bin by expression: 
bs = 0.25; %bin size;
mean_mutrate = zeros(height(obs_exp_muts),1);
mean_expression = zeros(height(obs_exp_muts),1);
flag = zeros(size(mean_mutrate));
binsize = zeros(size(mean_mutrate));
for b = 1:height(obs_exp_muts); %for each gene
    ge = obs_exp_muts.logexpr(b);
    gb = obs_exp_muts(obs_exp_muts.logexpr >= (ge - bs) & obs_exp_muts.logexpr <= (ge + bs),:); %pull genes close in expr
    mean_mutrate(b) = sum(gb.mis_by_gene) / sum(gb.length_AAs); %compute mutrate
    mean_expression(b) = mean(gb.logexpr); %compute average expression
    binsize(b) = height(gb);
    if height(gb) < 50;
        flag(b) = -1;
    end
end
obs_exp_muts.expected = mean_mutrate .* obs_exp_muts.length_AAs;
obs_exp_muts.Pval = cdf('Poisson',obs_exp_muts.mis_by_gene,obs_exp_muts.expected);
obs_exp_muts.obs_exp = (obs_exp_muts.mis_by_gene+1) ./ (obs_exp_muts.expected+1);
%not capable of being significant depl?
%obs_exp_muts.obs_exp = (obs_exp_muts.mis_by_gene - obs_exp_muts.expected) ./ obs_exp_muts.expected;
obs_exp_muts.P33 = cdf('Poisson',.33.*obs_exp_muts.expected,obs_exp_muts.expected);
obs_exp_muts.never = obs_exp_muts.P33 > 0.01;
obs_exp_muts.sigdep = obs_exp_muts.Pval < 0.001 & obs_exp_muts.obs_exp <.5;
obs_exp_muts.ranged = obs_exp_muts.obs_exp;
obs_exp_muts.ranged(obs_exp_muts.ranged > 2^2) = 2^2;
obs_exp_muts.ranged(obs_exp_muts.ranged < 2^-2) = 2^-2;
obs_exp_muts.invpval = 1-obs_exp_muts.Pval;
writetable(obs_exp_muts(:,{'mut_genes','invpval'}),horzcat(mutfilename(1:4),'_depl.rnk'),'FileType','text','delimiter','\t');
writetable(obs_exp_muts(:,{'mut_genes','sigdep','obs_exp'}),horzcat(mutfilename(1:4),'_obsexp.txt'),'FileType','text','delimiter','\t');

[~,edges] = histcounts(log2(obs_exp_muts.ranged(~obs_exp_muts.never)));
edges =2.^edges;
figure()
fig = gcf;
fig.InvertHardcopy='off';
fig.Units='inches'; 
fig.Position(3:4) = [2.25 1.5];
fig.PaperPositionMode='Auto';
fig.Renderer='painters';
hold on
h=histogram(obs_exp_muts.ranged(~obs_exp_muts.never),edges);
h.EdgeColor=[1 1 1];
h.FaceAlpha = 1;
h2=histogram(obs_exp_muts.ranged(~obs_exp_muts.never & obs_exp_muts.sigdep),edges);
h2.EdgeColor=[1 1 1];
h2.FaceAlpha = 1;
ax=gca;
ax.XScale='log';
ax.XTick = 2.^[-2 -1 0 1 2];
Label = ax.XTick;
Label(ax.XTick<1) = 1./ax.XTick(ax.XTick<1);
ax.XTickLabel=Label;
ax.TitleFontSizeMultiplier=1;
ax.LabelFontSizeMultiplier = 1; 
xlabel('Fold Depl''n | Fold Enrich')
ylabel('# Genes')
title(tumortype);
ax.Color = 'none';
ax.FontName = FontName;
ax.FontSize = FontSize;
ax.TickDir = 'out';
ax.TickLength = [0.04,0.03];
ax.LineWidth=linewidth;
fig.Color = 'none';
print(horzcat('histogram_',tumortype,'_depl_muts.svg'),'-dsvg');
fig.Color = [1 1 1];
