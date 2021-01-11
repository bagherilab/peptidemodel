function SETTINGS = PM_SETTINGS()

%% File locations.
SETTINGS.path        = '/homes/jessica_yu/PRJ_PeptideModel/';
SETTINGS.dataFile    = 'Data/PM_DataSets.xlsx';
SETTINGS.dictFile    = 'Data/PM_Dictionaries.xlsx';
SETTINGS.vectMat     = 'Results/PM_Vectorizations';
SETTINGS.dictMat     = 'Results/PM_Dictionaries';
SETTINGS.dataMat     = 'Results/PM_Data';
SETTINGS.resultsMat  = 'Results/PM_Results';
SETTINGS.analysisMat = 'Results/PM_Analysis';
SETTINGS.valMat      = 'Results/PM_Validation';
SETTINGS.plotSave    = 'Plots/PM_PLOT_';
SETTINGS.figSave     = 'Figures/PM_FIG_';
SETTINGS.txtSave     = 'Results/PM_';
SETTINGS.d3Save      = 'Results/D3/PM_';

%% Data properties.
SETTINGS.nFeats = 13*8;
SETTINGS.vectLengths = [8 13];
SETTINGS.numVects = [(2^8 - 1) (2^13 - 1)];
SETTINGS.receptors = {'STE2', 'HF10', 'PROM6', 'PROM7', 'MUT1', 'PROM3', 'TBBI2'};
SETTINGS.nReceps = length(SETTINGS.receptors);
SETTINGS.responses = {'EC50', 'EC50_log', 'EC50_inv', 'EC10', 'EC10_log', ...
    'EC10_inv', 'UM10', 'UM10_log', 'UM10_inv'};

%% Regression settings.
SETTINGS.comp = 5;
SETTINGS.map = reshape(1:104,13,[])';
SETTINGS.groups = {'P', 'R'};

%% Analysis settings.
SETTINGS.q2_cutoff = 0.25;
SETTINGS.q2_scaled_cutoff = 0.75;

%% Validation settings.
SETTINGS.threshold = 0.5;

%% Plot labeling.
SETTINGS.response_names = {'EC50','log(EC50)','inv(EC50)','EC10', ...
    'log(EC10)','inv(EC10)','UM10','log(UM10)','inv(UM10)'};
SETTINGS.group_names = {'Property', 'Residue'};
SETTINGS.group_labels = {{'P1','P2','P3','P4','P5','P6','P7','P8'};
    {'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13'}};
SETTINGS.scan = {'A','Y','D','K'};
SETTINGS.alpha = 'WHWLQLKPGQPMY';
    
for iGroup = 1:2
    for iRecep = 1:SETTINGS.nReceps
        SETTINGS.legend{(iGroup - 1)*SETTINGS.nReceps + iRecep} = ...
            [SETTINGS.groups{iGroup} ' - ' SETTINGS.receptors{iRecep}];
    end
end

%% Plot styling.
cG = [0.4660    0.6740    0.1880];
cR = [0.6350    0.0780    0.1840];
cB = [0         0.4470    0.7410];
cP = [ 71  46 138]/255;
cO = [255 156  68]/255;
c1 = [ 33 120  99]/255;
c2 = [240 136 136]/255;
c3 = [  0 163 189]/255;
c4 = [238 197 78]/255;

SETTINGS.recep_colors = [cG; cR; cB; c1; c2; c3; c4];
SETTINGS.group_styles = {'-', ':'};
SETTINGS.group_colors = [cP; cO];

%% Novel data.
SETTINGS.novel_names = {'alpha','A','B','C','D','E','F','G','H'};
SETTINGS.novel_seqs = {'WHWLQLKPGQPMY'
	'LHLLALKPGQPMY'
	'WHWLQLKPGEPLYGR'
	'LHLLALKPGQPLYGR'
	'LHLLAGQPGESLYGR'
	'-HALALKPGEPMY'
	'--ALALKPGEPMY'    
    'YHADQLKPGQPKY'
    'VHDLQLDPGQPLY'};
            
end