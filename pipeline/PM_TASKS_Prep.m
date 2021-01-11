function PM_TASKS_Prep(TASKID)

SETTINGS = PM_SETTINGS(); 

switch TASKID
%% Create binary vectorizations matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each vectorization for length N can be written as a binary string of
% length N, which can be quickly computed from the decimal representation
% to the binary representation.
    case 1
        filename = [SETTINGS.path SETTINGS.vectMat];
        fprintf('|| Creating vectorizations\n');
        for i = 1:length(SETTINGS.vectLengths)
            N = SETTINGS.vectLengths(i);
            s.(['v' num2str(N)]) = createVects(N);
        end
        save(filename, '-struct' , 's');
        
%% Convert dictionaries mapping amino acids to values %%%%%%%%%%%%%%%%%%%%%
% VHSE vectorization consists of eight properties for each of the 20 amino
% acids. Dictionaries convert residue key to the VHSE value.
    case 2
        filename = [SETTINGS.path SETTINGS.dictMat];
        fprintf('|| Processing dictionaries\n');
        s = processDicts(SETTINGS.dictFile);
        save(filename, '-struct' , 's');
        
%% Parse predictor data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each sequence of data is parsed from spreadsheet. The sequences must be 
% in the second column of the sheet and the names of the sequences must be
% in the first columne of the sheet. The ordering of sequences must be the
% same in all sheets; only the first sheet is used to parse predictors.
    case 3
        filename = [SETTINGS.path SETTINGS.dataMat];
        fprintf('|| Parsing predictor data\n');
        s = processPredictors(SETTINGS.dataFile, 8, SETTINGS.dictMat);
        save(filename, '-struct' , 's');
    
%% Parse response data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Three different response types are parsed from spreadsheet: response at
% 10 uM (UM10), half maximal response (EC50), and 10% maximal response
% (EC10). Binary response is also parsed for use in filtering data.
    case 4
        filename = [SETTINGS.path SETTINGS.dataMat];
        fprintf('|| Parsing response data ');
        for i = 1:SETTINGS.nReceps
            fprintf('.');
            receptor = SETTINGS.receptors{i};
            s = processResponses(SETTINGS.dataFile, receptor);
            save(filename, '-struct', 's', '-append');
        end
        fprintf('\n');
        
%% Process data for by property grouping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gets every combination of predictor and vectorization for the by property
% grouping and saves to a cell array. Data is appropriately sliced and
% normalized.
    case 5
        d = load(SETTINGS.vectMat);
        vects = d.v8;
        fprintf('|| Processing response data by P ');
        for iRecep = 1:SETTINGS.nReceps
            fprintf('\n | ');
            receptor = SETTINGS.receptors{iRecep};
            filename = [SETTINGS.path SETTINGS.dataMat '_P_' receptor];
            field = ['P_' receptor];
            s.(field).predictor = processDataX(SETTINGS.dataMat, ...
                receptor, 1, SETTINGS.responses, vects, SETTINGS.map);
            s.(field).response = processDataY(SETTINGS.dataMat, ...
                receptor, SETTINGS.responses);
            save(filename, '-struct', 's');
            clear s
        end
        fprintf('\n');
            
%% Process data for by residue grouping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gets every combination of predictor and vectorization for the by residue
% grouping and saves to a cell array. Data is appropriately sliced and
% normalized.
    case 6
        d = load(SETTINGS.vectMat);
        vects = d.v13;
        fprintf('|| Processing response data by R ');
        for iRecep = 1:SETTINGS.nReceps
            fprintf('\n | ');
            receptor = SETTINGS.receptors{iRecep};
            filename = [SETTINGS.path SETTINGS.dataMat '_R_' receptor];
            field = ['R_' receptor];
            s.(field).predictor = processDataX(SETTINGS.dataMat, ...
                receptor, 2, SETTINGS.responses, vects, SETTINGS.map);
            s.(field).response = processDataY(SETTINGS.dataMat, ...
                receptor, SETTINGS.responses);
            save(filename, '-struct', 's');
            clear s
        end
       fprintf('\n');
end

end