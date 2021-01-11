function PM_TASKS_Run_HPC(TASKID, INDEX)

%% HPC version of _Run_ tasks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tasks 1 and 2 are the parallelized versions of _Run_ tasks 1 and 2. The %
% vectorizations that are used are split into chunks of 25 indexed by the %
% job array number.                                                       %
%                                                                         %
% Tasks 3 and 4 are helper tasks to re-merge the chunks output by tasks 1 %
% and 2 into a single .mat file containing structures for each receptor.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SETTINGS.path = '/home/jsy331/Matlab/';
SETTINGS.vectMat = 'Results/PM_Vectorizations';
SETTINGS.dataMat = 'Results/PM_Data';
SETTINGS.resultsMat  = 'Results/PM_Results';
SETTINGS.receptors = {'STE2', 'HF10', 'PROM6', 'PROM7', 'MUT1', 'PROM3', 'TBBI2'};
SETTINGS.nReceps = length(SETTINGS.receptors);
SETTINGS.responses = {'EC50', 'EC50_log', 'EC50_inv', 'EC10', 'EC10_log', ...
    'EC10_inv', 'UM10', 'UM10_log', 'UM10_inv'};
SETTINGS.numVects = [(2^8 - 1) (2^13 - 1)];
SETTINGS.comp = 5;

switch TASKID
    case 1
        v = load(SETTINGS.vectMat);
        vects = v.v8;
        
        % Select subset of vectorizations based on job INDEX.
        start = (INDEX - 1)*25 + 1;
        stop = min(INDEX*25, SETTINGS.numVects(1));
        vects = vects(start:stop, :);
        
        fprintf('|| Processing response data by P \n');
        
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            fprintf(' | Receptor %s \n', receptor);
            
            % Run regressions.
            filename = [SETTINGS.path SETTINGS.dataMat '_P_' receptor];
            d = load(filename);
            d = d.(['P_' receptor]);
            d.predictor = d.predictor(start:stop, :);
            s = runRegress(d, SETTINGS.comp, vects, SETTINGS.responses);
            clear d
            
            % Save results.
            filename = [SETTINGS.path SETTINGS.resultsMat '_P_' receptor '_' num2str(INDEX)];
            save(filename, '-struct', 's');
        end
    case 2
        v = load(SETTINGS.vectMat);
        vects = v.v13;
        
        % Select subset of vectorizations based on job INDEX.
        start = (INDEX - 1)*25 + 1;
        stop = min(INDEX*25, SETTINGS.numVects(2));
        vects = vects(start:stop, :);
        
        fprintf('|| Processing response data by R \n');
        
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            fprintf(' | Receptor %s \n', receptor);
            
            % Run regressions.
            filename = [SETTINGS.path SETTINGS.dataMat '_R_' receptor];
            d = load(filename);
            d = d.(['R_' receptor]);
            d.predictor = d.predictor(start:stop, :);
            s = runRegress(d, SETTINGS.comp, vects, SETTINGS.responses);
            clear d
            
            % Save results.
            filename = [SETTINGS.path SETTINGS.resultsMat '_R_' receptor '_' num2str(INDEX)];
            save(filename, '-struct', 's');
        end   
    case 3
        receptor = SETTINGS.receptors{INDEX};
        fprintf(' | Receptor %s \n', receptor);
        s = struct();

        for iInd = 1:11
            filename = [SETTINGS.path SETTINGS.resultsMat '_P_' receptor '_' num2str(iInd) '.mat'];
            d = load(filename);

            for iResp = 1:length(SETTINGS.responses)
                response = SETTINGS.responses{iResp};
                fields = fieldnames(d.(response));

                for iField = 1:length(fields)
                    s.(response).(fields{iField}) = d.(response).(fields{iField});
                end
            end
        end

        filename = [SETTINGS.path SETTINGS.resultsMat '_P_' receptor];
        save(filename, '-struct', 's');
    case 4
        receptor = SETTINGS.receptors{INDEX};
        fprintf(' | Receptor %s \n', receptor);
        s = struct();

        for iInd = 1:328
            filename = [SETTINGS.path SETTINGS.resultsMat '_R_' receptor '_' num2str(iInd) '.mat'];
            d = load(filename);

            for iResp = 1:length(SETTINGS.responses)
                response = SETTINGS.responses{iResp};
                fields = fieldnames(d.(response));

                for iField = 1:length(fields)
                    s.(response).(fields{iField}) = d.(response).(fields{iField});
                end
            end
        end

        filename = [SETTINGS.path SETTINGS.resultsMat '_R_' receptor];
        save(filename, '-struct', 's');
end 
    
end