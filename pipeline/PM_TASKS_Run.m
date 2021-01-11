function PM_TASKS_Run(TASKID)

SETTINGS = PM_SETTINGS(); 

switch TASKID
%% Run PLSR for each receptor by property %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLSR is run for each by property groupings on all receptors indicated
% in SETTINGS for all response types.
    case 1
        v = load(SETTINGS.vectMat);
        vects = v.v8;
        fprintf('|| Processing response data by P \n');
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            fprintf(' | Receptor %s \n', receptor);
            
            % Run regressions.
            filename = [SETTINGS.path SETTINGS.dataMat '_P_' receptor];
            d = load(filename);
            s = runRegress(d.(['P_' receptor]), SETTINGS.comp, ...
                vects, SETTINGS.responses);
            clear d
            
            % Save results.
            filename = [SETTINGS.path SETTINGS.resultsMat '_P_' receptor];
            save(filename, '-struct', 's');
        end

%% Run PLSR for each receptor by residue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLSR is run for each by residue grouping on all receptors indicated
% in SETTINGS for all response types.
    case 2
        v = load(SETTINGS.vectMat);
        vects = v.v13;
        fprintf('|| Processing response data by R \n');
        
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            fprintf(' | Receptor %s \n', receptor);
            
            % Run regressions.
            filename = [SETTINGS.path SETTINGS.dataMat '_R_' receptor];
            d = load(filename);
            s = runRegress(d.(['R_' receptor]), SETTINGS.comp, ...
                vects, SETTINGS.responses);
            clear d
            
            % Save results.
            filename = [SETTINGS.path SETTINGS.resultsMat '_R_' receptor];
            save(filename, '-struct', 's');
        end        
end

end