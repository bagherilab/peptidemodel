function PM_TASKS_Shuffle(TASKID, INDEX)

SETTINGS = PM_SETTINGS(); 

switch TASKID     
%% Process data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        d = load(SETTINGS.vectMat);
        filename = [SETTINGS.path SETTINGS.dataMat '_Shuffled'];
        
        % Get response vector.
        response = processDataY(SETTINGS.dataMat, 'STE2', {'EC50'});
        response = response{1};
        
        % Shuffle response vectors (leaving the 101th one unshuffled).
        s.response{51,1} = response;
        n = length(response);
        for i = 1:50
            rng(i);
            s.response{i} = response(randperm(n));
        end
        
        % Create predictor matrices for each vectorization.
        fprintf('|| Processing data by P ');
        s.predictor_P = processDataX(SETTINGS.dataMat, 'STE2', 1, ...
            {'EC50'}, d.v8, SETTINGS.map);
        
        fprintf('\n|| Processing data by R ');
        s.predictor_R = processDataX(SETTINGS.dataMat, 'STE2', 2, ...
            {'EC50'}, d.v13, SETTINGS.map);
        
        fprintf('\n');
        save(filename, '-struct', 's');
        
%% Run PLSR for shuffled data by property (HPC version) %%%%%%%%%%%%%%%%%%%
    case 2
        filename = [SETTINGS.path SETTINGS.dataMat '_Shuffled'];
        v = load(SETTINGS.vectMat);
        vects = v.v8;
        fprintf('|| Regression by P \n');
        
        % Select subset of vectorizations based on job INDEX.
        start = (INDEX - 1)*25 + 1;
        stop = min(INDEX*25, SETTINGS.numVects(1));
        vects = vects(start:stop, :);

        % Run regressions.
        d = load(filename);
        data.predictor = d.predictor_P(start:stop, :);
        data.response = d.response;
        s = runShuffle(data, SETTINGS.comp, 51, vects);

        % Save results.
        filename = [SETTINGS.path SETTINGS.resultsMat '_Shuffled_P_' num2str(INDEX)];
        save(filename, '-struct', 's');

%% Run PLSR for shuffled data by residue (HPC version) %%%%%%%%%%%%%%%%%%%%
    case 3
        filename = [SETTINGS.path SETTINGS.dataMat '_Shuffled'];
        v = load(SETTINGS.vectMat);
        vects = v.v13;
        fprintf('|| Regression by R \n');
        
        % Select subset of vectorizations based on job INDEX.
        start = (INDEX - 1)*25 + 1;
        stop = min(INDEX*25, SETTINGS.numVects(2));
        vects = vects(start:stop, :);

        % Run regressions.
        d = load(filename);
        data.predictor = d.predictor_R(start:stop, :);
        data.response = d.response;
        s = runShuffle(data, SETTINGS.comp, 51, vects);

        % Save results.
        filename = [SETTINGS.path SETTINGS.resultsMat '_Shuffled_R_' num2str(INDEX)];
        save(filename, '-struct', 's');     

%% Merge by property results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4
        s = struct();
        fprintf('|| Merge regression by P \n');
        
        for iInd = 1:11
            filename = [SETTINGS.path SETTINGS.resultsMat '_Shuffled_P_' num2str(iInd) '.mat'];
            d = load(filename);

            for iRep = 1:51
                rep = ['y' num2str(iRep)];
                fields = fieldnames(d.(rep));

                for iField = 1:length(fields)
                    s.(rep).(fields{iField}) = d.(rep).(fields{iField});
                end
            end
        end

        filename = [SETTINGS.path SETTINGS.resultsMat '_Shuffled_P'];
        save(filename, '-struct', 's');
        
%% Merge by residue results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5
        s = struct();
        fprintf('|| Merge regression by R \n');

        for iInd = 1:328
            filename = [SETTINGS.path SETTINGS.resultsMat '_Shuffled_R_' num2str(iInd) '.mat'];
            d = load(filename);

            for iRep = 1:51
                rep = ['y' num2str(iRep)];
                fields = fieldnames(d.(rep));

                for iField = 1:length(fields)
                    s.(rep).(fields{iField}) = d.(rep).(fields{iField});
                end
            end
        end

        filename = [SETTINGS.path SETTINGS.resultsMat '_Shuffled_R'];
        save(filename, '-struct', 's'); 

%% Parse results by property s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 6
        filename = [SETTINGS.path SETTINGS.analysisMat '_Shuffled'];
        fprintf('|| Parsing by property regression \n');
        results_file = [SETTINGS.resultsMat '_Shuffled_P'];
        
        yids = {};
        for i = 1:51
            yids = [yids ['y' num2str(i)]];
        end
        s.P = parseResults(results_file, yids, SETTINGS.vectLengths(1), ...
                SETTINGS.q2_cutoff, SETTINGS.q2_scaled_cutoff);
            
        fprintf('\n');
        save(filename, '-struct', 's');
        
%% Parse results by residue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 7
        filename = [SETTINGS.path SETTINGS.analysisMat '_Shuffled'];
        fprintf('|| Parsing by residue regression \n');
        results_file = [SETTINGS.resultsMat '_Shuffled_R'];
        
        yids = {};
        for i = 1:51
            yids = [yids ['y' num2str(i)]];
        end
        s.R = parseResults(results_file, yids, SETTINGS.vectLengths(2), ...
                SETTINGS.q2_cutoff, SETTINGS.q2_scaled_cutoff);
            
        fprintf('\n');
        save(filename, '-struct', 's', '-append');
        
%% Calculate q2 distributions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 8
        d = load(SETTINGS.analysisMat);
        d2 = load([SETTINGS.analysisMat '_Shuffled']);
        filename = [SETTINGS.path SETTINGS.d3Save];
        bins = -0.5:0.05:0.5;
        
        for iGroup = 1:2
            group = SETTINGS.groups{iGroup};
            data = d.([group '_STE2']).q2(1,:);
            dataAvg = mean(d2.(group).q2);
                
            % Bin Q2 data.
            freqs = histcounts(data, bins);
            sFreqs = histcounts(dataAvg, bins);

            % Normalize by number of vectorizations to get frequency.
            freqs = freqs/SETTINGS.numVects(iGroup);
            sFreqs = sFreqs/SETTINGS.numVects(iGroup);
                
            % Save results to .csv file.
            csvwrite([filename 'SHUFFLE_' group '.csv'], [freqs; sFreqs]);
        end
end