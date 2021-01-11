function PM_TASKS_Validate(TASKID)

SETTINGS = PM_SETTINGS();

switch TASKID
%% Top feature selection with different thresholds %%%%%%%%%%%%%%%%%%%%%%%%
% Features are selected based on weighted frequency above a certain value.
% Various values for the location of this threshold are considered.
    case 1
        fprintf('|| Testing different feature selection thresholds \n');
        d = load(SETTINGS.analysisMat);
        [feats, featsByGroup] = makeFullFeatList();
        s.thresholds = 0.2:0.01:0.8;
        
        for iRecep = 1:SETTINGS.nReceps
            p = d.VIP_P.WM(iRecep,:);
            r = d.VIP_R.WM(iRecep,:);
            [avgSort, avgInds] = sort((r + p)/2, 'descend'); % sort by average
            
            % Try each threshold.
            for iThresh = 1:length(s.thresholds)
                avgThresh = find(avgSort < s.thresholds(iThresh), 1) - 1;
                avgTopInds = avgInds(1:avgThresh);
                
                % Save values.
                s.top_feat_vips{iThresh,iRecep} = avgSort(1:avgThresh);
                s.top_feat_inds{iThresh,iRecep} = avgTopInds;
                s.top_feats{iThresh,iRecep} = feats(avgTopInds);
                s.top_n(iThresh,iRecep) = avgThresh;
                
                % Get grouping characteristics of features.
                for iGroup = 1:2
                    group = SETTINGS.groups{iGroup};
                    histN = 0.5 + SETTINGS.vectLengths(iGroup);
                    groupFeats = featsByGroup.(group)(avgTopInds);
                    s.(group).top_feats{iThresh, iRecep} = groupFeats;
                    s.(group).counts{iRecep}(iThresh, :) = ...
                        histcounts(groupFeats, 0.5:histN);
                    s.(group).percs{iRecep}(iThresh, :) = ...
                        s.(group).counts{iRecep}(iThresh, :)/length(groupFeats);
                end
            end
        end
        
        save(SETTINGS.valMat, '-struct', 's');
    
%% Print out feature selection results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prints out the selected features to a text file.
    case 2 % Print out top feature lists.
        fprintf('|| Print out feature selection results \n');
        d = load(SETTINGS.valMat);
        filename = [SETTINGS.txtSave 'FEATURE_SELECTION.txt'];
        fid = fopen(filename, 'w');
        
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            fprintf(fid, ['RECEPTOR: ' receptor '\n']);
            
            ind = find(d.thresholds == SETTINGS.threshold);
            top_feats = d.top_feats{ind, iRecep};
            top_vips = d.top_feat_vips{ind, iRecep};
            
            for iFeat = 1:length(top_feats)
                fprintf(fid, '  [%.3f] %s\n', top_vips(iFeat), top_feats{iFeat});
            end
        end     

%% Process data for validation PLSR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to the original processing, this will create a structure where
% the only the top features are present in the predictor data set.
    case 3
        fprintf('|| Processing validation data\n');
        d = load(SETTINGS.valMat);
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            top_inds = d.top_feat_inds{d.thresholds == SETTINGS.threshold, iRecep};
            field = ['Data_' receptor];

            [s.(field).predictor, s.(field).Xnorm] = processDataVal(...
                SETTINGS.dataMat, receptor, top_inds, SETTINGS.responses);

            [s.(field).response, s.(field).Ynorm] = processDataY(...
                SETTINGS.dataMat, receptor, SETTINGS.responses);
        end
        save(SETTINGS.valMat, '-struct', 's', '-append');

%% Rerun PLSR using only the top features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each receptor, rerun PLSR using data of only the selected top 
% features against all nine response types.
    case 4
        fprintf('|| Running PLSR on validation data\n');
        d = load(SETTINGS.valMat);
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            temp = runRegress(d.(['Data_' receptor]), SETTINGS.comp, ...
                [1], SETTINGS.responses);
            
            % Parse out results into separate matrix.
            for iResp = 1:length(SETTINGS.responses)
                response = SETTINGS.responses{iResp};
                s.(['Results_' receptor]).(response) = temp.(response).x1;
            end
        end
        save(SETTINGS.valMat, '-struct', 's', '-append');

%% Parse out best q2 from top feature PLSR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse results from top feature only PLSR to get top q2 as well as the 
% top q2 obtained for each response using the grouping exhaustive method.
    case 5
        fprintf('|| Parsing results for best q2\n');
        
        % Get top q2 for each response from feature selected PLSR.
        d = load(SETTINGS.valMat);
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            for iResp = 1:length(SETTINGS.responses)
                response = SETTINGS.responses{iResp};
                q2 = max(d.(['Results_' receptor]).(response).q2);
                s.max_q2(iRecep, iResp) = q2;
            end
        end
        
        % Get top q2 and top response for each response from grouping approach.
        d = load(SETTINGS.analysisMat);
        for iGroup = 1:2
            group = SETTINGS.groups{iGroup};
            for iRecep = 1:SETTINGS.nReceps
                receptor = SETTINGS.receptors{iRecep};
                for iResp = 1:length(SETTINGS.responses)
                    data = d.([group '_' receptor]);
                    q2 = max(data.q2(iResp,:));
                    s.max_group_q2(iGroup, iRecep, iResp) = q2;
                end
                s.top_response(iRecep, iGroup) = data.top_resp;
            end
        end
        
        save(SETTINGS.valMat, '-struct', 's', '-append');

%% Get predictions for all peptides %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each receptor and response, determine the predicted response using
% the top PLSR model against each of the training peptides.
    case 6
        fprintf('|| Prediction responses for all peptides');
        d = load(SETTINGS.dataMat);
        for iRecep = 1:SETTINGS.nReceps
            fprintf('\n | ');
            receptor = SETTINGS.receptors{iRecep};
            for iResp = 1:length(SETTINGS.responses)
                fprintf('.');
                response = SETTINGS.responses{iResp};
                
                M = zeros(1, length(d.seqs));
                for iSeq = 1:length(d.seqs)
                    seq = d.seqs{iSeq};
                    M(iSeq) = getResponse(seq, iRecep, ...
                        iResp, receptor, response, SETTINGS.threshold, ...
                        SETTINGS.dictMat, SETTINGS.valMat);
                end
                
                s.(['Predicted_' receptor]).(response) = M;
            end
        end
        fprintf('\n');
        save(SETTINGS.valMat, '-struct', 's', '-append');

%% Get predictions for novel peptides %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each receptor and response, determine the predicted response using
% the top PLSR model against novel peptides.
    case 7
        fprintf('|| Prediction responses for novel peptides');
        for iRecep = 1:SETTINGS.nReceps
            fprintf('\n | ');
            receptor = SETTINGS.receptors{iRecep};
            for iResp = 1:length(SETTINGS.responses)
                fprintf('.');
                response = SETTINGS.responses{iResp};
                
                for iSeq = 1:length(SETTINGS.novel_seqs)
                    seq = SETTINGS.novel_seqs{iSeq};
                     M(iSeq) = getResponse(seq, iRecep, ...
                        iResp, receptor, response, SETTINGS.threshold, ...
                        SETTINGS.dictMat, SETTINGS.valMat);
                end
                
                s.(['Novel_' receptor]).(response) = M;
            end
        end
        fprintf('\n');
        save(SETTINGS.valMat, '-struct', 's', '-append');
      
%% Print out prediction results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prints out the predicted responses to receptors to a text file.
    case 8 
        fprintf('|| Print out novel prediction results \n');
        d = load(SETTINGS.valMat);
        filename = [SETTINGS.txtSave 'NOVEL_PREDICTIONS.txt'];
        fid = fopen(filename, 'w');
        
        for iSeq = 1:length(SETTINGS.novel_seqs)
                peptide = SETTINGS.novel_names{iSeq};
                fprintf(fid, ['PEPTIDE: ' peptide '\n']);
            for iResp = 1:length(SETTINGS.responses)
                response = SETTINGS.responses{iResp};
                fprintf(fid, ['   RESPONSE: ' response '\n']);
                for iRecep = 1:SETTINGS.nReceps
                    receptor = SETTINGS.receptors{iRecep};
                    pred = d.(['Novel_' receptor]).(response)(iSeq);
                    fprintf(fid, '      %s\t: %.3f\n', receptor, pred);
                end
            end
        end
        
%% Print out promiscuity index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prints out the promiscuity index of the receptors to a text file.
    case 9
        fprintf('|| Print out promiscuity index \n');
        d = load(SETTINGS.dataMat);
        filename = [SETTINGS.txtSave 'PROMISCUITY_INDEX.txt'];
        fid = fopen(filename, 'w');
        
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            pi = calcPromInd(d.(['response_' receptor]).EC50);
            fprintf(fid, '%s\t: %.3f\n', receptor, pi);
        end

%% Print out q2 values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prints out the q2 of the receptors for each grouping and the overall
% feature selected model to a text file.
    case 10
        fprintf('|| Print out q2 values \n');
        d = load(SETTINGS.valMat);
        filename = [SETTINGS.txtSave 'Q2.txt'];
        fid = fopen(filename, 'w');
        fprintf(fid, 'RECEP   RESPONSE  BY TOP   BY P    BY R    P    R\n');
        fprintf(fid, '-----   --------  ------  ------  ------  ---  ---\n');
        
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            for iResp = 1:length(SETTINGS.responses)
                response = SETTINGS.responses{iResp};
                q2t = d.max_q2(iRecep, iResp);
                q2p = d.max_group_q2(1, iRecep, iResp);
                q2r = d.max_group_q2(2, iRecep, iResp);
                
                if iResp == d.top_response(iRecep, 1)
                    p = '*';
                else
                    p = ' ';
                end
                
                if iResp == d.top_response(iRecep, 2)
                    r = '*';
                else
                    r = ' ';
                end
                
                fprintf(fid, '%+5s   %+8s  %+.3f  %+.3f  %+.3f   %s    %s\n', ...
                    receptor, response, q2t, q2p, q2r, p, r);
            end
        end
%% Print out frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 11
        d = load(SETTINGS.analysisMat);
        filename = [SETTINGS.txtSave 'GROUPING_FREQUENCIES.txt'];
        fid = fopen(filename, 'w');
        fprintf(fid, 'group,receptor,ind,val\n');
        
        for iGroup = 1:2
            group = SETTINGS.groups{iGroup};
            for iRecep = 1:SETTINGS.nReceps
                receptor = SETTINGS.receptors{iRecep};
                wfreq = d.vect_wfreq{iGroup, iRecep};
                
                for i = 1:length(wfreq)
                    fprintf(fid, '%s,%s,%d,%f\n', group, receptor, i - 1, wfreq(i));
                end
            end
        end 
%% Print out frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 12
        d = load(SETTINGS.valMat);
        filename = [SETTINGS.txtSave 'FEATURE_FREQUENCIES.txt'];
        fid = fopen(filename, 'w');
        fprintf(fid, 'group,receptor,ind,val\n');
        threshInd = find(d.thresholds == SETTINGS.threshold);
        m = d.top_feat_inds(threshInd,:);
        [~, featsByGroup] = makeFullFeatList();
        
        % Create binary matrix.
        M = {};
        for iGroup = 1:2
            group = SETTINGS.groups{iGroup};
            for iRecep = 1:SETTINGS.nReceps
                receptor = SETTINGS.receptors{iRecep};
                freq = d.(group).percs{iRecep}(threshInd,:);
                
                for i = 1:length(freq)
                    fprintf(fid, '%s,%s,%d,%f\n', group, receptor, i - 1, freq(i));
                end
            end
        end
end

end