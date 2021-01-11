function PM_TASKS_Analyze(TASKID)

SETTINGS = PM_SETTINGS(); 

switch TASKID
%% Parse regression results for by property grouping %%%%%%%%%%%%%%%%%%%%%%
% Parses the results structures from grouping exhaustive regression by 
% property. Q2 for each receptor/vectorization combination are tracked. 
% Then, top response is calculated and used to filter results for top q2 
% and vectorizations.
    case 1
        filename = [SETTINGS.path SETTINGS.analysisMat];
        fprintf('|| Parsing by property regression \n');
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            fprintf(' | Receptor %s ', receptor);
            results_file = [SETTINGS.resultsMat '_P_' receptor];
            s.(['P_' receptor]) = parseResults(results_file, ...
                SETTINGS.responses, SETTINGS.vectLengths(1), ...
                SETTINGS.q2_cutoff, SETTINGS.q2_scaled_cutoff);
            fprintf('\n');
        end
        save(filename, '-struct', 's');

%% Parse regression results for by residue grouping %%%%%%%%%%%%%%%%%%%%%%%
% Parses the results structures from grouping exhaustive regression by 
% residue. Q2 for each receptor/vectorization combination are tracked. 
% Then, top response is calculated and used to filter results for top q2 
% and vectorizations.
    case 2
        filename = [SETTINGS.path SETTINGS.analysisMat];
        fprintf('|| Parsing by residue regression \n');
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            fprintf(' | Receptor %s ', receptor);
            results_file = [SETTINGS.resultsMat '_R_' receptor];
            s.(['R_' receptor]) = parseResults(results_file, ...
                SETTINGS.responses, SETTINGS.vectLengths(2), ...
                SETTINGS.q2_cutoff, SETTINGS.q2_scaled_cutoff);
            fprintf('\n');
        end
        save(filename, '-struct', 's', '-append');

%% Parse top vectorizations into matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go through each element of each of the top vectorization for every
% receptor/grouping pair. Consolidates results into a binary matrix and a
% weighted matrix in which the binary values are weighted by scaled q2.
% Also calculated frequencies by summing down columns.
    case 3
        fprintf('|| Parsing top vectorizations \n');
        d = load(SETTINGS.analysisMat);
        for iGroup = 1:2
            group = SETTINGS.groups{iGroup};
            for iRecep = 1:SETTINGS.nReceps
                receptor = SETTINGS.receptors{iRecep};
                field = [group '_' receptor];
                data = d.(field);
                
                N = data.q2_thresh;
                v = data.top_vects_sorted(1:N);
                q2 = data.top_q2_scaled(1:N);
                m = []; wm = [];
                
                % Loop through each of the top vectorizations.
                for i = 1:length(v)
                    vec = sprintf(['%0' num2str(SETTINGS.vectLengths(iGroup)) 's'], ...
                        dec2bin(v(i)));
                    for j = 1:length(vec)
                        m(i,j) = str2num(vec(j));
                    end
                    wm(i,:) = m(i,:) * q2(i);
                end
                
                s.vect_mat{iGroup, iRecep} = m;
                s.vect_wmat{iGroup, iRecep} = wm; 
                s.vect_freq{iGroup, iRecep} = sum(m)/size(m,1);
                s.vect_wfreq{iGroup, iRecep} = sum(wm)/size(wm, 1);
            end
        end
        
        save(SETTINGS.analysisMat, '-struct', 's', '-append');

%% Parse out top vip features for top models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through each grouping and receptor and pull out the feature lists of
% the top vectorizations scored by VIP. Then sort the VIP lists and 
% identify the elbow in order to select out the most important features. 
% The occurance of each feature is then averaged across the vectorizations.
    case 4
        fprintf('|| Parsing top features \n');
        d = load(SETTINGS.analysisMat);
        for iGroup = 1:2
            M = zeros(SETTINGS.nReceps, SETTINGS.nFeats);
            WM = zeros(SETTINGS.nReceps, SETTINGS.nFeats);
            
            group = SETTINGS.groups{iGroup};
            for iRecep = 1:SETTINGS.nReceps
                receptor = SETTINGS.receptors{iRecep};
                field = [group '_' receptor];
                N = d.(field).q2_thresh;
                
                % Loop through each of the top N models for given receptor
                % and approach combination
                for i = 1:N
                    vipList = d.(field).VIPs{i};
                    [featInds, featVips] = getElbow(vipList);
                    
                    fullInds = getFullInds(d.vect_mat{iGroup, ...
                        iRecep}(i,:), iGroup, SETTINGS.map);
                    topInds = fullInds(featInds);

                    % Update matrix containing counts.
                    M(iRecep,topInds) = M(iRecep,topInds) + 1;
                    WM(iRecep,topInds) = WM(iRecep,topInds) + featVips'/featVips(1);
                end
                
                M(iRecep, :) = M(iRecep, :)/N;
                WM(iRecep, :) = WM(iRecep, :)/N;
            end
            
            s.(['VIP_' SETTINGS.groups{iGroup}]).M = M;
            s.(['VIP_' SETTINGS.groups{iGroup}]).WM = WM;
        end
        save(SETTINGS.analysisMat, '-struct', 's', '-append');

%% Calculate relative EC50 against response to alpha %%%%%%%%%%%%%%%%%%%%%%
% Loop through each of four scan residues and positions along the peptide
% to calculate the log fold change in EC50 from that against alpha factor.
    case 5
        d = load(SETTINGS.dataMat);
        for iRecep = 1:SETTINGS.nReceps
            receptor = SETTINGS.receptors{iRecep};
            field = ['response_' receptor];
            E0 = d.(field).EC50(end);
            for iScan = 1:length(SETTINGS.scan)
                for iRes = 1:13
                    seq = SETTINGS.alpha; % copy native sequence
                    seq(iRes) = SETTINGS.scan{iScan}; % replace position with scan residue
                    if strcmp(seq, SETTINGS.alpha) % check if resulting sequences is wt
                        s.foldChangeWT{iRecep}(iScan,iRes) = 1;
                    else
                        s.foldChangeWT{iRecep}(iScan,iRes) = 0;
                    end
                    ind = find(ismember(d.seqs, seq)); % find index of sequence
                    s.foldChange{iRecep}(iScan,iRes) = ...
                        d.(field).EC50(ind)/E0; % add relative EC50
                end
            end
        end
        save(SETTINGS.analysisMat, '-struct', 's', '-append');
end

end

function [E, vips] = getElbow(vip_list)
[sorted, inds] = sort(vip_list, 'descend'); % sort by decreasing vip
n = length(sorted);
x = sorted.^2; % x distance
y = (1:n)/n*max(sorted);
y = y.^2'; % y distance
z = sqrt(x + y);
[~, ind] = min(z);
E = inds(1:ind);
vips = sorted(1:ind);
end