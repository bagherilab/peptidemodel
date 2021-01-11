function PM_CONTROLLER(TASKID)
%% This function is the main point for running all data generation,
%% processing, analysis, and plot/figure generation tasks.

PRINT = 1;

switch TASKID
    
%% PART 1: Load and prepare experimental data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following tasks will load and prepare experimental data and place 
% them into matlab structures, which are then saved.
    case 1
        fprintf('RUNNINAG TASK 1: Load and prepare experimental data\n');
        PM_TASKS_Prep(1); % make vectorizations list
        PM_TASKS_Prep(2); % create dictionary mappings
        PM_TASKS_Prep(3); % parse predictor data
        PM_TASKS_Prep(4); % parse response data for each receptor
        PM_TASKS_Prep(5); % process data for by property grouping
        PM_TASKS_Prep(6); % process data for by residue grouping

%% PART 2: Run grouping exhaustive PLSR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following tasks will run grouping exhaustive PLSR for each receptor
% on each of nine responses (three types,three variations). Groupings are
% either By Property or By Residue.
    case 2
        fprintf('RUNNING TASK 2: Run grouping exhaustive PLSR\n');
        PM_TASKS_Run(1); % run by property PLSR
        PM_TASKS_Run(2); % run by residue PLSR

%% PART 3: Parse and analyze results from PLSR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following tasks will parse the results of grouping exhaustive PLSR
% structures into matrices. These results are then analyzed for top
% response type, top performing vectorizations, and top features of the
% top performing vectorizations.
    case 3
        fprintf('RUNNING TASK 3: Parse and analyze results from PLSR\n');
        PM_TASKS_Analyze(1); % process results from by property
        PM_TASKS_Analyze(2); % process results from by residue
        PM_TASKS_Analyze(3); % convert top vectorizations to matrix
        PM_TASKS_Analyze(4); % parse top features based on VIP
        PM_TASKS_Analyze(5); % calculate fold change in EC50
        
%% PART 4: Validate results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following tasks will select the top features based on average
% feature score and then run PLSR on data sets containing only those
% features.
    case 4
        fprintf('RUNNING TASK 4: Validation of results\n');
        PM_TASKS_Validate(1); % try different feature selection thresholds   
        PM_TASKS_Validate(2); % print out feature selection results
        PM_TASKS_Validate(3); % process data for top features only PLSR
        PM_TASKS_Validate(4); % run top features only PLSR
        PM_TASKS_Validate(5); % parse out top q2 
        PM_TASKS_Validate(6); % get predictions for all peptides
        PM_TASKS_Validate(7); % get predictions for novel peptides
        PM_TASKS_Validate(8); % print out novel prediction results
        PM_TASKS_Validate(9); % print out promiscuity index calculations
        PM_TASKS_Validate(10); % print out q2 values
end

end