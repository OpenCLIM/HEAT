% generate_launch_file.m
% Using a series of GUI inputs, create a file containing all of the
% necessary data to launch HEAT.

%% Select whether to run test
setup_list = {'Setup test','Simulation'};
ind_setup = listdlg('PromptString','0. Run setup test or simulation?','ListString',setup_list,'SelectionMode','single','ListSize',[240 60],'Name','HEAT basic setup guide');
setup = setup_list(ind_setup);

if ind_setup == 1
    % Set some basic experiment details
    inputs.ExptName = 'setup_test';
    inputs.OverwriteExpt = 1;
    inputs.Domain = 'UK';
    inputs.SaveDerivedOutput = 1;
    inputs.OverwriteDerivedOutput = 1;
    
    % Select dataset(s) to use
    inputs.DataType = {'UKCP18','HadUKGrid'};
    inputs.Variable = {'sWBGT'};
    inputs.Dataset = {'GCM-01','GCM-04','60km'};
    inputs.MMM = 1;
    
    % Subsetting of data for analysis
    inputs.TemporalRange = [19900101, 20191230];
    inputs.AnnSummer = 'JJA';
    inputs.SpatialRange = [51,55.5;-11, -5];
    % inputs.Region = []; % 'Scotland','North East','North West','Yorkshire and
    %     % the Humber','East Midlands','West Midlands','East of England','Greater
    %     % London','South East','South West','Wales','Northern Ireland','Isle of
    %     % Man','RoI','Ireland','GB'
    
    % Output types
    inputs.ExtremeMeanPctile = [95 99];
    inputs.PlotAll = 0;
    inputs.OutputType = {'map'};
    inputs.SaveFigures = 1;
    inputs.WorkflowOutput = 'ARCADIA';
    clear setup setup_list ind_setup
    
    disp('Created input file to run a setup test of HEAT.')
    disp('Run setup test by typing:')
    disp('>> HEAT(inputs)')
    disp(' ')
    disp('Setup test should take ~ 1 minute and produce:')
    disp('1. 22 derived netCDF files of sWBGT for UKCP18 GCMs 01 and 04 in the derived data folder')
    disp('2. Two map figures of sWBGT data over Ireland in the output folder, setup_test')
    disp('3. 144 ARCADIA-ready .csv files in the output folder, setup_test')
    disp(' ')
    
    
    return
end


%% Select global or UK
dom_list = {'Global','UK'};
ind_dom = listdlg('PromptString','1. Load global or UK data?','ListString',dom_list,'SelectionMode','single','ListSize',[240 60],'Name','HEAT basic setup guide');
Domain = dom_list(ind_dom);

%% Select whether to use observations, reanalysis or model data
if ind_dom == 1
    data_list = {'Model data (UKCP18 or CMIP5)','Observations (Berkeley Earth temperature data)','Reanalysis (ERA5)'};
    tf_data = zeros(1,3);

else
    data_list = {'Model data (UKCP18 or CMIP5)','Observations (HadUK-Grid)','Reanalysis (ERA5)'};
    tf_data = zeros(1,3);

end
ind_data = listdlg('PromptString','2. Select data product to use.','ListString',data_list,'ListSize',[240 60],'Name','HEAT basic setup guide');

for i = 1:length(ind_data)
    tf_data(ind_data(i)) = 1;
end

%% Select variable
% Global observations of humidity are unavailable, therefore heat stress
% cannot be calculated
if ind_dom == 1 && ind_data == 2
    var_list = {'Mean temperature','Max. temperature','Min. temperature'};
else
    var_list = {'Mean temperature','Max. temperature','Min. temperature','Vapour pressure','Heat stress variable'};
end
ind_var = listdlg('PromptString','3. Select variable to load.','ListString',var_list,'ListSize',[240 80],'Name','HEAT basic setup guide');


%% Select model family/resolution to use

if ind_dom == 1
    
    % If using model data
    if tf_data(1) == 1
        res_list = {'60km'};
    end
    
    if tf_data(2) == 1
        res_list = cat(2,res_list,'TBC (obs. raw)');
    end
    
    if tf_data(3) == 1
        res_list = cat(2,res_list,'0.25� (reanal. raw)');
    end
    
    
else
    % If using model data
    if tf_data(1) == 1
        res_list = {'60km','12km','2.2km'};
    end
    
    if tf_data(2) == 1
        res_list = cat(2,res_list,'1km (obs. only)');
    end
    
    if tf_data(3) == 1
        res_list = cat(2,res_list,'0.25� (reanal. raw)');
    end
end

ind_res = listdlg('PromptString','3. Select resolution to load.','ListString',res_list,'ListSize',[240 80],'Name','HEAT basic setup guide');

tf_res = zeros(1,length(res_list));

for i = 1:length(ind_res)
    tf_res(ind_res(i)) = 1;
end


%% Select individual model simulation
if tf_data(1) == 1
    
    sim_list1 = {};
    sim_list2 = {};
    sim_list3 = {};
    
    if tf_res(1) == 1
        sim_list1 = {'GCM-01','GCM-02','GCM-03','GCM-04','GCM-05','GCM-06','GCM-07','GCM-08','GCM-09','GCM-10','GCM-11','GCM-12','GCM-13','GCM-14','GCM-15',...
            'CMIP5-16','CMIP5-19','CMIP5-21','CMIP5-23','CMIP5-24','CMIP5-25','CMIP5-27','CMIP5-28'};
    end
    
    if tf_res(2) == 1
        sim_list2 = {'RCM-01','RCM-04','RCM-05','RCM-06','RCM-07','RCM-08','RCM-09','RCM-10','RCM-11','RCM-12','RCM-13','RCM-15'};
    end
    
    if tf_res(3) == 1
        sim_list3 = {'CPM-01','CPM-04','CPM-05','CPM-06','CPM-07','CPM-08','CPM-09','CPM-10','CPM-11','CPM-12','CPM-13','CPM-15'};
    end
    
    sim_list = cat(2,sim_list1,sim_list2,sim_list3);
    
    ind_sim = listdlg('PromptString','4a. Select simulation to load.','ListString',sim_list,'ListSize',[240 240],'Name','HEAT basic setup guide');

        
end


%% Write output


