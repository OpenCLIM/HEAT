function [] = HEAT_step2(inputs2,data,xyz,Dataset,Variable,ExptName)

% Set directory paths
init_HEAT

disp('Running extremes analysis')

% Temporally subset data as required
[data,xyz.dates] = subset_temporal(data,xyz.dates,inputs2.Years,inputs2.AnnSummer);

% Spatially subset data as required


% Calculate percentile if required
if isfield(inputs2,'Pctile')
    % Create empty array for output
    data_pctiles = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs2.Pctile));
    
    for p = 1:length(inputs2.Pctile)
        data_pctiles(:,:,p) = prctile(data,inputs2.Pctile(p),3);
    end
end

% Calculate extreme mean if required
if isfield(inputs2,'ExtremeMeanPctile')
    % Create empty array for output
    data_ExMean = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs2.ExtremeMeanPctile));
    
    for p = 1:length(inputs2.ExtremeMeanPctile)
        
        Txx_temp = data >= prctile(data,inputs2.ExtremeMeanPctile(p),3);
        Txx = nan(size(Txx_temp));
        Txx(Txx_temp == 1) = 1;
        data_Txx = data .* Txx;
        data_ExMean(:,:,p) = squeeze(nanmean(data_Txx,3));
    end
end


% Generate requested output
disp('Producing output')
for o = 1:length(inputs2.OutputType)
    OutputType = string(inputs2.OutputType(o));
    
    % Generate map if requested
    if strcmp(OutputType,'map')
        figure
        UK_subplot(data_ExMean,[Dataset,' ',Variable, ' extreme mean >95th percentile'])
        %                             figure
        %                             UK_subplot(data_pctiles,[Dataset,' ',Variable, ' 95th percentile'])
    end
    
    
    % Generate ARCADIA input csv if requested
    if strcmp(OutputType,'ARCADIA')
        
        % If interested in every single grid cell
        if strcmp(inputs2.Region,'all')
            
            % csv needs saved for every grid box: go through each lat-long
            for i = 1:length(data(:,1,1))
                for j = 1:length(data(1,:,1))
                    
                    % Set csv output file name
                    csv_name = [Outputdir,'/',ExptName,'/',num2str(i),num2str(j),'_',ExptName,'_',Variable,'.csv'];
                    
                    % If csv has been created already
                    if exist(csv_name,'file')
                        % Load existing csv
                        data_csv = csvread(csv_name);
                        % Extract individual grid cell for current dataset
                        data_csv_temp = squeeze(data(i,j,:))';
                        % Concatenate this with the existing csv and resave
                        data_csv = cat(1,data_csv,data_csv_temp);
                        csvwrite(csv_name,data_csv)
                    else
                        % Otherwise extract individual grid cell for current
                        % dataset and create new csv
                        data_csv = squeeze(data(i,j,:))';
                        csvwrite(csv_name,data_csv)
                    end
                end
            end
            
        else % Otherwise take a regional mean
            regmean = calc_reg_mean(data,inputs2.Region);
            
            % Set csv output file name
                    csv_name = [Outputdir,'/',ExptName,'/',inputs2.Region,'_',ExptName,'_',Variable,'.csv'];
                    
                    % If csv has been created already
                    if exist(csv_name,'file')
                        % Load existing csv
                        data_csv = csvread(csv_name);

                        % Concatenate this with the existing csv and resave
                        data_csv = cat(1,data_csv,regmean');
                        csvwrite(csv_name,data_csv)
                    else
                        % Otherwise create new csv
                        csvwrite(csv_name,regmean')
                    end
        end
        
    end
    
end
disp('-----')


