function [] = HEAT_step3(inputs,DataType,Variable)
% HEAT_step3(inputs,DataType,Variable)
%
% Run the third step of HEAT to output workflow data. This is data that is
% used by other OpenCLIM models in workflows on DAFNI. Currently, the only
% workflow that HEAT feeds into is the ARCADIA heat-mortality model. Other
% models can be added as required.
%
% It uses the same loading process as HEAT_step2, whereby each dataset
% (e.g. simulation) is loaded in turn either from raw or derived data for
% the correct spatial and temporal domain. There is an option to output
% regional means for the ARCADIA model ? if a region is specified then
% whole spatial domain is loaded before taking the regional mean.
%
% ARCADIA currently only requires daily mean temperature or sWBGT as far
% as I know. Daily max. and daily min. temperatures can be exported also.
%
% TO DO: Check if ARCADIA always requires a full annual cycle of data. Katie's
% document on Teams suggests so, in which case the AnnSummer subsetting
% should be removed and the calendar for UKCP18 models needs converted from
% 360 to 365 days (I can use the ISIMIP2b method for this).

%% Generate ARCADIA-ready csv if requested
disp(' ')
disp('Producing workflow output for ARCADIA')
disp('-----')

            
            %% Output ARCADIA data for grid cell or regional mean
            % If interested in every single grid cell
            if ~isfield(inputs,'OutputRegion')
                
                
                % csv needs saved for every grid box: go through each lat-long
                for i = 1:length(data(:,1,1))
                    for j = 1:length(data(1,:,1))
                        
                        % Only save if the point is on land
                        if LSM(i,j) == 1
                            
                            % Set csv output file name
                            if strcmp(DataType,'HadUKGrid')
                                csv_name = [Outputdir,'/',inputs.ExptName,'/',num2str(i),'_',num2str(j),'_HadUK-Grid-',Dataset,'_',Variable,'.csv'];
                            elseif strcmp(DataType,'UKCP18')
                                csv_name = [Outputdir,'/',inputs.ExptName,'/',num2str(i),'_',num2str(j),'_UKCP18-',res(1:end-1),'_',Variable,'.csv'];
                            end
                            
                            
                            %                         % If csv has been created already
                            %                         if exist(csv_name,'file')
                            %                             % Load existing csv
                            %                             data_csv = csvread(csv_name);
                            %                             % Extract individual grid cell for current dataset
                            %                             data_csv_temp = squeeze(data(i,j,:))';
                            %                             % Concatenate this with the existing csv and resave
                            %                             data_csv = cat(1,data_csv,data_csv_temp);
                            %                             csvwrite(csv_name,data_csv)
                            %                         else
                            %                             % Otherwise extract individual grid cell for current
                            %                             % dataset and create new csv
                            %                             data_csv = squeeze(data(i,j,:))';
                            %                             csvwrite(csv_name,data_csv)
                            %                         end
                            
                            % Save the data as a .csv file
                            dlmwrite(csv_name,squeeze(data(i,j,:))','-append','newline','pc','delimiter',',','precision',4);
                        end
                    end
                end
                
                
            else % Otherwise take a regional mean
                regmean = calc_reg_mean(data,inputs.Region);
                
                
                % Set csv output file name
                csv_name = [Outputdir,'/',inputs.ExptName,'/',inputs.Region,'_',inputs.ExptName,'_',Variable,'.csv'];
                
%                 % If csv has been created already
%                 if exist(csv_name,'file')
%                     % Load existing csv
%                     data_csv = csvread(csv_name);
%                     
%                     % Concatenate this with the existing csv and resave
%                     data_csv = cat(1,data_csv,regmean');
%                     csvwrite(csv_name,data_csv)
%                 else
%                     % Otherwise create new csv
%                     csvwrite(csv_name,regmean')
%                 end
                
                % Save the data as a .csv file
                dlmwrite(csv_name,regmean','-append','newline','pc','delimiter',',','precision',4);
                
            end
        end
        
    end
end
disp('Step 3 complete. These files have been created:')
ls(['/data/outputs/',inputs.ExptName])
end

