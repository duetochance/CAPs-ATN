%% ARGS
% 
%     - handles : dataset
%     - num_seeds: number of seeds, can be 1 or 2 or 3 (max)


%% setSeed 
% v 0.2

function [Seed] = setSeed(handles,num_seeds,seed_name)
if num_seeds == 1
    fprintf("\n ############\n->> Single seeds version \n##################\n");
    %disp(handles.brain_info{1}.dim);
    
    %disp(handles.brain_info{1}.mat);
    % Multiselection is on, so that as many as three files can be picked
   % [filename_seed,pathname_seed]=uigetfile({'*.*','All Files'},...
       % 'Select Seed File...','MultiSelect','on');
       %##### Modify this for get another seed
    if seed_name==1
        filename_seed='seedROI_DMN_PCC_Pievani2017JAD.nii';
    else
        filename_seed='amygdalar_dx_ROI.nii';
    end
    pathname_seed='/home/sbombieri/ROIs/';
    disp(filename_seed);
    if ~isequal(filename_seed,0) || ~isequal(pathname_seed,0)
   
        % In the case in which only one seed file is entered ('char' type),
        % we convert into an array
        if ischar(filename_seed)
            filename_seed = {filename_seed};
        end
        disp(filename_seed);
        % If we enter that statement, it means that we have only one seed type
        % across subjects (that is, no subject-specific data)
        if length(filename_seed) <= 3

            % Number of entered seeds
            handles.n_seed = length(filename_seed);

            for myindex = 1:length(filename_seed)

                disp(['Processing seed number ',num2str(myindex),'...']);

                File_seed = fullfile(pathname_seed, filename_seed{myindex});

                % The NIFTI file for the seed is read, converted into the
                % proper resolution, masked and made logical
                tmp_seedh = spm_vol(File_seed);
                tmp_seed = spm_read_vols(tmp_seedh);

                tmp_seed2 = CAP_V2V(tmp_seed,tmp_seedh.dim,...
                tmp_seedh.mat,handles.brain_info{1}.dim,handles.brain_info{1}.mat);

                tmp = tmp_seed2(:);
                tmp = logical(tmp(handles.mask));

                % If the file is of suitable dimensions
                if islogical(tmp) && size(tmp,2) == 1 && size(tmp,1) == sum(handles.mask)

                    % Then we put it in the handles, enable the plotting button, and
                    % make the seed selection button green
                    handles.seed(:,myindex) = tmp;

                  
                else
                    errordlg('The file you entered appears to be of wrong dimensions...');
                    
                end

            end
        end
    else
        disp('## SELECT CORRECT SEED ##');
    end
   Seed=handles.seed;
   
   %%########## Multilple seeds version ###############
else  
    fprintf("########################\n ->> Multilple seeds version \n##################");
    %disp(handles.brain_info{1}.dim);
    
    %disp(handles.brain_info{1}.mat);
    % Multiselection is on, so that as many as three files can be picked
   % [filename_seed,pathname_seed]=uigetfile({'*.*','All Files'},...
       % 'Select Seed File...','MultiSelect','on');
    filename_seed={'seedROI_DMN_PCC_Pievani2017JAD.nii','amygdalar_dx_ROI.nii'};
    pathname_seed='/home/sbombieri/ROIs/';
    disp(filename_seed);
    if ~isequal(filename_seed,0) || ~isequal(pathname_seed,0)
   
        % In the case in which only one seed file is entered ('char' type),
        % we convert into an array
        if ischar(filename_seed)
            filename_seed = {filename_seed};
        end
        disp(filename_seed);
        % If we enter that statement, it means that we have only one seed type
        % across subjects (that is, no subject-specific data)
        if length(filename_seed) <= 3

            % Number of entered seeds
            handles.n_seed = length(filename_seed);

            for myindex = 1:length(filename_seed)

                disp(['Processing seed number ',num2str(myindex),'...']);

                File_seed = fullfile(pathname_seed, filename_seed{myindex});

                % The NIFTI file for the seed is read, converted into the
                % proper resolution, masked and made logical
                tmp_seedh = spm_vol(File_seed);
                tmp_seed = spm_read_vols(tmp_seedh);

                tmp_seed2 = CAP_V2V(tmp_seed,tmp_seedh.dim,...
                tmp_seedh.mat,handles.brain_info{1}.dim,handles.brain_info{1}.mat);

                tmp = tmp_seed2(:);
                tmp = logical(tmp(handles.mask));

                % If the file is of suitable dimensions
                if islogical(tmp) && size(tmp,2) == 1 && size(tmp,1) == sum(handles.mask)

                    % Then we put it in the handles, enable the plotting button, and
                    % make the seed selection button green
                    handles.seed(:,myindex) = tmp;

                  
                else
                    errordlg('The file you entered appears to be of wrong dimensions...');
                    
                end

            end
        end
    else
        disp('## PLEASE SELECT CORRECT SEED ##');
    end
   Seed=handles.seed; 
end

end
