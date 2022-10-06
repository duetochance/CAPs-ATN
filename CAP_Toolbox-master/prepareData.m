%% Prepare data for computation in terms of getting directories 
% with the volumes already in 3D format
% Data_OI char array of directories where is possible to locate subjects
% and volumes
%  ARGS:
%
%   - Data_OI : list of all subjects created by unfold&rename.py
%   
%   - type : [1] or [0]; 1=Patients, 0=Healty or onepop
%   
%   
%   data investiagation

%% The Good One
function [DATASET] = prepareData(Data_OI,typecomp)
%#############################################################    HEALTY CONTROLS    #############################################################
imax=size(Data_OI,1);
prog=0;
CAP_DIR='/home/sbombieri/CODE/Code/CAP_Toolbox-master';
if strcmp(typecomp,'healty')==1
    %Consts

    fprintf(1,'\n##### Preparing data for HEALTY subjects #####\n')
    % preliminary preps
    % Loads and sets the brain underlay used for plotting purposes
    handles.prefix='vol'; %our new prefix
    Underlay = load_nii('Underlay.nii');
    Underlay_mat = [Underlay.hdr.hist.srow_x; Underlay.hdr.hist.srow_y; Underlay.hdr.hist.srow_z; 0 0 0 1];
    Underlay_dim = Underlay.hdr.dime.dim;
    Underlay_dim = Underlay_dim(2:4);
    handles.Underlay_info.dim = Underlay_dim;
    handles.Underlay_info.mat = Underlay_mat;
    clear Underlay;
    clear Underlay_dim;
    clear Underlay_mat;
    load('brain.mat');
    assignin('base','brain', brain);
    handles.brain = brain;

    handles.n_datasets = 0;
    %

    DATASET = {};
    handles.n_subjects = size(Data_OI,1);
    
    if strcmp(mask,'y')==1
        %%%  MASK

        ToMask = Data_OI(1,:);

        % Current file path
        tmp = ToMask(1,:);

        % Selects all the functional files that match our criteria
        % (prefix 'vol' meaning MNI space data)
        FFiles = cellstr(spm_select('List',fullfile(tmp),['^' handles.prefix '.*\.' 'nii' '$']));
        fprintf('Masking for controls resolution..\n');
        % We read the header of the first one to get the data
        % resolution
        a = spm_vol(fullfile(CAP_DIR,'DefaultData','Default_mask.nii'));
        handles_brain_info = spm_vol(fullfile(tmp,FFiles{1}));
        DATASET.brain_info=handles_brain_info;

        b = spm_read_vols(a);
        b(b < 0.9) = 0;
        b(b >= 0.9) = 1;

        maskf = CAP_V2V(b,a.dim,a.mat,handles_brain_info{1}.dim,handles_brain_info{1}.mat);


        %check mask, save data for now...
        DATASET.MASK_INFO.a=a;
        DATASET.MASK_INFO.MASK=maskf;
        DATASET.MASK_INFO.b=b;

        %fprintf('\n####### SAVED MASK INFOS #########\n')



        % Filling accordingly
        handles.mask = logical(maskf(:));
        DATASET.mask=handles.mask;
    end
    %%% DATA
    DATASET.numberOfSubjects=size(Data_OI,1);
    disp("Loading data");
    fprintf(1,"Completeness:     0%%");
    try
            % We now want to update the FD and TC variables by going through all
            % the subjects to add to the considered group...
            for i = 1:size(Data_OI,1)
                prog = ( 100*(i/imax) );
                
                %disp(['Currently preparing the data from run ',num2str(i),'...']);
                %fprintf(' \n ')
                % As before, the "vol" prefix is looked for
                FFiles = cellstr(spm_select('List',fullfile(Data_OI(i,:)),['^' handles.prefix '.*\.' 'nii' '$']));

                % Functional files are read one after the other to build
                % tmp_data
                tmp_data = [];

                for t = 1:length(FFiles)

                    tmp1 = spm_vol(fullfile(Data_OI(i,:),FFiles{t}));
                    tmp2 = spm_read_vols(tmp1{1});
                    tmp3 = tmp2(:);

                    tmp_data = [tmp_data;tmp3(handles.mask)'];
                end

                % Z-scoring is performed within the toolbox
                %tmp_data = detrend(tmp_data);
                tmp_data = zscore(tmp_data);

                % The ready-to-analyse data is put in TC
                handles.TC{i} = tmp_data;
                clear tmp_data
                fprintf(1,' \b\b\b\b%3.0f%% ',prog);
                try
                    % Look for the text file with motion parameters (should be the
                    % first text file found)
                    MFile = cellstr(spm_select('List',fullfile(Data_OI(i,:)),['.*\.' 'txt' '$']));
                    %disp('FD file: ',MFile);
                    % Computes framewise displacement and fills the FD matrix
                    % accordingly
                    file=cell2mat(fullfile(Data_OI(i,:),MFile{1}));
                    handles.FD{1}(:,i) = univr_CAP_ComputeFD(file);
                   % handles.FD(:,i) = CAP_ComputeFD(fullfile(Data_OI(i,:),MFile{1}));

                catch error
                    handles.FD{1}(:,i) = zeros(length(FFiles),1);
                    fprintf(2,'Identifier: %s\n',error.identifier);
                    fprintf(2,'Message from error: %s\n',error.message);
                    errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    error.stack(1).name, error.stack(1).line, error.message);
                    fprintf(2, '%s\n', errorMessage);
                    disp(['Error at run ',num2str(i),'...']);
                    disp('Could not process motion text file; assuming zero movement...');
                end

            try
                if handles.n_datasets == 0
                        DATASET.SubjSize.VOX = size(handles.TC{1},2);
                        DATASET.SubjSize.TP = size(handles.TC{1},1);
                end
            catch error
                    disp('#############');
                 errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    error.stack(1).name, error.stack(1).line, error.message);
                    fprintf(2, '%s\n', errorMessage);
                    disp(['Error at run ',num2str(i),'...']);
            end

            % If we are loading the first dataset, we convert the underlay
            % to the resolution of the functional data for plotting


                % The brain variable now contains a good resolution
                % underlay that can directly be overlapped with the
                % functional data
                % Loads and sets the brain underlay used for plotting purposes

                handles.brain = CAP_V2V(handles.brain,handles.Underlay_info.dim,handles.Underlay_info.mat,handles_brain_info{1}.dim,handles_brain_info{1}.mat);

            end

    catch e
            errordlg('Could not load subject population!');
            fprintf(2,'The identifier was:\n%s',e.identifier);
            fprintf(2,'There was an error! The message was:\n%s',e.message);
            disp(['Error at run ',num2str(i),'...']);
    end %end try
       % STORE DATA TO DATASET STRUCT

                    DATASET.FD=handles.FD;
                    DATASET.SubjSize.VOX = size(handles.TC{1},2);
                    DATASET.SubjSize.TP = size(handles.TC{1},1);
                    DATASET.TC=handles.TC;
                    DATASET.brain=handles.brain;
                    clear handles;
                    clear tmp;
                    clear tmp1; clear tmp2; clear tmp3;
                    
end
disp("    END Loading PATIENTS    ");
    % #############################################################    PATIENTS   #############################################################
    
if  strcmp(typecomp,'patients')==1     
    %Consts

    fprintf(1,'\n##### Preparing data for PATIENTS #####\n')
    % preliminary preps
    % Loads and sets the brain underlay used for plotting purposes
    handles.prefix='vol'; %our new prefix
    Underlay = load_nii('Underlay.nii');
    Underlay_mat = [Underlay.hdr.hist.srow_x; Underlay.hdr.hist.srow_y; Underlay.hdr.hist.srow_z; 0 0 0 1];
    Underlay_dim = Underlay.hdr.dime.dim;
    Underlay_dim = Underlay_dim(2:4);
    handles.Underlay_info.dim = Underlay_dim;
    handles.Underlay_info.mat = Underlay_mat;
    clear Underlay;
    clear Underlay_dim;
    clear Underlay_mat;
    load('brain.mat');
    assignin('base','brain', brain);
    handles.brain = brain;
    handles.mask=mask;
    handles.n_datasets = 0;
    %

    DATASET = {};
    handles.n_subjects = size(Data_OI,1);
    %%%  MASK

    ToMask = Data_OI(1,:);

    % Current file path
    tmp = ToMask(1,:);

    % Selects all the functional files that match our criteria
    % (prefix 'vol' meaning MNI space data)
    %FFiles = cellstr(spm_select('List',fullfile(tmp),['^' handles.prefix '.*\.' 'nii' '$']));
    fprintf('Masking for patients resolution..\n');
    % We read the header of the first one to get the data
    % resolution
    %a = spm_vol(fullfile(CAP_DIR,'DefaultData','Default_mask.nii'));
    %handles_brain_info = spm_vol(fullfile(tmp,FFiles{1}));
    %DATASET.brain_info=handles_brain_info;

    %b = spm_read_vols(a);
    %b(b < 0.9) = 0;
    %b(b >= 0.9) = 1;

    %maskf = CAP_V2V(b,a.dim,a.mat,handles_brain_info{1}.dim,handles_brain_info{1}.mat);


    %check mask, save data for now...
    %DATASET.MASK_INFO.a=a;
    %DATASET.MASK_INFO.MASK=maskf;
    %DATASET.MASK_INFO.b=b;

    %disp('####### SAVED MASK INFOS #########')



    % Filling accordingly
    %handles.mask = logical(maskf(:));
    %DATASET.mask=handles.mask;
    %%% DATA
    DATASET.numberOfSubjects=size(Data_OI,1);
    disp("Loading data");
    disp("Completeness:     0%");
    try
            % We now want to update the FD and TC variables by going through all
            % the subjects to add to the considered group...
            for i = 1:size(Data_OI,1)
                %fprintf(' \n ')
                prog = ( 100*(i/imax) );
                
                %disp(['Currently preparing the data from run ',num2str(i),'...']);
                fprintf(' \n ')
                % As before, the "vol" prefix is looked for
                FFiles = cellstr(spm_select('List',fullfile(Data_OI(i,:)),['^' handles.prefix '.*\.' 'nii' '$']));

                % Functional files are read one after the other to build
                % tmp_data
                tmp_data = [];

                for t = 1:length(FFiles)

                    tmp1 = spm_vol(fullfile(Data_OI(i,:),FFiles{t}));
                    tmp2 = spm_read_vols(tmp1{1});
                    tmp3 = tmp2(:);

                    tmp_data = [tmp_data;tmp3(handles.mask)'];
                end

                % Z-scoring is performed within the toolbox
                %tmp_data = detrend(tmp_data);
                tmp_data = zscore(tmp_data);

                % The ready-to-analyse data is put in TC
                handles.TC{i} = tmp_data;
                clear tmp_data
                fprintf(1,' \b\b\b%3.0f%% ',prog);
                try
                    % Look for the text file with motion parameters (should be the
                    % first text file found)
                    MFile = cellstr(spm_select('List',fullfile(Data_OI(i,:)),['.*\.' 'txt' '$']));
                    %disp('FD file: ',MFile);
                    % Computes framewise displacement and fills the FD matrix
                    % accordingly
                    file=cell2mat(fullfile(Data_OI(i,:),MFile{1}));
                    handles.FD{1}(:,i) = univr_CAP_ComputeFD(file);
                   % handles.FD(:,i) = CAP_ComputeFD(fullfile(Data_OI(i,:),MFile{1}));

                catch error
                    handles.FD{1}(:,i) = zeros(length(FFiles),1);
                    fprintf(1,'Identifier: %s\n',error.identifier);
                    fprintf(1,'Message from error: %s\n',error.message);
                    errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    error.stack(1).name, error.stack(1).line, error.message);
                    fprintf(1, '%s\n', errorMessage);
                    disp(['Error at run ',num2str(i),'...']);
                    disp('Could not process motion text file; assuming zero movement...');
                end

            try
                if handles.n_datasets == 0
                        DATASET.SubjSize.VOX = size(handles.TC{1},2);
                        DATASET.SubjSize.TP = size(handles.TC{1},1);
                end
            catch error
                    disp('#############');
                 errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    error.stack(1).name, error.stack(1).line, error.message);
                    fprintf(2, '%s\n', errorMessage);
                    disp(['Error at run ',num2str(i),'...']);
            end

            % If we are loading the first dataset, we convert the underlay
            % to the resolution of the functional data for plotting


                % The brain variable now contains a good resolution
                % underlay that can directly be overlapped with the
                % functional data
                % Loads and sets the brain underlay used for plotting purposes

                %handles.brain = CAP_V2V(handles.brain,handles.Underlay_info.dim,handles.Underlay_info.mat,handles_brain_info{1}.dim,handles_brain_info{1}.mat);

            end

    catch e
            errordlg('Could not load subject population!');
            fprintf(2,'The identifier was:\n%s',e.identifier);
            fprintf(2,'There was an error! The message was:\n%s',e.message);
            disp(['Error at run ',num2str(i),'...']);
    end %end try
       % STORE DATA TO DATASET STRUCT

                    DATASET.FD=handles.FD;
                    DATASET.SubjSize.VOX = size(handles.TC{1},2);
                    DATASET.SubjSize.TP = size(handles.TC{1},1);
                    DATASET.TC=handles.TC;
                    %DATASET.brain=handles.brain;
                    clear handles;
                    clear tmp;
                    clear tmp1; clear tmp2; clear tmp3;
    
    
end

% ##############################  STEPS #######################
if  strcmp(typecomp,'step')==1     
    %Consts

    fprintf(1,'\n##### Loading portion of subjects  #####\n')
    % preliminary preps
    % Loads and sets the brain underlay used for plotting purposes
    handles.prefix='vol'; %our new prefix
    Underlay = load_nii('Underlay.nii');
    Underlay_mat = [Underlay.hdr.hist.srow_x; Underlay.hdr.hist.srow_y; Underlay.hdr.hist.srow_z; 0 0 0 1];
    Underlay_dim = Underlay.hdr.dime.dim;
    Underlay_dim = Underlay_dim(2:4);
    handles.Underlay_info.dim = Underlay_dim;
    handles.Underlay_info.mat = Underlay_mat;
    clear Underlay;
    clear Underlay_dim;
    clear Underlay_mat;
    load('brain.mat');
    assignin('base','brain', brain);
    handles.brain = brain;
    
    handles.n_datasets = 0;
    %

    DATASET = {};
    handles.n_subjects = size(Data_OI,1);
    ToMask = Data_OI(1,:);
    tmp = ToMask(1,:);
    FFiles = cellstr(spm_select('List',fullfile(tmp),['^' handles.prefix '.*\.' 'nii' '$']));
       
    % Current file path
    
    handles_brain_info = spm_vol(fullfile(tmp,FFiles{1}));
    DATASET.brain_info=handles_brain_info;
    %%%  MASK
    
      
    fprintf('\n>> Porting mask...\n')



    % Selects all the functional files that match our criteria
    % (prefix 'vol' meaning MNI space data)
    fprintf('Masking for subjects resolution..\n');
    % We read the header of the first one to get the data
    % resolution
    a = spm_vol(fullfile(CAP_DIR,'DefaultData','Default_mask.nii'));


    b = spm_read_vols(a);
    b(b < 0.9) = 0;
    b(b >= 0.9) = 1;

    maskf = CAP_V2V(b,a.dim,a.mat,handles_brain_info{1}.dim,handles_brain_info{1}.mat);


    %check mask, save data for now...
    DATASET.MASK_INFO.a=a;
    DATASET.MASK_INFO.MASK=maskf;
    DATASET.MASK_INFO.b=b;




    % Filling accordingly
    handles.mask = logical(maskf(:));
    DATASET.mask=handles.mask;
    %fprintf('\n####### SAVED MASK INFOS #########\n')

        
 
   
    %%% DATA
    DATASET.numberOfSubjects=size(Data_OI,1);
    disp("Loading data");
    disp("Completeness:     0%");
    try
            % We now want to update the FD and TC variables by going through all
            % the subjects to add to the considered group...
            for i = 1:size(Data_OI,1)
                %fprintf(' \n ')
                prog = ( 100*(i/imax) );
                
                %disp(['Currently preparing the data from run ',num2str(i),'...']);
                
                % As before, the "vol" prefix is looked for
                FFiles = cellstr(spm_select('List',fullfile(Data_OI(i,:)),['^' handles.prefix '.*\.' 'nii' '$']));

                % Functional files are read one after the other to build
                % tmp_data
                tmp_data = [];

                for t = 1:length(FFiles)

                    tmp1 = spm_vol(fullfile(Data_OI(i,:),FFiles{t}));
                    tmp2 = spm_read_vols(tmp1{1});
                    tmp3 = tmp2(:);

                    tmp_data = [tmp_data;tmp3(handles.mask)'];
                end

                % Z-scoring is performed within the toolbox
                %tmp_data = detrend(tmp_data);
                tmp_data = zscore(tmp_data);

                % The ready-to-analyse data is put in TC
                handles.TC{i} = tmp_data;
                clear tmp_data
                fprintf(1,'\b\b\b\b\b%3.0f%%\n ',prog);
                try
                    % Look for the text file with motion parameters (should be the
                    % first text file found)
                    MFile = cellstr(spm_select('List',fullfile(Data_OI(i,:)),['.*\.' 'txt' '$']));
                    %disp('FD file: ',MFile);
                    % Computes framewise displacement and fills the FD matrix
                    % accordingly
                    file=cell2mat(fullfile(Data_OI(i,:),MFile{1}));
                    handles.FD{1}(:,i) = univr_CAP_ComputeFD(file);
                   % handles.FD(:,i) = CAP_ComputeFD(fullfile(Data_OI(i,:),MFile{1}));

                catch error
                    handles.FD{1}(:,i) = zeros(length(FFiles),1);
                    fprintf(1,'Identifier: %s\n',error.identifier);
                    fprintf(1,'Message from error: %s\n',error.message);
                    errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    error.stack(1).name, error.stack(1).line, error.message);
                    fprintf(1, '%s\n', errorMessage);
                    disp(['Error at run ',num2str(i),'...']);
                    disp('Could not process motion text file; assuming zero movement...');
                end

            try
                if handles.n_datasets == 0
                        DATASET.SubjSize.VOX = size(handles.TC{1},2);
                        DATASET.SubjSize.TP = size(handles.TC{1},1);
                end
            catch error
                    disp('#############');
                 errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    error.stack(1).name, error.stack(1).line, error.message);
                    fprintf(2, '%s\n', errorMessage);
                    disp(['Error at run ',num2str(i),'...']);
            end

            % If we are loading the first dataset, we convert the underlay
            % to the resolution of the functional data for plotting


                % The brain variable now contains a good resolution
                % underlay that can directly be overlapped with the
                % functional data
                % Loads and sets the brain underlay used for plotting purposes

                %handles.brain = CAP_V2V(handles.brain,handles.Underlay_info.dim,handles.Underlay_info.mat,handles_brain_info{1}.dim,handles_brain_info{1}.mat);

            end

    catch e
            errordlg('Could not load subject population!');
            fprintf(2,'The identifier was:\n%s',e.identifier);
            fprintf(2,'There was an error! The message was:\n%s',e.message);
            disp(['Error at run ',num2str(i),'...']);
    end %end try
       % STORE DATA TO DATASET STRUCT

                    DATASET.FD=handles.FD;
                    DATASET.SubjSize.VOX = size(handles.TC{1},2);
                    DATASET.SubjSize.TP = size(handles.TC{1},1);
                    DATASET.TC=handles.TC;
                    %DATASET.brain=handles.brain;
                    clear handles;
                    clear tmp;
                    clear tmp1; clear tmp2; clear tmp3;
    
% ################################ ONEPOP #############################
% parte per scirpt con solo 1 popolazione
if strcmp(typecomp,'onepop')==1
    fprintf(" Setting up things for final evaluation...\n This may take a while, sorry :( ");
    %Consts
   
    fprintf(1,'\n##### Preparing data for onepop #####\n')
    % preliminary preps
    % Loads and sets the brain underlay used for plotting purposes
    handles.prefix='vol'; %our new prefix
    Underlay = load_nii('Underlay.nii');
    Underlay_mat = [Underlay.hdr.hist.srow_x; Underlay.hdr.hist.srow_y; Underlay.hdr.hist.srow_z; 0 0 0 1];
    Underlay_dim = Underlay.hdr.dime.dim;
    Underlay_dim = Underlay_dim(2:4);
    handles.Underlay_info.dim = Underlay_dim;
    handles.Underlay_info.mat = Underlay_mat;
    clear Underlay;
    clear Underlay_dim;
    clear Underlay_mat;
    load('brain.mat');
    assignin('base','brain', brain);
    handles.brain = brain;

    handles.n_datasets = 0;
    %

    DATASET = {};
    handles.n_subjects = size(Data_OI,1);
    %%%  MASK

    ToMask = Data_OI(1,:);

    % Current file path
    tmp = ToMask(1,:);

    % Selects all the functional files that match our criteria
    % (prefix 'vol' meaning MNI space data)
    FFiles = cellstr(spm_select('List',fullfile(tmp),['^' handles.prefix '.*\.' 'nii' '$']));
    fprintf('Porting mask to subjects resolution..');
    % We read the header of the first one to get the data
    % resolution
    a = spm_vol(fullfile(CAP_DIR,'DefaultData','Default_mask.nii'));
    handles_brain_info = spm_vol(fullfile(tmp,FFiles{1}));
    DATASET.brain_info=handles_brain_info;

    b = spm_read_vols(a);
    b(b < 0.9) = 0;
    b(b >= 0.9) = 1;

    maskf = CAP_V2V(b,a.dim,a.mat,handles_brain_info{1}.dim,handles_brain_info{1}.mat);


    %check mask, save data for now...
    DATASET.MASK_INFO.a=a;
    DATASET.MASK_INFO.MASK=maskf;
    DATASET.MASK_INFO.b=b;

    fprintf('\n####### SAVED MASK INFOS #########\n')



    % Filling accordingly
    handles.mask = logical(maskf(:));
    DATASET.mask=handles.mask;
    %%% DATA
    DATASET.numberOfSubjects=size(Data_OI,1);
    disp("Loading data");
    disp("Completeness:     0%");
    try
            % We now want to update the FD and TC variables by going through all
            % the subjects to add to the considered group...
            %fprintf(1,'\nCurrently preparing the data from run 1');
            for i = 1:size(Data_OI,1)
                %fprintf(' \n ')
                prog = ( 100*(i/imax) );
                %fprintf(1,'\b\b\b\b%3.0f%%\n \b',prog);
                
               % if i>1
                    %fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bCurrently preparing the data from run %i\n',num2str(i));
               % end
                %disp(['Currently preparing the data from run ',num2str(i),'...']);
                %fprintf(' \n ')
                % As before, the "vol" prefix is looked for
                FFiles = cellstr(spm_select('List',fullfile(Data_OI(i,:)),['^' handles.prefix '.*\.' 'nii' '$']));

                % Functional files are read one after the other to build
                % tmp_data
                tmp_data = [];

                for t = 1:length(FFiles)

                    tmp1 = spm_vol(fullfile(Data_OI(i,:),FFiles{t}));
                    tmp2 = spm_read_vols(tmp1{1});
                    tmp3 = tmp2(:);

                    tmp_data = [tmp_data;tmp3(handles.mask)'];
                end

                % Z-scoring is performed within the toolbox
                %tmp_data = detrend(tmp_data);
                tmp_data = zscore(tmp_data);

                % The ready-to-analyse data is put in TC
                handles.TC{i} = tmp_data;
                clear tmp_data
                fprintf(1,'\b\b\b\b\b%3.0f%%\n',prog);
                try
                    % Look for the text file with motion parameters (should be the
                    % first text file found)
                    MFile = cellstr(spm_select('List',fullfile(Data_OI(i,:)),['.*\.' 'txt' '$']));
                    %disp('FD file: ',MFile);
                    % Computes framewise displacement and fills the FD matrix
                    % accordingly
                    file=cell2mat(fullfile(Data_OI(i,:),MFile{1}));
                    handles.FD{1}(:,i) = univr_CAP_ComputeFD(file);
                   % handles.FD(:,i) = CAP_ComputeFD(fullfile(Data_OI(i,:),MFile{1}));

                catch error
                    handles.FD{1}(:,i) = zeros(length(FFiles),1);
                    fprintf(1,'Identifier: %s\n',error.identifier);
                    fprintf(1,'Message from error: %s\n',error.message);
                    errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    error.stack(1).name, error.stack(1).line, error.message);
                    fprintf(1, '%s\n', errorMessage);
                    disp(['Error at run ',num2str(i),'...']);
                    disp('Could not process motion text file; assuming zero movement...');
                end

            try
                if handles.n_datasets == 0
                        DATASET.SubjSize.VOX = size(handles.TC{1},2);
                        DATASET.SubjSize.TP = size(handles.TC{1},1);
                end
            catch error
                    disp('#############');
                 errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    error.stack(1).name, error.stack(1).line, error.message);
                    fprintf(2, '%s\n', errorMessage);
                    disp(['Error at run ',num2str(i),'...']);
            end

            % If we are loading the first dataset, we convert the underlay
            % to the resolution of the functional data for plotting


                % The brain variable now contains a good resolution
                % underlay that can directly be overlapped with the
                % functional data
                % Loads and sets the brain underlay used for plotting purposes

                handles.brain = CAP_V2V(handles.brain,handles.Underlay_info.dim,handles.Underlay_info.mat,handles_brain_info{1}.dim,handles_brain_info{1}.mat);

            end

    catch e
            errordlg('Could not load subject population!');
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'There was an error! The message was:\n%s',e.message);
    end %end try
       % STORE DATA TO DATASET STRUCT

                    DATASET.FD=handles.FD;
                    DATASET.SubjSize.VOX = size(handles.TC{1},2);
                    DATASET.SubjSize.TP = size(handles.TC{1},1);
                    DATASET.TC=handles.TC;
                    DATASET.brain=handles.brain;
                    clear handles;
                    clear tmp;
                    clear tmp1; clear tmp2; clear tmp3;    
    
end


end %end func





