%% This is an example script to run the CAPs analyses without the GUI
% In this script, we assume one population of subjects only
%% The Good One
clear all;
clc;
tt1=datetime(); %start time
savedir='F:\Tirocinio_HDD\Data';
fprintf('\n######## Start ########\nPART II\n---- Reading folders ----\n');
%opts=delimitedTextImportOptions('LeadingDelimitersRule','ignore');
%Data_OI=readtable('F:\Tirocinio_HDD\Data\ADNI\subjectlist.txt',opts);
%###################
%trovare il modo di dividere i sani dai malati, forse nel python che genera
%subjectlist? magari creando 2 subjectlist, una per i sani e una per i
%malati
%###################
%Data_OI=importdata('F:\Tirocinio_HDD\Data\ADNI\subjectlist.txt'); %for overall/final evaluation
%disp('Subjects in this run: ');
%disp(Data_OI);
%Data_OI = spm_select(Inf,'dir','Select the directories containing the functional data to analyse within the assessed group...');   
%fprintf("* For a total of %i subjects.\n",size(Data_OI,1));
%DATASET_Total = prepareData(Data_OI,'all'); %for overall/final evaluation

% Be sure to provide to these two the correct subjectlist.txt, one for
% patients and one for healty controls:

Data_OI_Patients=importdata('F:\Tirocinio_HDD\Data\ADNI\subjectlist_PATIENTS.txt');
disp('>> List of Patients <<');
disp(Data_OI_Patients);
Data_OI_Controls=importdata('F:\Tirocinio_HDD\Data\ADNI\subjectlist_CONTROLS.txt');
disp('>> List of Controls <<');
disp(Data_OI_Controls);

fprintf("* For a total of %i subjects.\n",size(Data_OI_Patients,1)+size(Data_OI_Controls,1));
%1=Patients, 0=Healty or onepop | 1=onepop, 2=twopop, 3=all dataset for final
DATASET_HEALTY=prepareData(Data_OI_Controls,'healty');
disp('Controls import complete.');

DATASET_PATIENTS=prepareData(Data_OI_Patients,'patients',DATASET_HEALTY.mask); %serve solo la maschera, come nella GUI una va bene per tutti, in realtà non è maschera solo dei sani ma maschera in generale, in questo momento la si salva nella struct dei sani, in futuro da cambiare e fare a parte
disp('Patients import complete.');



%test
%return
%end test, comment when ready

%solo su onepop
%DATASET.ROI_INFO{1}=load_nii('F:\Tirocinio_HDD\Data\ADNI\ROIs\seedROI_DMN_PCC_Pievani2017JAD.nii');%

%SEED=setSeed(DATASET_HEALTY,2); %alla fine del dataset gli serve solo brain_info e la maschera

UnivrMetrics_PATIENTS={};
UnivrMetrics_HEALTY={};
UnivrParameters={};
UnivrSpatioTemporalSelection={};

disp('All data imported correctly.');
disp('---- End of loading phase ----');


%% 1. Loading the data files

% Data: cell array, each cell of size n_TP x n_masked_voxels
%DATASET_HEALTY.TC = DATASET_HEALTY.TC;

%DATASET_PATIENTS.TC = DATASET_PATIENTS.TC;

% Mask: n_voxels x 1 logical vector
mask = DATASET_HEALTY.mask;

% Header: the header (obtained by spm_vol) of one NIFTI file with proper
% data dimension and .mat information
brain_info = DATASET_HEALTY.brain_info;

% Framewise displacement: a n_TP x n_subj matrix with framewise
% displacement information

FD1 = DATASET_HEALTY.FD{1};

FD2 = DATASET_PATIENTS.FD{1};

% Seed: a n_masked_voxels x n_seed logical vector with seed information
Seed =  setSeed(DATASET_HEALTY,2);

fprintf("\n -data loaded correctly\n");

%% 2. Specifying the main parameters

% Threshold above which to select frames
T = 15;

% Selection mode ('Threshold' or 'Percentage')
SelMode = 'Percentage';

% Threshold of FD above which to scrub out the frame and also the t-1 and
% t+1 frames (if you want another scrubbing setting, directly edit the
% code)
Tmot = 0.1; %prima 0.5

% Type of used seed information: select between 'Average','Union' or
% 'Intersection'
SeedType = 'Average';

% Contains the information, for each seed (each row), about whether to
% retain activation (1 0) or deactivation (0 1) time points
SignMatrix = [1,0];
if(SignMatrix(1)==1)
    ret='Activation';
    disp('-> CAPs: selected ACTIVATION');
else
    ret='Deactivation';
    disp('-> CAPs: selected DEACTIVATION');
end
% Percentage of positive-valued voxels to retain for clustering
Pp = 100;

% Percentage of negative-valued voxels to retain for clustering
Pn = 100;

% Number of repetitions of the K-means clustering algorithm
n_rep = 50;

% Percentage of frames to use in each fold of consensus clustering
Pcc = 80;

% Number of folds we run consensus clustering for
N = 25;

%Save to our struct
UnivrParameters.mainParameters.T=T;
UnivrParameters.mainParameters.SelectBy=SelMode;
UnivrParameters.mainParameters.Tmot=Tmot;
UnivrParameters.mainParameters.SeedType=SeedType;
UnivrParameters.mainParameters.TypeOfRetained=ret;
UnivrParameters.mainParameters.Pp=Pp;
UnivrParameters.mainParameters.Pn=Pn;
UnivrParameters.mainParameters.n_repetition_Kmeans=n_rep;
UnivrParameters.mainParameters.Pcc=Pcc;
UnivrParameters.mainParameters.N=N;
disp('-> Parameters ok. ');
%% 2.1 Recap parameters
% show a recap of inserted parameters.
disp("*----------*");
fprintf(1,"Current parameters:\n->Tmot (M)=%f\n->T=%d\n->Number of repetitions of the K-means clustering algorithm= %d\n->Percentage of positive-valued voxels to retain for clustering= %d\n->Percentage of negative-valued voxels to retain for clustering= %d\n->Percentage of frames to use in each fold of consensus clustering= %d\n",Tmot,T,n_rep,Pp,Pn,Pcc);
disp("*----------*");

%% 3. Selecting the frames to analyse    
fprintf('\n->Selecting the frames to be analyzed.\n');
% Xon will contain the retained frames, and Indices will tag the time
% points associated to these frames, for each subject (it contains a
% subfield for retained frames and a subfield for scrubbed frames)
[Xon1,~,Indices1] = CAP_find_activity(DATASET_HEALTY.TC,Seed,T,FD1,Tmot,SelMode,SeedType,SignMatrix);
[Xon2,~,Indices2] = CAP_find_activity(DATASET_PATIENTS.TC,Seed,T,FD2,Tmot,SelMode,SeedType,SignMatrix);




fprintf('-> Done <-\n');
%% 4 Finding plausible number of clusters
fprintf("\n PART III \n");
%% 4.1 Consensus clustering (if wished to determine the optimum K)
fprintf("Clustering\n");
% This specifies the range of values over which to perform consensus
% clustering: if you want to run parallel consensus clustering processes,
% you should feed in different ranges to each call of the function
% K_range = 2:20; %original
time1=datetime();
disp('Finding correct number K for clusters');
UnivrPlots={};

consensusk=0; %skip consensus if set to 0 

K_range=6; %max number of cluster you want to check 
% Have each of these run in a separate process on the server =)
%[Consensus] = CAP_ConsensusClustering(Xon1,K_range,'items',Pcc/100,N,'correlation');  %Original


if (consensusk~=0)
    % This specifies the range of values over which to perform consensus
    % clustering: if you want to run parallel consensus clustering processes,
    % you should feed in different ranges to each call of the function
    %K_range = 2:3;
    UnivrPlots.ConsensusPeformed='True';
    fprintf('## Consensus ##\n');
    t1=datetime();
    
    
    % Have each of these run in a separate process on the server =)
    %Viene lanciato solo su i sani in questo caso
    subsampleFraction=Pcc/100; %original
    fprintf(1,'Checking max number of cluster = %d\n',K_range);
    fprintf(1,'Performing consensus on %.f%% of total.\n',Pcc);
   % cluster=parcluster("LocalProfile1");
    %parfor ( i = 1:(K_range-1),cluster)
    %    [Consensus(:,:,i)] = CAP_ConsensusClustering(Xon1,i+1,'items',subsampleFraction,N,'correlation');
   % end
   for i = 1:(K_range-1)
        [Consensus(:,:,i)] = CAP_ConsensusClustering(Xon1,i+1,'items',subsampleFraction,N,'correlation');
    end
   
   
    t2=datetime();
    fprintf(1,'--> Consensus Clustering took: %s\n',t2-t1);
    clear subsampleFraction;
    clear t1;
    clear t2;
% Calculates the quality metrics
% [~,Qual] = ComputeClusteringQuality(Consensus,[]); %this before
    [CDF,Qual] = ComputeClusteringQuality(Consensus,[]);
    %bar(2:K_range,1-Qual);
% Qual should be inspected to determine the best cluster number(s)
%##### 4.1.1 Inspect Qual matrix #####
 % 
  %  bar(2:K_range,1-Qual);
    % ##### 4.1.2 Silouhette #####
    UnivrPlots.Silhouette.SilhouettePerformed='True';
    [eva,indices,sampledFrames,Silh]=compute_silhouette(Xon1,K_range,Pcc/100,N);
    UnivrPlots.Silhouette.Silhouette=eva;
    %UnivrPlots.Silhouette.Silhouette2=eva2;
    UnivrPlots.Silhouette.Indices=indices;
    UnivrPlots.Silhouette.sampledData=sampledFrames;
    %K_opt = nan;
    %close all;

    UnivrPlots.Consensus.Consensus=Consensus;
    UnivrPlots.Consensus.Quality=Qual;
else
    disp("Consensus clustering skipped");
    K_opt=5; %number of default cluster if you skipped this section
    UnivrPlots.Consensus.ConsensusPeformed='False';
    UnivrPlots.Silhoulette.SilhouettePerformed='False';
end
time2=datetime();
fprintf(1,'Total time CONSENSUS: %s\n',time2-time1);

%% 4.2 Inspecting both 
if (consensusk~=0)



    close all;
  
    univrPlots(K_range,UnivrPlots,Qual,eva);

    UnivrPlots.K_range=K_range;

    fprintf(1,"--->Suggested K from silhouette (evalclusters) = %d",eva.OptimalK);
    %Ask user number K
    prompt="Select the number K of clusters you want: ";
    K_opt=input(prompt);
    UnivrSpatioTemporalSelection.NumberOfClusters=K_opt;
    clear custom_cm;
    clear tmp_plot;
    clear prompt;
    close all;
% You should fill this with the actual value 
%K_opt = 6;

end


%% 5. Clustering into CAPs
clear Seed;
fprintf('--> number of clusters selected: %i \n',K_opt);
fprintf('\n## Clustering into CAPs... ##\n');
[CAP,~,~,idx1,CorrDist] = Run_Clustering(cell2mat(Xon1),...
        K_opt,mask,brain_info{1},Pp,Pn,n_rep,[],SeedType);

CAPs={};
CAPs.CAP_HEALTY=CAP;
CAPs.idx_HEALTY=idx1;
UnivrSpatioTemporalSelection.CAPs=CAPs;
    
%% 6. Assignment of the frames from population 2

% Parameter that governs the stringency of assignment: if Ap = 5%, we
% assign a frame to a CAP if spatial correlation exceeds the 5th percentile
% of the distribution of spatial correlations between the CAP, and its
% constituting frames
Ap = 5;

idx2 = CAP_AssignFrames(CAP,cell2mat(Xon2),CorrDist,Ap)';
CAPs.CAP_PATIENTS=CAP; %messe io
CAPs.idx_PATIENTS=idx2; %messe io



%% 7. Computing metrics
fprintf('--> Computing metrics...\n');
% The TR of your data in seconds
TR = 3;

[ExpressionMap1,Counts1,Entries1,Avg_Duration1,Duration1,TransitionProbabilities1,...
    From_Baseline1,To_Baseline1,Baseline_resilience1,Resilience1,Betweenness1,...
    InDegree1,OutDegree1,SubjectEntries1] = Compute_Metrics_simpler(idx1,...
    Indices1.kept.active,Indices1.scrubbedandactive,K_opt,TR);

[ExpressionMap2,Counts2,Entries2,Avg_Duration2,Duration2,TransitionProbabilities2,...
    From_Baseline2,To_Baseline2,Baseline_resilience2,Resilience2,Betweenness2,...
    InDegree2,OutDegree2,SubjectEntries2] = Compute_Metrics_simpler(idx2,...
    Indices2.kept.active,Indices2.scrubbedandactive,K_opt,TR);

fprintf("DONE metrics.\n");
%% 7.1 Save metrics & other data
fprintf('Saving metrics...\n');
%PATIENTS
UnivrMetrics_PATIENTS.ExpressionMap=ExpressionMap2;
UnivrMetrics_PATIENTS.Counts=Counts2;
UnivrMetrics_PATIENTS.Entries=Entries2;
UnivrMetrics_PATIENTS.Avg_Duration=Avg_Duration2;
UnivrMetrics_PATIENTS.Duration=Duration2;
UnivrMetrics_PATIENTS.TransitionProbabilities=TransitionProbabilities2;
UnivrMetrics_PATIENTS.From_Baseline=From_Baseline2;
UnivrMetrics_PATIENTS.To_Baseline=To_Baseline2;
UnivrMetrics_PATIENTS.Baseline_resilience=Baseline_resilience2;
UnivrMetrics_PATIENTS.Resilience=Resilience2;
UnivrMetrics_PATIENTS.Betweenness=Betweenness2;
UnivrMetrics_PATIENTS.InDegree=InDegree2;
UnivrMetrics_PATIENTS.OutDegree=OutDegree2;
UnivrMetrics_PATIENTS.SubjectEntries=SubjectEntries2;
UnivrMetrics_PATIENTS.TR=TR;

%GENERAL
UnivrMetrics.TR=TR;
UnivrParameters.maxNumberOfClusters=K_range;
UnivrParameters.NumberOfClusters=K_opt;
UnivrParameters.Ap=Ap;
UnivrParameters.CorrDist=CorrDist;

UnivrData={};
%UnivrData.Dataset=DATASET_Total;
%UnivrData.Seed=Seed;
UnivrData.FD_PATIENTS=FD2;
UnivrData.FD_HEALTY=FD2;
UnivrData.brain_info=brain_info;
%UnivrData.mask=mask;
%UnivrData.TC_PATIENTS=DATASET_PATIENTS.TC;
%UnivrData.TC_HEALTY=DATASET_HEALTY.TC;
UnivrData.SeedType=SeedType;
UnivrData.ListMaxNumberOfClusters=K_range;

UnivrSpatioTemporalSelection.Indices_HEALTY=Indices1;
UnivrSpatioTemporalSelection.Indices_PATIENTS=Indices2;
%UnivrData.Subjects=Data_OI; %da aggiungere Patients+Controls

UnivrParameters.SaveFolder=savedir;


%HEALTY CONTROLS
UnivrMetrics_HEALTY.TR=TR;
UnivrMetrics_HEALTY.ExpressionMap=ExpressionMap1;
UnivrMetrics_HEALTY.Counts=Counts1;
UnivrMetrics_HEALTY.Entries=Entries1;
UnivrMetrics_HEALTY.Avg_Duration=Avg_Duration1;
UnivrMetrics_HEALTY.Duration=Duration1;
UnivrMetrics_HEALTY.TransitionProbabilities=TransitionProbabilities1;
UnivrMetrics_HEALTY.From_Baseline=From_Baseline1;
UnivrMetrics_HEALTY.To_Baseline=To_Baseline1;
UnivrMetrics_HEALTY.Baseline_resilience=Baseline_resilience1;
UnivrMetrics_HEALTY.Resilience=Resilience1;
UnivrMetrics_HEALTY.Betweenness=Betweenness1;
UnivrMetrics_HEALTY.InDegree=InDegree1;
UnivrMetrics_HEALTY.OutDegree=OutDegree1;
UnivrMetrics_HEALTY.SubjectEntries=SubjectEntries1;


fprintf("\nDONE saving metrics.\n");

%% 7. Save CAPs nifti
%Not present in original script (missing)
%adapting from CAP_TB.m main file 

fprintf('-->Now saving to nifti\n');

savename=datestr(datetime('today'));
fprintf(1,"\nOutput images will have names as %s\n",savename);
CAPToNIFTI(CAP,mask,brain_info{1},savedir,savename);

disp('->Completed converting to Nifti.');

%% 8. Save CAPs for every subject

savedir_subjects='F:\Tirocinio_HDD\Data\Results\';
mkdir(savedir_subjects);
%getCAPS_forEachSubject();
%     -idx : uno tra idx1 o idx2
%     -DATASET: uno tra _PATIENTS o _HEALTY
%     -Data_OI: come DATASET
%     -savedir: directory where to save files
%     -Xon: Xon1 o Xon2

%PATIENTS
getCAPS_forEachSubject_twopop(idx2,DATASET_PATIENTS,Data_OI_Patients,savedir_subjects,Xon2,K_opt,"patients",CAP,mask,brain_info{1});

%HEALTY CONTROLS
getCAPS_forEachSubject_twopop(idx1,DATASET_HEALTY,Data_OI_Controls,savedir_subjects,Xon1,K_opt,"healty",CAP,mask,brain_info{1});



%disp("-> For plot consensus and silhouette type 'eval(UnivrPlots.figure);' <-");
fprintf('\n-> Figure of silhouette & consensus clustering is saved in the current directory, named "Max_clust.fig" <-');

fprintf(1,"\n######## All ended correctly ########\n");
%% 9
%Save data to current directory
fprintf('\n###### Saving data to RESULTS.mat #####\n');
Results=struct();
Results.Data=UnivrData;
clear UnivrData;
Results.Metrics=UnivrMetrics;
clear UnivrMetrics;
Results.Plots=UnivrPlots;
clear UnivrPlots;
Results.SpatioTemporalSelection=UnivrSpatioTemporalSelection;
clear UnivrSpatioTemporalSelection;
Results.Parameters=UnivrParameters;
Results.SavedOnDate=datetime(now,'ConvertFrom','datenum');
tt2=datetime();
totTime=tt2-tt1;
Results.TotalTime=totTime;
save(fullfile(savedir,'RESULTS.mat'),'Results','-v7.3');
fprintf('Results will available at: %s\n',savedir);



fprintf(1,'\n############All process took: %s ############\n',totTime);





%% 10 Clear workspace
clear savedir_subjects;
close all;
%clear Results;
clear savedir;
disp("Cleaning workspace");
clear CAP;
clear idx1;
clear ExpressionMap1;
clear Counts1;
clear Entries1;
clear Avg_Duration1;
clear Duration1;
clear TransitionProbabilities1;
clear From_Baseline1;
clear To_Baseline1;
clear Baseline_resilience1;
clear Resilience1;
clear Betweenness1;
clear InDegree1;
clear OutDegree1;
clear SubjectEntries1;

clear idx2;
clear ExpressionMap2;
clear Counts2;
clear Entries2;
clear Avg_Duration2;
clear Duration2;
clear TransitionProbabilities2;
clear From_Baseline2;
clear To_Baseline2;
clear Baseline_resilience2;
clear Resilience2;
clear Betweenness2;
clear InDegree2;
clear OutDegree2;
clear SubjectEntries2;
clear K_opt;

clear TR;
clear Tmot;
clear mask;
clear brain_info;
clear DATASET_HEALTY.TC;
clear DATASET_PATIENTS.TC;
clear FD1;
clear FD2;

clear Silh;
clear savename;
clear idx;
clear varname;
clear Xon;
clear xontmp;
clear tosavedir;
clear sub;
clear sub_tmp;
clear SignMatrix;
clear SelMode;
clear Pcc;
clear Pn;
clear n_rep;
clear N;
clear m;
clear meanxontmp;
clear i;
clear FD;
clear cl;
clear brain;
clear K_Range;
clear T;
clear Pp;
clear CAPs;
clear Indices;
clear consensusk;
clear idx_sub;
clear Data_OI;
clear K_range;
clear TC;
clear SeedType;
clear DATASET;
clear ret;
clear OUT_SUBJECTCAPS;
clear Qual;
clear eva;
clear Consensus;
clear s;
clear Ap;
clear CDF;
clear tt1;
clear tt2;
clear totTime;
clear Xon1;
clear Xon2;
clear time1;
clear time2;
%clear sampledFrames;
clear CorrDist;
clear Data_OI_Controls;
clear Data_OI_Patients;
clear UnivrMetrics;
clear UnivrMetrics_HEALTY;
clear UnivrMetrics_PATIENTS;
clear Indices1;
clear Indices2;
clear Seed;
%clear DATASET_Total;



















