%% The Good One
clear all;
clc;
disp("# TWO POP VERSION ");
tt1=datetime(); %start time
savedir='/home/sbombieri/RESULTS';
fprintf('\n######## Start ########\nPART II\n---- Reading folders ----\n');
Data_OI_Patients=importdata('/home/sbombieri/DATA/TauAmyFDG/subjectlist_PATIENTS.txt');
disp('>> List of Patients <<');
disp(Data_OI_Patients);
Data_OI_Controls=importdata('/home/sbombieri/DATA/TauAmyFDG/subjectlist_CONTROLS.txt');
disp('>> List of Controls <<');
disp(Data_OI_Controls);

UnivrMetrics_PATIENTS={};
UnivrMetrics_HEALTY={};
UnivrParameters={};
UnivrSpatioTemporalSelection={};
% The TR of your data in seconds
TR = 3;  %default for now...

%### SEED SELECTION OPTION
seed_num=1; %default
seed_num_type=1; %1 or 2, seed version
prompt="Select the seed:\n \t-> PCC: 1\n \t-> Amygdala: 2\n OR if y want to run with the 2 seeds version : -> 3\n";
seed_num=input(prompt);
if seed_num==3
 seed_num_type=2;
end



%###### max K
K_range=5; %max number of cluster you want to check 
prompt="\nSelect max number of clusters u looking for: (will be used as default if consensus is performed without a display) \n-> ";
K_range=input(prompt);
default_k=K_range;

%######SCREEN INFO:
is_screen=0; %if screen is present set to 1




fprintf("\nTotal subjects: %i \n",size(Data_OI_Patients,1)+size(Data_OI_Controls,1));
clear Data_OI;

%%%##### CONSENSUS
consensusk=0; %skip consensus if set to 0
prompt="Perform consensus? 1 (yes) | 0 (no)\t";
consensusk=input(prompt);


% CONTROLS
disp("CONTROLS");
%find best group (BETA)
N = size(Data_OI_Controls,1);
K = 1:ceil(sqrt(N));
D = K(rem(N,K)==0);
d = [D sort(N./D)];
%d=divisors(size(Data_OI_Controls,1));
%group=d(:,2);

%other way of divisors:
fprintf("Divisors: \n");
disp(d);
prompt="Select number of subjects for each group to load: ";
group_controls=input(prompt);


%group=5;

% PATIENTS
disp("PATIENTS");
%find best group (BETA)
N = size(Data_OI_Patients,1);
K = 1:ceil(sqrt(N));
D = K(rem(N,K)==0);
d = [D sort(N./D)];
%d=divisors(size(Data_OI_Patients,1));
%group=d(:,2);

%other way of divisors:
fprintf("Divisors: \n");
disp(d);
prompt="Select number of subjects for each group to load: ";
group_patients=input(prompt);

clear d;

%%%%%%%%%%%%Maxiter 
maxiter=300; %max iterations of kmeans

isMask=0; %0= to be calculated
run=0;


Xon1=[];
Indices1={};
Indices1.scrubbed=[];
Indices1.scrubbedandactive=[];
Indices1.kept={};
Indices1.kept.active=[];

Xon2=[];
Indices2={};
Indices2.scrubbed=[];
Indices2.scrubbedandactive=[];
Indices2.kept={};
Indices2.kept.active=[];




%% 2. Specifying the main parameters

% Threshold above which to select frames
T = 30;

% Selection mode ('Threshold' or 'Percentage')
SelMode = 'Percentage';

% Threshold of FD above which to scrub out the frame and also the t-1 and
% t+1 frames (if you want another scrubbing setting, directly edit the
% code)
Tmot = 0.1;

% Type of used seed information: select between 'Average','Union' or
% 'Intersection'
SeedType = 'Average';

% Contains the information, for each seed (each row), about whether to
% retain activation (1 0) or deactivation (0 1) time points
SignMatrix = [1,0];

% Percentage of positive-valued voxels to retain for clustering
Pp = 100;

% Percentage of negative-valued voxels to retain for clustering
Pn = 100;

% Number of repetitions of the K-means clustering algorithm
n_rep = 50;

% Percentage of frames to use in each fold of consensus clustering
Pcc = 80;

% Number of folds we run consensus clustering for
N = 50;

%% 2.1 Recap parameters
% show a recap of inserted parameters.
if run==0
    disp("*----------*");
    fprintf(1,"Current parameters:\n->Tmot (M)=%f\n->T=%d\n->Number of repetitions of the K-means clustering algorithm= %d\n->Percentage of positive-valued voxels to retain for clustering= %d\n->Percentage of negative-valued voxels to retain for clustering= %d\n->Percentage of frames to use in each fold of consensus clustering= %d\n",Tmot,T,n_rep,Pp,Pn,Pcc);
    disp("*----------*");
    run=1;
end


%% START
%% ##Controls##
disp("Processing Controls..");
for i=1:group_controls:size(Data_OI_Controls,1)
   fprintf('\nProcessing: %i subjects.\n',size(Data_OI_Controls(i:1:i+group_controls-1),1));
   disp(Data_OI_Controls(i:1:i+group_controls-1));
   

   tmpTC = prepareData(Data_OI_Controls(i:1:i+group_controls-1),'step'); 
   
   if isMask==0
    mask=tmpTC.mask;
    SEED=setSeed(tmpTC,seed_num_type,seed_num); 
    isMask=1;
    brain_info=tmpTC.brain_info;
   end
    [XonTMP,~,IndicesTMP] = CAP_find_activity(tmpTC.TC,SEED,T,tmpTC.FD{1,1},Tmot,SelMode,SeedType,SignMatrix);
   
    Xon1=[Xon1,XonTMP];
    Indices1.scrubbed=[Indices1.scrubbed,IndicesTMP.scrubbed];
    Indices1.scrubbedandactive=[Indices1.scrubbedandactive,IndicesTMP.scrubbedandactive];
    Indices1.kept.active=[Indices1.kept.active,IndicesTMP.kept.active];
    
    %clear tmpTC;
    VOX=tmpTC.SubjSize;
    
    
end

%Patients will be loaded in a second step 
%% 4. Consensus clustering (if wished to determine the optimum K)
time1=datetime();
disp('Finding correct number K for clusters');
UnivrPlots={};


%se non viene fatto il consensus dovrebbe venire usato come K_opt

if (consensusk~=0)
    fprintf(1,'Checking max number of cluster = %d\n',K_range);
    % This specifies the range of values over which to perform consensus
    % clustering: if you want to run parallel consensus clustering processes,
    % you should feed in different ranges to each call of the function
    %K_range = 2:3;
    UnivrPlots.ConsensusPeformed='True';
    fprintf('## Consensus ##\n');
    t1=datetime();
    % Have each of these run in a separate process on the server =)
    for i = 2:K_range
        %[Consensus(:,:,i)] = %CAP_ConsensusClustering(Xon1,i+1,'items',Pcc/100,N,'correlation',maxiter);
        [Consensus(:,:,i)] = %CAP_ConsensusClustering(Xon1,i,'items',Pcc/100,N,'correlation',maxiter);
    end
    
   % parfor i=1:K_range
    %    [Consensus(:,:,i)] =CAP_ConsensusClustering(Xon1,i+1,'items',Pcc/100,N,'correlation',maxiter);
   % end
    t2=datetime();
    fprintf(1,'Consensus Clustering took: %s\n',t2-t1);
    clear t1;
    clear t2;
    %%%%%old:

    %K_range = 2:3;
    %[Consensus] = CAP_ConsensusClustering(Xon,K_range,'items',Pcc/100,N,'correlation');
    %%%%%%

    %%%%%%
    %
    % CHECK Qual after parallel modification
    %
    %%%%%%
    % Calculates the quality metrics
    [CDF,Qual] = ComputeClusteringQuality(Consensus,[]);
    %bar(2:K_range,1-Qual);
    % Qual should be inspected to determine the best cluster number(s)

    % You should fill this with the actual value :::: must inspect Qual
    % before... or find a method to automate this passage

    
        %% 4.1 Inspect Qual matrix
      %  bar(2:K_range,1-Qual);
           %% 4.1 Silouhette
        UnivrPlots.Silhouette.SilhouettePerformed='True';
        [eva,indices,sampledFrames,Silh]=compute_silhouette(Xon1,K_range,Pcc/100,N,maxiter);
        UnivrPlots.Silhouette.Silhouette=eva;
        %UnivrPlots.Silhouette.Silhouette2=eva2;
        UnivrPlots.Silhouette.Indices=indices;
        UnivrPlots.Silhouette.sampledData=sampledFrames;
        %K_opt = nan;
        %close all;

        UnivrPlots.Silhouette.S.X=Silh.X;
        UnivrPlots.Silhouette.S.idxs=Silh.clust;
        clear Silh;
        UnivrPlots.Consensus.Consensus=Consensus;
        UnivrPlots.Consensus.Quality=Qual;

else
    disp("Consensus clustering skipped");
    K_opt=K_range; %number of default cluster if you skipped this section
    UnivrPlots.Consensus.ConsensusPeformed='False';
    UnivrPlots.Silhoulette.SilhouettePerformed='False';
    
end
time2=datetime();
fprintf(1,'Total time: %s\n',time2-time1);
%Move this to end
%clear Qual;
%clear Consensus;
if (consensusk~=0)

    if is_screen==1
        %% 4.2 Inspecting both 
        close all;

        univrPlots(K_range,UnivrPlots,Qual,eva);

        UnivrPlots.K_range=K_range;

        fprintf(1,"--->Suggested K from silhouette (evalclusters) = %d",eva.OptimalK);

        %% 4.3 Ask number k to user
        %Ask user number K
        prompt="Select the number K of clusters you want: ";
        K_opt=input(prompt);
        clear custom_cm;
        clear tmp_plot;
        clear prompt;
        close all;
    
    else
        fprintf(1,"Not possible to show results, saving only the images.");
        univrPlots_no_show(K_range,UnivrPlots,Qual,eva);
        UnivrPlots.K_range=K_range;
        K_opt= default_k;
    
    end
    
end





UnivrSpatioTemporalSelection.NumberOfClusters=K_opt;
%% 5. Clustering into CAPs
fprintf("-> Clustering\n");
[CAP,~,~,idx1,CorrDist] = Run_Clustering(cell2mat(Xon1),K_opt,mask,brain_info{1},Pp,Pn,n_rep,[],SeedType,maxiter);
    
%[CAP,~,~,idx1,CorrDist] = Run_Clustering_UNIVR(cell2mat(Xon1),K_opt,n_rep,[],SeedType);

CAPs={};
CAPs.CAP_HEALTY=CAP;
CAPs.idx_HEALTY=idx1;

disp("Done");
%% 6. Assignment of the frames from population 2
%% ##Patients##
disp("Processing patients..");
for i=1:group_patients:size(Data_OI_Patients,1)
   fprintf('\nProcessing: %i subjects.\n',size(Data_OI_Patients(i:1:i+group_patients-1),1));
   disp(Data_OI_Patients(i:1:i+group_patients-1));
   

   tmpTC = prepareData(Data_OI_Patients(i:1:i+group_patients-1),'step'); 
   
   if isMask==0
    mask=tmpTC.mask;
    SEED=setSeed(tmpTC,seed_num_type,seed_num); 
    isMask=1;
    brain_info=tmpTC.brain_info;
   end
    [XonTMP,~,IndicesTMP] = CAP_find_activity(tmpTC.TC,SEED,T,tmpTC.FD{1,1},Tmot,SelMode,SeedType,SignMatrix);
   
    Xon2=[Xon2,XonTMP];
    Indices2.scrubbed=[Indices2.scrubbed,IndicesTMP.scrubbed];
    Indices2.scrubbedandactive=[Indices2.scrubbedandactive,IndicesTMP.scrubbedandactive];
    Indices2.kept.active=[Indices2.kept.active,IndicesTMP.kept.active];
    
    %clear tmpTC;
    VOX=tmpTC.SubjSize;
    
    
end
clear tmpTC;
clear isMask;
clear run;
clear XonTMP;


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

Indices1.kept.active=logical(Indices1.kept.active);
Indices1.scrubbedandactive=logical(Indices1.scrubbedandactive);
Indices1.scrubbed=logical(Indices1.scrubbed);

Indices2.kept.active=logical(Indices2.kept.active);
Indices2.scrubbedandactive=logical(Indices2.scrubbedandactive);
Indices2.scrubbed=logical(Indices2.scrubbed);


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
%UnivrData.FD_PATIENTS=FD2;
%UnivrData.FD_HEALTY=FD2;

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

%% 7.3 Save CAPs nifti
%Not present in original script (missing)
%adapting from CAP_TB.m main file 

fprintf('-->Now saving to nifti\n');

savename=datestr(datetime('today'));
fprintf(1,"\nOutput images will have names as %s\n",savename);
CAPToNIFTI(CAP,mask,brain_info{1},savedir,savename);

disp('->Completed converting to Nifti.');

%% 8. Save CAPs for every subject

savedir_subjects='/home/sbombieri/RESULTS/';
mkdir(savedir_subjects);
%getCAPS_forEachSubject();
%     -idx : uno tra idx1 o idx2
%     -DATASET: uno tra _PATIENTS o _HEALTY
%     -Data_OI: come DATASET
%     -savedir: directory where to save files
%     -Xon: Xon1 o Xon2

%PATIENTS
disp("SAVING PATIENTS CONTROLS CAPs....");
getCAPS_forEachSubject_twopop_MOD(idx2,VOX,Data_OI_Patients,savedir_subjects,Xon2,K_opt,"patients",CAP,mask,brain_info{1});

%HEALTY CONTROLS
disp("SAVING HEALTY CONTROLS CAPs....");
getCAPS_forEachSubject_twopop_MOD(idx1,VOX,Data_OI_Controls,savedir_subjects,Xon1,K_opt,"healty",CAP,mask,brain_info{1});



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
clear prompt;
clear ans;
clear IndicesTMP;
clear SEED;

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


