%% VERSION 2.0
clear all;
clc;
savedir='/home/sbombieri/RESULTS';  
fprintf('##### Starting process #####\n');
% ###### VARIABLES ######
Data_OI=importdata('/home/sbombieri/DATA/TauAmyFDG/subjectlist.txt');
fprintf("Total subjects: %i \n",size(Data_OI,1));
%find best group (BETA)
%%d=divisors(size(Data_OI,1));
%group=d(:,2);

%######SCREEN INFO:
is_screen=0; %if screen is present set to 1

%###### max K
K_range=6; %max number of cluster you want to check 
prompt="Select max number of clusters u looking for: ";
K_range=input(prompt);

%##### CONSENSUS
prompt="Perform consensus? (0=NO, 1=YES)";
consensusk=0; %skip consensus if set to 0
consensusk=input(prompt);

%other way of divisors:
%fprintf("Divisors: \n");
%disp(d);
prompt="Select number of subjects for each group to load: ";
group=input(prompt);

clear d;
%group=5;
isMask=0; %0= to be calculated
run=0;

Xon=[];
Indices={};
Indices.scrubbed=[];
Indices.scrubbedandactive=[];
Indices.kept={};
Indices.kept.active=[];
UnivrMetrics={};
UnivrParameters={};
UnivrSpatioTemporalSelection={};
tt1=datetime();
%%  ###   START   ###

%
%per il momento scelto gruppi di 5 perchè ho una lista di 45 soggetti
for i=1:group:size(Data_OI,1)
   %disp(Data_OI(i:1:i+group-1));

   fprintf('\nProcessing: %i subjects.\n',size(Data_OI(i:1:i+group-1),1));
   disp(Data_OI(i:1:i+group-1));
   

   tmpTC = prepareData(Data_OI(i:1:i+group-1),'step'); 
   
   if isMask==0
    mask=tmpTC.mask;
    SEED=setSeed(tmpTC,1); 
    isMask=1;
    brain_info=tmpTC.brain_info;
   end

       

   
   
   
   
   
       %% 2. Specifying the main parameters

    % Threshold above which to select frames
    T = 30;

    % Selection mode ('Threshold' or 'Percentage')
    SelMode = 'Percentage';

    % Threshold of DATASET.FD{1} above which to scrub out the frame and also the t-1 and
    % t+1 frames (if you want another scrubbing setting, directly edit the
    % code)
    % ---> Is M in GUI
    Tmot = 0.1; %0.5 prima

    % Type of used seed information: select between 'Average','Union' or
    % 'Intersection'
    SeedType = 'Average';

    % Contains the information, for each seed (each row), about whether to
    % retain activation (1 0) or deactivation (0 1) time points

    %rilanciare cper vedere più kMax

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
    N = 50;


    disp('-> Parameters ok. ');
    %% 2.1 Recap parameters
    % show a recap of inserted parameters.
    if run==1
        disp("*----------*");
        fprintf(1,"Current parameters:\n->Tmot (M)=%f\n->T=%d\n->Number of repetitions of the K-means clustering algorithm= %d\n->Percentage of positive-valued voxels to retain for clustering= %d\n->Percentage of negative-valued voxels to retain for clustering= %d\n->Percentage of frames to use in each fold of consensus clustering= %d\n",Tmot,T,n_rep,Pp,Pn,Pcc);
        disp("*----------*");
        run=1;
    end
    %% 3. Selecting the frames to analyse    
    disp('->Selecting the frames to analyse');
    % Xon will contain the retained frames, and Indices will tag the time
    % points associated to these frames, for each subject (it contains a
    % subfield for retained frames and a subfield for scrubbed frames)
    [XonTMP,~,IndicesTMP] = CAP_find_activity(tmpTC.TC,SEED,T,tmpTC.FD{1,1},Tmot,SelMode,SeedType,SignMatrix);
    Xon=[Xon,XonTMP];
    Indices.scrubbed=[Indices.scrubbed,IndicesTMP.scrubbed];
    Indices.scrubbedandactive=[Indices.scrubbedandactive,IndicesTMP.scrubbedandactive];
    Indices.kept.active=[Indices.kept.active,IndicesTMP.kept.active];
    
    
    VOX=tmpTC.SubjSize;

end

clear IndicesTMP;
clear XonTMP;
clear tmpTC;
%% 4. Consensus clustering (if wished to determine the optimum K)
time1=datetime();
disp('Finding correct number K for clusters');
UnivrPlots={};


%se non viene fatto il consensus dovrebbe venire usato come K_opt
fprintf(1,'Checking max number of cluster = %d\n',K_range);
if (consensusk~=0)
    % This specifies the range of values over which to perform consensus
    % clustering: if you want to run parallel consensus clustering processes,
    % you should feed in different ranges to each call of the function
    %K_range = 2:3;
    UnivrPlots.ConsensusPeformed='True';
    fprintf('## Consensus ##\n');
    t1=datetime();
    % Have each of these run in a separate process on the server =)
    for i = 1:K_range-1
        [Consensus(:,:,i)] = CAP_ConsensusClustering(Xon,i+1,'items',Pcc/100,N,'correlation');
    end
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
    [eva,indices,sampledFrames,Silh]=compute_silhouette(Xon,K_range,Pcc/100,N);
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


    %% 4.2 Inspecting both 
    close all;
  
    univrPlots(K_range,UnivrPlots,Qual,eva);
%     %subplot(3,1,3);
%     %plot(eva2);
%     title("evalcusters(X_ss,'kmeans','silhouette','KList',1:6)");
%     ylim([0,1]);
%     xlim([2,K_range]);
%     subplot(2,1,2);
%     plot(eva);
%     title("evalclusters(X_ss,tClusters,'silhouette')");
%     ylim([0,1]);
%     xlim([2,K_range]);
%     subplot(2,1,1);
%     tmp_plot = bar(2:K_range,1-Qual);
%     xlabel('Cluster number K');
%     ylabel('Stability');
%     xlim([2-0.6,K_range+0.6]);
%     ylim([0,1]);
%     %set('Box','off');    
%    custom_cm = cbrewer('seq','Reds',25);     
%    colormap(custom_cm(6:25,:));
%    UnivrPlots.colormap=custom_cm;
    %UnivrPlots.figure=
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
    
end






UnivrSpatioTemporalSelection.NumberOfClusters=K_opt;
%% 5. Clustering into CAPs
disp('Clustering into CAPs...');

[CAP,~,~,idx] = Run_Clustering(cell2mat(Xon),K_opt,mask,brain_info{1},Pp,Pn,n_rep,[],SeedType);
%[CP2,Disp,Std_Clusters,idx,d,sfrac]
CAPs={};
CAPs.CAP=CAP;
CAPs.idx=idx;


%% 6. Computing metrics
disp('Computing metrics...');
% The TR of your data in seconds
TR = 3;  %default for now...
Indices.kept.active=logical(Indices.kept.active);
Indices.scrubbedandactive=logical(Indices.scrubbedandactive);
Indices.scrubbed=logical(Indices.scrubbed);
[ExpressionMap,Counts,Entries,Avg_Duration,Duration,TransitionProbabilities,...
    From_Baseline,To_Baseline,Baseline_resilience,Resilience,Betweenness,...
    InDegree,OutDegree,SubjectEntries] = Compute_Metrics_simpler(idx,...
    Indices.kept.active,Indices.scrubbedandactive,K_opt,TR);


%% 7. Save CAPs nifti
%Not present in original script (missing)
%adapting from CAP_TB.m main file 
disp('Done');
disp('Saving to nifti');

savename=datestr(datetime('today'));
fprintf(1,"Output images will have names as %s\n",savename);
CAPToNIFTI(CAP,mask,brain_info{1},savedir,savename);

disp('->Completed converting to Nifti.');

%% 8. Save CAPs for every subject


%run('getCAPS_forEachSubject.m');
savedir_subjects='/home/sbombieri/RESULTS/SUBJECTS/';
getCAPS_forEachSubject(idx,VOX,Data_OI,savedir_subjects,Xon,K_opt,CAP,mask,brain_info{1});


%% 6.1 Save metrics & other data
UnivrMetrics.ExpressionMap=ExpressionMap;
UnivrMetrics.Counts=Counts;
UnivrMetrics.Entries=Entries;
UnivrMetrics.Avg_Duration=Avg_Duration;
UnivrMetrics.Duration=Duration;
UnivrMetrics.TransitionProbabilities=TransitionProbabilities;
UnivrMetrics.From_Baseline=From_Baseline;
UnivrMetrics.To_Baseline=To_Baseline;
UnivrMetrics.Baseline_resilience=Baseline_resilience;
UnivrMetrics.Resilience=Resilience;
UnivrMetrics.Betweenness=Betweenness;
UnivrMetrics.InDegree=InDegree;
UnivrMetrics.OutDegree=OutDegree;
UnivrMetrics.SubjectEntries=SubjectEntries;
UnivrMetrics.TR=TR;

UnivrSpatioTemporalSelection.CAPs=CAPs;
%UnivrSpatioTemporalSelection.Xon=Xon;
UnivrSpatioTemporalSelection.Indices=Indices;


UnivrData.ListMaxNumberOfClusters=K_range;
UnivrSpatioTemporalSelection.Indices=Indices;
UnivrData.Subjects=Data_OI;

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
UnivrParameters.maxNumberOfClusters=K_range;
UnivrParameters.NumberOfClusters=K_opt;
%% 9
%Save data to current directory
fprintf('###### Saving data to RESULTS.mat #####\n');
Results=struct();
Results.Data=UnivrData;
Results.Metrics=UnivrMetrics;
Results.Plots=UnivrPlots;
Results.SpatioTemporalSelection=UnivrSpatioTemporalSelection;
Results.Parameters=UnivrParameters;
Results.SavedOnDate=datetime(now,'ConvertFrom','datenum');
save(fullfile(savedir,'RESULTS.mat'),'Results','-v7.3');
fprintf('Results will available at: %s\n',savedir);


tt2=datetime();
fprintf(1,'All process took: %s\n',tt2-tt1);

%% 9.1 Clear workspace
clear UnivrData;
clear UnivrPlots;
clear UnivrMetrics;
clear UnivrParameters;
clear UnivrSpatioTemporalSelection;

%clear Results;
clear savedir;
close all;
disp("Cleaning workspace");
clear CAP;
clear idx;
clear ExpressionMap;
clear Counts;
clear Entries;
clear Avg_Duration;
clear Duration;
clear TransitionProbabilities;
clear From_Baseline;
clear To_Baseline;
clear Baseline_resilience;
clear Resilience;
clear Betweenness;
clear InDegree;
clear OutDegree;
clear SubjectEntries;
clear K_opt;
clear SEED;
clear TR;
clear Tmot;
clear DATASET.mask;
clear DATASET.brain_info;

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
clear DATASET.FD{1};
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
clear DATASET.TC;
clear DATASET.ROI{1,1}Type;
clear DATASET;
clear ret;
clear OUT_SUBJECDATASET.TCAPS;
clear Qual;
clear eva;
clear Consensus;
clear sampledFrames;
clear time1;
clear time2;
clear tt1;
clear tt2;
clear SeedType;
clear savedir_subjects;

%disp("-> For plot consensus and silhouette type 'eval(UnivrPlots.figure);' <-");
disp('-> Figure of silhouette & consensus clustering is saved in the current directory, named "Max_clust.fig" <-');

fprintf(1,"\n######## All ended correctly ########\n");
