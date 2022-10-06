%% Init 
% SUBJECTCAPS is a variable similar to CAP of toolbox but instead of
% containing general CAPs (i.e. CAP1, CAP2 ...) it contains the CAP for
% each subject,so each row represents the CAPn (where n is number of CAP)
% for that specific subject
% 

%% Args
%     -idx : uno tra idx1 o idx2
%     -DATASET: uno tra _PATIENTS o _HEALTY
%     -Data_OI: come DATASET
%     -savedir: directory where to save files
%     -Xon: Xon1 o Xon2
%     -K_opt
%     -type: "healty" or "patients"

function [] = getCAPS_forEachSubject_twopop_MOD(idx,DATASET,Data_OI,savedir,Xon,K_opt,type,CAP,mask,brain_info)
warning('off');
disp(CAP);
    if type=="healty"
        %DATASET.SubjSize.VOX=DATASET_SubjSize.VOX;
        savedirCLASS=fullfile(savedir,"HEALTY");
        fprintf(1,"-> Creating HEALTY folder %s\n",savedirCLASS);
        mkdir(savedirCLASS)
    else
        %DATASET.SubjSize.VOX=DATASET_SubjSize.VOX;
        savedirCLASS=fullfile(savedir,"PATIENTS");
        fprintf(1,"-> Creating PATIENTS folder %s\n",savedirCLASS);
        mkdir(savedirCLASS)
    end
    SUBJECTCAPS=zeros(K_opt,DATASET.VOX);
    OUT_SUBJECTCAPS=struct;
    %for each subject take the clusters present
    %%%
    % Xon is named HeavyOutputs.SpatioTemporalSelection.ClusteredFrames by
    % toolbox
    %%%
    tmp_idx=idx;
    for sub=1:size(Data_OI,1)
        disp('#####');
        fprintf(1,'Subject no. %d\n',sub);
        sub_tmp=Xon{1,sub};
        idx_sub=tmp_idx(1:size(sub_tmp,2));
        tmp_idx(1:size(sub_tmp,2))=[];
        fprintf(1,'sub_tmp : %d x %d\n',size(sub_tmp));
        fprintf(1,'idx_sub : %d x %d\n',size(idx_sub));
        disp('#####');
        %done with assignments
        for i=1:K_opt
            cl=idx_sub==i; 
            xontmp=sub_tmp(:,cl);
            meanxontmp=mean(xontmp,2);
            m=meanxontmp.';
            SUBJECTCAPS(i,:)=m;


        end
        disp('--------------------');
        %save for workspace
        varname=genvarname(append('Subj_',num2str(sub)));
        OUT_SUBJECTCAPS.varname=SUBJECTCAPS;
        %save NIFTI
        tosavedir=fullfile(savedirCLASS,varname);
        disp(tosavedir);
        mkdir(tosavedir);
        fprintf(1,"-> savedir= %s\n",tosavedir);
        fprintf(1,"-> varname= %s\n",varname);
         
        CAPToNIFTI(CAP,mask,brain_info,tosavedir,varname);



    end
    disp('Completed and saved to Results/');










    %% clear some now useleess variables

    clear tmp_idx;
    warning('on');
end
