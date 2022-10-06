
function [eva,tClusters,X_ss,toSilhouette] = compute_silhouette(X,maxK,Subsample_fraction,N,maxiter)
%compute_silhouette 
%   Compute silhouette
disp('## Silhouette ##');
% Number of subjects (from CAP_ConsensusClustering function from the
% toolbox)
n_subjects = length(X);

n_items = 0;

frames_index = cell(1,n_subjects);

for s = 1:n_subjects
    frames_index{s} = n_items + (1:size(X{s},2));
    n_items = n_items + size(X{s},2);
end

% Number of dimensions
n_dims = size(X{1},1);

%% Setting up things
%% Applying silhouette

% Number of items to subsample
n_items_ss = floor(Subsample_fraction*n_items); %qui viene fuori il 20, 0.8*25=20 , n_items sono i frame ritenuti dopo la selezione iniziale
%disp([n_items_ss,Subsample_fraction,n_items]);
% Does the subsampling

[X_ss,~] = datasample((cell2mat(X))',n_items_ss,1,'Replace',false);

sizesampled=size(X_ss,1);
tClusters=zeros(sizesampled,maxK);

t1=datetime();




% ############## EVALCLUSTERS (Silhouette) ##########
% If criterion is 'CalinskiHarabasz', 'DaviesBouldin', or 'silhouette', you
% can also specify clust as a n-by-K matrix containing the proposed
% clustering solutions. n is the number of observations in the sample data,
% and K is the number of proposed clustering solutions. Column j contains
% the cluster indices for each of the N points in the jth clustering
% solution.
%eva = evalclusters(x,clust,criterion)
%eva2=evalclusters(X_ss2,'kmeans','silhouette','klist',1:maxK);
for i=1:maxK
    tClusters(:,i) = kmeans(X_ss,i,'Replicates',N,'Start','uniform','maxiter',maxiter);
end

eva=evalclusters(X_ss,tClusters,'silhouette');
t2=datetime();

fprintf(1,'Evalclusters kmeans time: %s\n',t2-t1);

%########## SILHOUETTE  ##########
clust=kmeans(X_ss,maxK,'Replicates',N,'Start','uniform','maxiter',maxiter);
toSilhouette={};
toSilhouette.X=X_ss;
toSilhouette.clust=clust;
%% Tests
%[CP2,Disp,Std_Clusters,idx,d,sfrac] = Run_Clustering(Xsampled,maxK,mask,brain_info{1},Pp,Pn,n_rep,[],SeedType);
end



