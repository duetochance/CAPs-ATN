function [h]= univrPlots(K_range,UnivrPlots,Qual,eva)
    close all;
    
    %% EVALCLUSTERS & CONSENSUS
    h(1)=figure('Name','Consensus & Evalclusters with silhouette');
    movegui('east')
    title("evalcusters(X_ss,'kmeans','silhouette','KList',1:6)");
    ylim([0,1]);
    xlim([2,K_range]);
    subplot(2,1,2);
    plot(eva);
    title("evalclusters(X_ss,tClusters,'silhouette')");
    ylim([0,1]);
    xlim([2,K_range]);
    subplot(2,1,1);
    %tmp_plot = bar(2:K_range,1-Qual);
    bar(Qual);
    xlabel('Cluster number K');
    ylabel('Stability');
    xlim([2-0.6,K_range+0.6]);
    ylim([0,1]);
    %set('Box','off');


    %% SILHOUETTE
    h(2)=figure('Name','Silhouette');
    movegui('west')
    L=max(UnivrPlots.Silhouette.Silhouette.InspectedK);
    %for i=1:L
    %subplot(L,1,i);
    %silhouette(UnivrPlots.Silhouette.Silhouette.X,UnivrPlots.Silhouette.Indices(:,i),'Euclidean');
    silhouette(UnivrPlots.Silhouette.S.X,UnivrPlots.Silhouette.S.idxs,'Euclidean');
    %silhouette(UnivrPlots.Silhouette.Silhouette.X,UnivrPlots.Silhouette.Indices,'Euclidean');
    title(sprintf('Run with K-max=%d', K_range));
    xlim([0,1]);
    %end
    sgtitle('Silhouettes');
    savefig(h,'Max_clust');
    
    
    
end


