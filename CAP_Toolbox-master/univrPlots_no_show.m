function [h]= univrPlots_no_show(K_range,UnivrPlots,Qual,eva)
    
    %FOLDER to save
    %savedir='/home/sbombieri/RESULTS/';
    %% EVALCLUSTERS & CONSENSUS
    ghostfigure1 = figure('Name','Consensus & Evalclusters with silhouette | 1',"Visible",false);
    title("evalcusters(X_ss,'kmeans','silhouette','KList',1:6)");
    ylim([0,1]);
    xlim([2,K_range]);
    ghostplot1 = plot(eva);
    ghostfigure2 = figure('Name','Consensus & Evalclusters with silhouette | 2',"Visible",false);
    title("evalclusters(X_ss,tClusters,'silhouette')");
    ylim([0,1]);
    xlim([2,K_range]);
    ghostplot2  = bar(Qual);
    xlabel('Cluster number K');
    ylabel('Stability');
    xlim([2-0.6,K_range+0.6]);
    ylim([0,1]);
    saveas(ghostfigure1,'Max_clust1_CONSENSUS.png');
    saveas(ghostfigure2,'Max_clust2_EVALCLUSTERS.png');


    %% SILHOUETTE
    ghostfigure3 = figure('Name','Silhouette',"Visible",false);
    L=max(UnivrPlots.Silhouette.Silhouette.InspectedK);
    ghostsil = silhouette(UnivrPlots.Silhouette.S.X,UnivrPlots.Silhouette.S.idxs,'Euclidean');
    title(sprintf('Run with K-max=%d', K_range));
    xlim([0,1]);
    %sgtitle('Silhouettes');
    savefig('Max_clust3_SILHOUETTE.fig');
    
    
    
end


