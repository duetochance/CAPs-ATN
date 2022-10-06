function [FD] = univr_CAP_ComputeFD(motfile_name)
    %fprintf(1,'*FD from file\n --> %s',motfile_name);
    Mot = textread(motfile_name);
    Mot = Mot(:,1:6);

    % Converts the rotational components into [mm]
    Mot(:,1:3) = 50*Mot(:,1:3); % prima (:,4:6) rotazioni-traslazioni in fsl

    % Computes FD
    FD = sum(abs([0 0 0 0 0 0; diff(Mot)]),2);

end