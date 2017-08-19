

function [PCA_Properties]= PCA_analysis(Event_Properties,Network_Properties);

%% PCA on event properties

%Import data
event_dur=Event_Properties.noNaN.duration;
MAX_PKS=Event_Properties.noNaN.peak;
event_amp=Event_Properties.noNaN.amplitude;
event_mean=Event_Properties.noNaN.mean;
event_AUC=Event_Properties.noNaN.AUC;
event_width=Event_Properties.noNaN.width;

%PCA for each event 
for u=1:size(MAX_PKS,2); for uu=1:size(MAX_PKS{u},2); 
allproperties{u}(uu,:)=[event_dur{u}(uu)...
MAX_PKS{u}(uu)...
event_amp{u}(uu)...
event_mean{u}(uu)...
event_AUC{u}(uu)...
event_width{u}(uu)];  
end
end
allproperties_mat=cell2mat(allproperties');

% Do the PCA for each event
[coeff_e,score_e,latent_e,~,explained_e] = pca(allproperties_mat);
% Calculate eigenvalues and eigenvectors of the covariance matrix
covarianceMatrix_e = cov(allproperties_mat);
[V_e,D_e] = eig(covarianceMatrix_e);
% "coeff" are the principal component vectors. These are the eigenvectors of the covariance matrix. 
% Multiply the original data by the principal component vectors to get the projections of the original data on the
% principal component vector space. This is also the output "score". 
dataInPrincipalComponentSpace_e = allproperties_mat*coeff_e;
%Plot first 2 PCs
vbls = {'dur','peak','ampl','meandf','AUC','width'};
%figure; biplot(coeff_e(:,1:2),'scores',score_e(:,1:2),'varlabels',vbls);

PCA_Properties.Event.dataPCS=dataInPrincipalComponentSpace_e;
PCA_Properties.Event.score=score_e;
PCA_Properties.Event.coeff=coeff_e;
PCA_Properties.Event.latent=latent_e;
PCA_Properties.Event.explained=explained_e;
PCA_Properties.Event.covarianceMatrix=covarianceMatrix_e;
PCA_Properties.Event.varlabels=vbls;
PCA_Properties.Event.allvariables=allproperties_mat;

%% PCA for each cell: event + network properties
%Import data
freq=Network_Properties.frequency;
AUC_rate=Network_Properties.AUC_rate;
corr_pairs=Network_Properties.nb_corr_pairs;


for u=1:size(MAX_PKS,2); 
allproperties_cells(u,:)=[mean(event_dur{u})...
mean(MAX_PKS{u})...
mean(event_amp{u})...
mean(event_mean{u})...
mean(event_AUC{u})...
mean(event_width{u})...
freq(u)...
AUC_rate(u)...
corr_pairs(u)]; 
end

% Do the PCA for each cell
[coeff_c,score_c,latent_c,~,explained_c] = pca(allproperties_cells);
%Weighted PCA
%[coeff_c,score_c,latent_c,~,explained_c] = pca(allproperties_cells, 'VariableWeights','variance');
% Calculate eigenvalues and eigenvectors of the covariance matrix
covarianceMatrix_c = cov(allproperties_cells);
[V_c,D_c] = eig(covarianceMatrix_c);
% "coeff" are the principal component vectors. These are the eigenvectors of the covariance matrix. 
% Multiply the original data by the principal component vectors to get the projections of the original data on the
% principal component vector space. This is also the output "score". 
dataInPrincipalComponentSpace_c = allproperties_cells*coeff_c;
%Plot first 2 PCs
vbls = {'dur','peak','ampl','meandf','AUC','width','freq','AUC rate','corr_pairs'};
%figure; biplot(coeff_c(:,1:2),'scores',score_c(:,1:2),'varlabels',vbls);


PCA_Properties.Cell.dataPCS=dataInPrincipalComponentSpace_c;
PCA_Properties.Cell.score=score_c;
PCA_Properties.Cell.coeff=coeff_c;
PCA_Properties.Cell.latent=latent_c;
PCA_Properties.Cell.explained=explained_c;
PCA_Properties.Cell.covarianceMatrix=covarianceMatrix_c;
PCA_Properties.Cell.varlabels=vbls;
PCA_Properties.Cell.allvariables=allproperties_cells;
end