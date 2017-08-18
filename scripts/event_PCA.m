

%% PCA on event properties

Event_Properties.duration=event_dur;
Event_Properties.peak=MAX_PKS;
Event_Properties.amplitude=event_amp;
Event_Properties.mean=event_mean;
Event_Properties.AUC=event_AUC;
Event_Properties.width=event_width;

%PCA for each cell with mean values for each event
for u=1:size(peakidx,2); 
 allproperties_cells(u,:)=[EventsProperties.mean.duration(u)...
 EventsProperties.mean.peak(u)...
 EventsProperties.mean.amplitude(u)...
 EventsProperties.mean.meandf(u)...
EventsProperties.mean.area(u)...
 EventsProperties.mean.risetime(u)...       
 EventsProperties.mean.decaytime(u)...
 EventsProperties.mean.halfwidth(u)];     
end; 
allproperties_cells=allproperties_cells(~any(isnan(allproperties_cells),2),:);

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
vbls = {'dur','peak','ampl','meandf','area','rise','decay','halfwidth'};
figure; biplot(coeff_c(:,1:2),'scores',score_c(:,1:2),'varlabels',vbls);

%PCA for each event
for u=1:size(peakidx,2); for uu=1:size(peakidx{u},2); 
 allproperties{u}(uu,:)=[EventsProperties.duration{u}(uu)...
 EventsProperties.peak{u}(uu)...
 EventsProperties.amplitude{u}(uu)...
 EventsProperties.meandf{u}(uu)...
 EventsProperties.area{u}(uu)...
 EventsProperties.risetime{u}(uu)...       
 EventsProperties.decaytime{u}(uu)...
 EventsProperties.halfwidth{u}(uu)];     
end; end;
allproperties_mat=cell2mat(allproperties');
allproperties_mat=allproperties_mat(~any(isnan(allproperties_mat),2),:);
EventsProperties.nbrevents.analyzed=size(allproperties_mat,1);


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
vbls = {'dur','peak','ampl','meandf','area','rise','decay','halfwidth'};
figure; biplot(coeff_e(:,1:2),'scores',score_e(:,1:2),'varlabels',vbls);

PCAProperties.Event.dataPCS=dataInPrincipalComponentSpace_e;
PCAProperties.Event.score=score_e;
PCAProperties.Event.coeff=coeff_e;
PCAProperties.Event.latent=latent_e;
PCAProperties.Event.explained=explained_e;
PCAProperties.Event.covarianceMatrix=covarianceMatrix_e;
PCAProperties.Event.varlabels=vbls;
PCAProperties.Event.allvariables=allproperties_mat;


PCAProperties.Cell.dataPCS=dataInPrincipalComponentSpace_c;
PCAProperties.Cell.score=score_c;
PCAProperties.Cell.coeff=coeff_c;
PCAProperties.Cell.latent=latent_c;
PCAProperties.Cell.explained=explained_c;
PCAProperties.Cell.covarianceMatrix=covarianceMatrix_c;
PCAProperties.Cell.varlabels=vbls;
PCAProperties.Cell.allvariables=allproperties_cells;