

function [hist_events]= histo_events(Event_Properties,options);

%% Histogram
%Import data
%all analyzed
event_dur=Event_Properties.noNaN.duration;
MAX_PKS=Event_Properties.noNaN.peak;
event_amp=Event_Properties.noNaN.amplitude;
event_mean=Event_Properties.noNaN.mean;
event_AUC=Event_Properties.noNaN.AUC;
event_width=Event_Properties.noNaN.width;

%Kake table
for u=1:size(MAX_PKS,2); for uu=1:size(MAX_PKS{u},2); 
allproperties{u}(uu,:)=[event_dur{u}(uu)...
MAX_PKS{u}(uu)...
event_amp{u}(uu)...
event_mean{u}(uu)...
event_AUC{u}(uu)...
event_width{u}(uu)];  
end
end

if options.exclude==true
%excluded
event_dur_exc=Event_Properties.Excluded.event_dur;
MAX_PKS_exc=Event_Properties.Excluded.MAX_PKS;
event_amp_exc=Event_Properties.Excluded.event_amp;
event_mean_exc=Event_Properties.Excluded.event_mean;
event_AUC_exc=Event_Properties.Excluded.event_AUC;
event_width_exc=Event_Properties.Excluded.event_width;
for u=1:size(MAX_PKS_exc,2); for uu=1:size(MAX_PKS_exc{u},2); 
allproperties_exc{u}(uu,:)=[event_dur_exc{u}(uu)...
MAX_PKS_exc{u}(uu)...
event_amp_exc{u}(uu)...
event_mean_exc{u}(uu)...
event_AUC_exc{u}(uu)...
event_width_exc{u}(uu)];  
end
end
allproperties_hist_exc=cell2mat(allproperties_exc');
allproperties_hist_exc=allproperties_hist_exc(~any(isnan(allproperties_hist_exc),2),:);
end
allproperties_hist=cell2mat(allproperties');
allproperties_hist=allproperties_hist(~any(isnan(allproperties_hist),2),:);

figure
nb_para=size(allproperties_hist,2);
legendinfo = {'event duration','event peak','event amplitude','event mean dF/F','event AUC','event width'};
for i=1:nb_para
subplot(nb_para,1,i)
histogram(allproperties_hist(:,i))
hold on;
if options.exclude==true
histogram(allproperties_hist_exc(:,i))
legend('analyzed events', 'excluded events')
end
title(legendinfo{i})
ylabel('Number of events');
end


%% PCA
%On analyzed events
[coeff_e,score_e,latent_e,~,explained_e] = pca(allproperties_hist);
% Calculate eigenvalues and eigenvectors of the covariance matrix
covarianceMatrix_e = cov(allproperties_hist);
[V_e,D_e] = eig(covarianceMatrix_e);
% "coeff" are the principal component vectors. These are the eigenvectors of the covariance matrix. 
% Multiply the original data by the principal component vectors to get the projections of the original data on the
% principal component vector space. This is also the output "score". 
dataInPrincipalComponentSpace_e = allproperties_hist*coeff_e;
%Plot first 2 PCs
vbls = {'dur','peak','ampl','meandf','AUC','width'};
%figure; biplot(coeff_e(:,1:2),'scores',score_e(:,1:2),'varlabels',vbls);

if options.exclude==true
%On excluded events
[coeff_e_exc,score_e_exc,latent_e_exc,~,explained_e_exc] = pca(allproperties_hist_exc);
% Calculate eigenvalues and eigenvectors of the covariance matrix
covarianceMatrix_e_exc = cov(allproperties_hist_exc);
[V_e_exc,D_e_exc] = eig(covarianceMatrix_e_exc);
% "coeff" are the principal component vectors. These are the eigenvectors of the covariance matrix. 
% Multiply the original data by the principal component vectors to get the projections of the original data on the
% principal component vector space. This is also the output "score". 
dataInPrincipalComponentSpace_e_exc = allproperties_hist_exc*coeff_e_exc;
%Plot first 2 PCs
%figure; biplot(coeff_e(:,1:2),'scores',score_e(:,1:2),'varlabels',vbls);
end

figure;
scatter3(score_e(:,1),score_e(:,2),score_e(:,3));
hold on;
if options.exclude==true
scatter3(score_e_exc(:,1),score_e_exc(:,2), score_e_exc(:,3));
legend('analyzed events', 'excluded events')
end
title('PCA')



hist_events.PCA.dataPCS=dataInPrincipalComponentSpace_e;
hist_events.PCA.score=score_e;
hist_events.PCA.coeff=coeff_e;
hist_events.PCA.latent=latent_e;
hist_events.PCA.explained=explained_e;
hist_events.PCA.covarianceMatrix=covarianceMatrix_e;
hist_events.PCA.varlabels=vbls;
hist_events.allvariables=allproperties_hist;
if options.exclude==true
hist_events.Excluded.PCA.dataPCS=dataInPrincipalComponentSpace_e_exc;
hist_events.Excluded.PCA.score=score_e_exc;
hist_events.Excluded.PCA.coeff=coeff_e_exc;
hist_events.Excluded.PCA.latent=latent_e_exc;
hist_events.Excluded.PCA.explained=explained_e_exc;
hist_events.Excluded.PCA.covarianceMatrix=covarianceMatrix_e_exc;
hist_events.Excluded.PCA.varlabels=vbls;
hist_events.Excluded.allvariables=allproperties_hist_exc;
end

end