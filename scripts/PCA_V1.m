


% Do the PCA for traces
[coeff_c,score_c,latent_c,~,explained_c] = pca(C_df);
%Weighted PCA
%[coeff_c,score_c,latent_c,~,explained_c] = pca(allproperties_cells, 'VariableWeights','variance');

% Calculate eigenvalues and eigenvectors of the covariance matrix
covarianceMatrix_c = cov(C_df);
[V_c,D_c] = eig(covarianceMatrix_c);
% "coeff" are the principal component vectors. These are the eigenvectors of the covariance matrix. 
% Multiply the original data by the principal component vectors to get the projections of the original data on the
% principal component vector space. This is also the output "score". 
dataInPrincipalComponentSpace_c = C_df*coeff_c;
%Plot first 2 PCs
%figure; biplot(coeff_c(:,1:2),'scores',score_c(:,1:2));

time_cdf=(1:size(C_df,1))*(1/30.2061910850898);
time_behavior=(1:size(texture,1))*(1/10000);

figure;
for i=1:5
subplot(1,5,i) 
 plot(time_cdf, (score_c(:,i)));
 hold on; plot(time_behavior, (tex_binary(:,:)));
end

 