
function [Imaging]=baselinesub(Imaging, options);

%% Import 
windwith=options.windwith;
C_df{1}=Imaging.trace;
C_df{2}=Imaging.trace_restricted;


%% Yout = msbackadj(X, Intensities)
for t=1:2 %full trace and restricted 
for i=1:size(C_df{t},2);
x{t}=(1:length(C_df{t}))';
C_df_sub{t}(:,i) = msbackadj(x{t}, C_df{t}(:,i),'StepSize',windwith);
end
end

%% Display figure
if options.dispfig==true,
c2plot=options.c2plot;

C_df_sub_c2plot{1} = msbackadj(x{1}, C_df{1}(:,c2plot),'StepSize', windwith, 'ShowPlot', true);
figure; 
subplot(2,1,1) 
plot(C_df{1}(:,c2plot));
subplot(2,1,2) 
plot(C_df_sub_c2plot{1});

end

%% Make structure
Imaging.trace_baselinesub=C_df_sub{1};
Imaging.trace_restricted_baselinesub=C_df_sub{2};

Imaging.options.msbackadj=options.msbackadj;
Imaging.options.windwith=options.windwith;

end
