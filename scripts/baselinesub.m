
function [Imaging]=baselinesub(Imaging);

windwith=Imaging.options.baselinesub.windwith;
c2plot=Imaging.options.celltoplot;
C_df{1}=Imaging.trace;
C_df{2}=Imaging.trace_restricted;

for t=1:2 %full trace and restricted 
for i=1:size(C_df{t},2);
x{t}=(1:length(C_df{t}))';
C_df_sub{t}(:,i) = msbackadj(x{t}, C_df{t}(:,i),'StepSize',windwith);
end
end
C_df_sub_c2plot{1} = msbackadj(x{1}, C_df{1}(:,c2plot),'StepSize', windwith, 'ShowPlot', true);
figure; plot(C_df_sub_c2plot{1});
figure; plot(C_df{1}(:,c2plot));
Imaging.trace_baselinesub=C_df_sub{1};
Imaging.trace_restricted_baselinesub=C_df_sub{2};
end
