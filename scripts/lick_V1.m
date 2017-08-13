%Lick response
Behavior_all=Behavior;
Behavior=Behavior_all{4};


position=Behavior.position;
lick=Behavior.lick;


lick_binary=diff(lick)==1; 
lick_binary=[lick_binary; 0];

goal=10 %20 cm
goalbelt=(max(position)-goal)
goal_ind=find(position>=goalbelt);
lickgoal=lick_binary(goal_ind);


nb_lickgoal=sum(lickgoal);
nb_licktot=sum(lick_binary);
fraction_lickgoal=nb_lickgoal/nb_licktot

tic;
save(strcat(directory_name,'\behavior'), '-v7.3');
toc;


