load('model_final.mat');
load('rxnTisMatccle2023.mat')
run ccle2023info.m;
model=model_final;
objectiveAbbr = printObjective(model);
%%
tissue_t={'ACH-000681';'ACH-000758'};
%%
tissueids_t=contains(ccleconditions,tissue_t);
c_id_t=find(tissueids_t==1);
%%
epsilon=getCobraSolverParams('LP', 'feasTol')*100;
%%
c1_t=find(rxnTisMatccle(:,c_id_t(1))==1);
c2_t=find(rxnTisMatccle(:,c_id_t(2))==1);
%% swiftcore
order=1:length(model.rxns);
all=order';
smat=model.S;
revs=find(model.lb<0);
irrevs=setdiff(all,revs);
run penaltyweights.m;
full=order';
full(revs)=1;
full(irrevs)=0;
%consistent = swiftcc(model.S,full,'gurobi');
%swiftcoreconsistent=intersect(consistent,rxncore);
result = blocked(model.S,full,'gurobi');
modelnewS.S=model.S;
modelnewS.c=model.c;
modelnewS.lb=model.lb;
modelnewS.ub=model.ub;
modelnewS.rxns=model.rxns;
modelnewS.rev=full;
component = partition(modelnewS,'gurobi', 'swift');
new=result.x(length(model.mets):(length(model.mets)+length(model.rxns)))
blockednumbers=find(new==-1);
newfull(blockednumbers)=1;
unblockednumbers=setdiff(all,blockednumbers);
newfull(unblockednumbers)=0;
flux = core(modelnewS,newfull,weights,'gurobi');
find(flux==0);
epsilon=getCobraSolverParams('LP', 'feasTol')*100;
%%
corematrix = {c1_t;c2_t};
for i = 1:length(corematrix)
corem1=[731;773;954;2271;2277;7488;5190;corematrix{i,1}];
[modelswiftcore(i), reconInd, LP] = swiftcore(model,corem1,weights,epsilon,false,'gurobi');
end
%%
modelb1=modelswiftcore(1)
modelb1.lb=modelb1.lb.*1e3;
modelb1.ub=modelb1.ub.*1e3;
%model1=modelb1;
%dmemnormalconstraints
%dmemtransportnormal
%modelout=model2;
%[model1Out, removedRxnInd, keptRxnInd] = checkDuplicateRxn(modelout,'S');
modelA549final2=modelb1;
%%
save('modelA549final2','modelA549final2');
%%
modelb1=modelswiftcore(2)
modelb1.lb=modelb1.lb.*1e3;
modelb1.ub=modelb1.ub.*1e3;
%model1=modelb1;
%dmemnormalconstraints
%dmemtransportnormal
%modelout=model2;
%[model1Out, removedRxnInd, keptRxnInd] = checkDuplicateRxn(modelb1,'S');
modelHK74final2=modelb1;
%%
save('modelHK28final2','modelHK74final2');