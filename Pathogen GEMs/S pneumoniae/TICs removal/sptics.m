%% remove TICs from the SP model
load('SPD39_host.mat')
[SPD39_host, removedRxnInd2, keptRxnInd2] = checkDuplicateRxn(SPD39_host,'FR');
[SPD39_host, removedRxnInd2, keptRxnInd2] = checkDuplicateRxn(SPD39_host,'S');
fba=optimizeCbModel(SPD39_host)
[rxnInLoopfinal, Nfinal, loopInfofinal] = findMinNull(SPD39_host);
find(rxnInLoopfinal(:,1)==1);
backward1a=SPD39_host.rxns(ans);
find(rxnInLoopfinal(:,2)==1);
forward1a=SPD39_host.rxns(ans);
both1a=intersect(forward1a,backward1a);
%%
block={'GLUDy'
'KARI'
'KARI_3hmoa'
'LAORP'
'PCOP'
'ribnp'
'ACt6'
'ADEt2'
'NDPK8'};
blockid = findRxnIDs(SPD39_host,block(1:8));
SPD39_host.lb(blockid)=0;
SPD39_host.ub(blockid)=0;
fba1=optimizeCbModel(SPD39_host)
%%
[rxnInLoopfinal, Nfinal, loopInfofinal] = findMinNull(SPD39_host);
find(rxnInLoopfinal(:,1)==1);
backward1a=SPD39_host.rxns(ans);
find(rxnInLoopfinal(:,2)==1);
forward1a=SPD39_host.rxns(ans);
both1a=intersect(forward1a,backward1a);
%%
blockid = findRxnIDs(SPD39_host,backward1a);
SPD39_host.lb(blockid)=0;
fba2=optimizeCbModel(SPD39_host)
%%
forward={'PYK'
'pykd'};
forid = findRxnIDs(SPD39_host,forward);
SPD39_host.ub(forid)=0;
fba3=optimizeCbModel(SPD39_host)