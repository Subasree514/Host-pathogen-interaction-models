%%
%load('Streptococcus pneumoniae D39.mat');
%commonrsmets_sp={'cu2[c]';'fe2[c]';'fe3[c]';'gthrd[c]';'h2o2[c]';'h2o2[e]';'h2s[c]';'mthgxl[c]';'no3[c]';'so3[c]';'so4[c]';'trdox[c]';'trdrd[c]'};
%rsrxns1=findRxnsFromMets(model,{'o2s[c]','o2s[e]','h2o2[c]','h2o2[e]','no[c]','no[e]','gthox[c]','gthox[e]','gthrd[c]','gthrd[e]','h2s[c]','h2s[e]','no2[c]','no3[c]'});
%rsrxns2=findRxnsFromMets(model,commonrsmets_sp);
%%
%load('Streptococcus_pneumoniae_G54.mat')
%commonrsmets_sp={'cu2[c]';'fe2[c]';'fe3[c]';'gthrd[c]';'h2o2[c]';'h2o2[e]';'h2s[c]';'mthgxl[c]';'no3[c]';'so3[c]';'so4[c]';'trdox[c]';'trdrd[c]'};
%rsrxns3=findRxnsFromMets(model,{'o2s[c]','o2s[e]','h2o2[c]','h2o2[e]','no[c]','no[e]','gthox[c]','gthox[e]','gthrd[c]','gthrd[e]','h2s[c]','h2s[e]','no2[c]','no3[c]'});
%rsrxns4=findRxnsFromMets(model,commonrsmets_sp);
%%
clear
load('iIT341_modified.mat')
model2=iIT341;
dmem_tr={'cys__L[c] + h2o[c] -> h2s[c] + nh4[c] + pyr[c]'
'acser[c] + trdrd[c] + tsul[c] <=> ac[c] + cys__L[c] + h[c] + so3[c] + trdox[c]'
'dhor__S[c] + o2[c] -> h2o2[c] + orot[c]'
'atp[c] + fe3[e] + h2o[c] -> adp[c] + fe3[c] + h[c] + pi[c] '
'4 fe2[c] + 4 h[c] + o2[c] -> 4 fe3[c] + 2 h2o[c]'
'cgly[c] + glu__L[c] -> gthrd[c] + h2o[c]'
'h2o[c] + hcys__L[c] <=> 2obut[c] + h2s[c] + nh4[c]'
'h[e] + so4[e] <=> h[c] + so4[c]'
'h2o2[c] + trdrd[c] -> 2 h2o[c] + trdox[c]'};
rxnNamesnew={'CYSDS'
'CYSS3r'
'DHORDi'
'FE3abc'
'FERO'
'GGTA'
'HCYSHSL'
'SO4t2'
'THIORDXi'};
dmemsubsytem={'Methionine and cysteine metabolism'
'Methionine and cysteine metabolism'
'Pyrimidine synthesis'
'Transport, extracellular'
'Vitamin B2 metabolism'
'Glutathione metabolism'
'Sulfur metabolism'
'Transport, extracellular'
'ROS detoxification'};
for iIT341 = 1:length(rxnNamesnew)
model2 = addReaction(model2,rxnNamesnew{iIT341,1},'reactionFormula',dmem_tr{iIT341,1},'subSystem',dmemsubsytem{1,1});
end
%%
model2 = addReaction(model2,'h2o2_tr','reactionFormula','h2o2[e] <=> h2o2[c]');
model2 = addReaction(model2,'h2s_tr','reactionFormula','h2s[e] <=> h2s[c]');
model2 = addReaction(model2,'o2s_tr','reactionFormula','o2s[e] <=> o2s[c]');
model2=addExchangeRxn(model2,{'o2s[e]';'h2s[e]';'h2o2[e]'});
%%
rsexchanges={'EX_o2s[e]';'EX_h2s[e]';'EX_h2o2[e]'};
rsids=findRxnIDs(model2,rsexchanges);
model2.lb(rsids)=0;
iIT341=model2;
save('iIT341_modified','iIT341')