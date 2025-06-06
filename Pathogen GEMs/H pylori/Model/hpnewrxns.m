load('iIT341_modified.mat')
model2=iIT341;
%% add the following reactions to the iIT341 model
rxns_new={'cys__L[c] + h2o[c] -> h2s[c] + nh4[c] + pyr[c]'
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
rxnssubsystem={'Methionine and cysteine metabolism'
'Methionine and cysteine metabolism'
'Pyrimidine synthesis'
'Transport, extracellular'
'Vitamin B2 metabolism'
'Glutathione metabolism'
'Sulfur metabolism'
'Transport, extracellular'
'ROS detoxification'};

for i = 1:length(rxnNamesnew)
model2 = addReaction(model2,rxnNamesnew{i,1},'reactionFormula',rxns_new{i,1},'subSystem',rxnssubsystem{i,1});
end

%% update the following metabolites and add their exchanges
model2=addExchangeRxn(model2,{'o2s[e]';'h2s[e]';'h2o2[e]'});
rsexchanges={'EX_o2s[e]';'EX_h2s[e]';'EX_h2o2[e]'};
rsids=findRxnIDs(model2,rsexchanges);
model2.lb(rsids)=0;

%% add transport reactions for the three metabolites
model2 = addReaction(model2,'h2o2_tr','reactionFormula','h2o2[e] <=> h2o2[c]');
model2 = addReaction(model2,'h2s_tr','reactionFormula','h2s[e] <=> h2s[c]');
model2 = addReaction(model2,'o2s_tr','reactionFormula','o2s[e] <=> o2s[c]');

iIT341=model2;
save('iIT341_modified','iIT341')