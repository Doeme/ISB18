function fid=write_escher_file(model, v, fname)
fid=fopen(fname, 'wt');
rxns=model.rxns;
% Tweak reaction names (fixes most but not all discrepancies)
rxns=strrep(model.rxns,'(e)','_e');
rxns=strrep(rxns,'-D','_D');
rxns=strrep(rxns,'-L','_L');
fprintf(fid,'ID,stage\n');
if fid == -1, return, end
for i=1:length(model.rxns)
    fprintf(fid,'%s,%1.8f\n', rxns{i}, v(i));
end
fclose(fid);
