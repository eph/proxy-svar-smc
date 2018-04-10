newData1 = importdata(data_file);
% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
YYdata = newData1.data;
text = newData1.textdata;
clear data textdata
nDate = datenum(text(2:end,1));
[~,i_var] = ismember(i_var_str,text(1,2:end));
[~,i_instr] = ismember(i_var_instr,text(1,2:end));


%***********************************************/
% RETRIEVE POSITIONS OF VARIABLES IN THE DATASET/
%***********************************************/


[~,i_transf] = ismember(i_var_transf,i_var_str);

%************************************************/
% RETRIEVE POSITION OF FIRST AND LAST OBSERVATION/
%************************************************/

sample_init = datenum(str_sample_init, 'yyyy-mm-dd');
sample_end = datenum(str_sample_end, 'yyyy-mm-dd');
sample_iv_init = datenum(str_iv_init, 'yyyy-mm-dd');

[~, sample_init_row] = ismember(sample_init,nDate,'rows');
[~, sample_end_row] = ismember(sample_end,nDate,'rows');
[~, sample_iv_row] = ismember(sample_iv_init,nDate,'rows');

%****************************************************/
% SELECT APPROPRIATE ROWS AND COLUMNS OF DATA MATRIX /
%****************************************************/
YY = YYdata(sample_init_row:sample_end_row,i_var);
mm = YYdata(sample_iv_row:sample_end_row,i_instr);
