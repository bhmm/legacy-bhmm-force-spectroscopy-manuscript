% Generate LaTeX table.
nlengths = NL;

string = '\begin{tabular}{|c|c|';
for length_index = 1:NL
  string = sprintf('%sc', string);
end
string = sprintf('%s|}',string);
disp(string)
disp('\hline');

% header line
disp(sprintf('property & true value & \\multicolumn{%d}{c|}{observation length} \\\\', nlengths));
string = ' &';
for length_index = 1:NL
  string = sprintf('%s & $10^%d$', string, round(log10(L_n(length_index))));
end
string = sprintf('%s \\\\ \\hline', string);
disp(string)

% \pi_i
property = 'Pi';
for i = 1:nstates
  string = sprintf('$\\pi_{%d}$ & $%.3f$', i, true_model.Pi(i));
  for length_index = 1:NL
    string = sprintf('%s & $%.3f$', string, model_store{length_index}.Pi(i));
  end
  string = sprintf('%s \\\\', string);
  disp(string)
end  

disp('\hline');

% T_ij
property = 'Tij';
for i = 1:nstates
  for j = 1:nstates
    string = sprintf('$T_{%d%d}$ & $%.3f$', i, j, true_model.Tij(i,j));
    for length_index = 1:NL
      string = sprintf('%s & $%.3f$', string, model_store{length_index}.Tij(i,j));
    end
    string = sprintf('%s \\\\', string);
    disp(string)
  end
end  

disp('\hline');

% mu_i
property = 'mu_i';
for i = 1:nstates
  string = sprintf('$\\mu_{%d}$ & $%.3f$', i, true_model.states{i}.mu);
  for length_index = 1:NL
    string = sprintf('%s & $%.3f$', string, model_store{length_index}.states{i}.mu);
  end
  string = sprintf('%s \\\\', string);
  disp(string)
end  

disp('\hline');

% sigma_i
property = 'sigma_i';
for i = 1:nstates
  string = sprintf('$\\sigma_{%d}$ & $%.3f$', i, true_model.states{i}.sigma);
  for length_index = 1:NL
    string = sprintf('%s & $%.3f$', string, model_store{length_index}.states{i}.sigma);
  end
  string = sprintf('%s \\\\', string);
  disp(string)
end  

disp('\hline');

disp('\end{tabular}');
