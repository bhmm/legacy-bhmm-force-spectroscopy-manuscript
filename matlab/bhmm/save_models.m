function save_models(filename, models)
% Save critical model data (without state trajectories) to specified file.
%
% save_models(filename)
%
% ARGUMENTS
%  filename - file to save models to

% Determine number of models.
nmodels = length(models);

% Create storage for compressed models.
packed_models = [];

for model = 1:nmodels
end

