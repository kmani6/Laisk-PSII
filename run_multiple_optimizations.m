function run_multiple_optimizations(analysis_name, nsamples, nparallel,randomseed)


if nargin <4
    rng('shuffle')
else
    rng(randomseed)
end
seeds = rand(nsamples,1);
parfor (i = 1:nsamples,nparallel)
    get_optimal_parameters(analysis_name,seeds(i))
end

end