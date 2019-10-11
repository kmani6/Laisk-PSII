function AnalysisNames(tspan,analyses_name)

    for i = 1:length(analyses_name)
        [Sol] = LaiskControlFluorescence(tspan,analyses_name{i});
    end

end 