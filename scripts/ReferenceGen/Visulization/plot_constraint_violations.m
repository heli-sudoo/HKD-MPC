function plot_constraint_violations(vl)
% vl -> strcture array of constraint violations
% vl.econstr -> equality constraint violation
% vl.pconstr -> path constraint violation
% vl.tconstr -> terminal constraint violation

if ~isempty(vl(1).econstr)
    figure
    plot(0:length(vl)-1, [vl(:).econstr]);
    xlabel("outer loop iteration");
    ylabel('equality constraint');
end
if ~isempty(vl(1).pconstr)
    figure
    plot(0:length(vl)-1, [vl(:).pconstr]);
    xlabel("outer loop iteration");
    ylabel('path constraint');
end
if ~isempty(vl(1).tconstr)
    plot(0:length(vl)-1, [vl(:).tconstr]);
    xlabel("outer loop iteration");
    ylabel('terminal constraint');
end
end