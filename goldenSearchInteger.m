function out = goldenSearchInteger(fun, lb, ub, options)

if ~exist('options','var')
    tol = 1;
    extendBoundary = [false, false];
    ftol = 1e-3;
    verbose = true;
else
    tol = options.tol;
    extendBoundary = options.extendBoundary;
    ftol = options.ftol;
    verbose = options.verbose;

end

evalpoints = zeros(4,1);
invGr = 2/(1+sqrt(5)); %inverse of golden ratio
evalpoints(1) = lb;
evalpoints(4) = ub;
evalpoints(2) = lb + floor((1-invGr)*(ub - lb));
evalpoints(3) = lb + floor(invGr*(ub - lb));

val = zeros(4,1);
for i = 1:length(evalpoints)
    val(i) = fun(evalpoints(i));
end
[~,minId] = min(val);
if minId == 4 
    if extendBoundary(2) == true
        while minId == 4
            evalpoints(2:3) = evalpoints(3:4);
            val(2:3) = val(3:4);
            ub = ceil(ub/invGr);
            evalpoints(4) = ub;
            val(4) = fun(evalpoints(4));
            [~,minId] = min(val);
        end
        if verbose == true
            fprintf("%s\n","Extending upper bound");
            fprintf("%s %d\n","New upper bound:",ub);
        end
    end
end
if minId == 1
    if extendBoundary(1) == true
        while minId == 1
            evalpoints(2:3) = evalpoints(1:2);
            val(2:3) = val(1:2);
            lb = floor(lb*invGr);
            evalpoints(1) = lb;
            val(1) = fun(evalpoints(1));
            [~,minId] = min(val);
        end
        if verbose == true
            fprintf("%s\n","Extending lower bound");
            fprintf("%s %d\n","New lower bound:",lb);
        end
    end
end

iter = 0;
while (ub - lb) > 3*tol %while the lb and ub are at least 3 integers away
    iter = iter + 1;
    if verbose == true
        fprintf("%s %d\n","Iteration:",iter);
    end
    if minId == 1 || minId == 2
        ub = evalpoints(3);
        evalpoints(4) = ub;
        val(4) = val(3);
        evalpoints(3) = evalpoints(2);
        val(3) = val(2);
        evalpoints(2) = lb + floor((1-invGr)*(ub - lb));
        val(2) = fun(evalpoints(2));
        if verbose == true
            fprintf("ub = %d, val at ub = %.3f\n",ub, val(4));
            fprintf("evalpoints(3) = %d, val(3) = %.3f\n",evalpoints(3),val(3));
            fprintf("evalpoints(2) = %d, val(2) = %.3f\n",evalpoints(2),val(2));
        end
        if abs(val(3)-val(2) )< ftol
            break
        end
        
    elseif minId == 3 || minId == 4
        lb = evalpoints(2);
        evalpoints(1) = lb;
        val(1) = val(2);
        evalpoints(2) = evalpoints(3);
        val(2) = val(3);
        evalpoints(3) = lb + floor(invGr*(ub - lb));
        val(3) = fun(evalpoints(3));
        if verbose == true
            fprintf("lb = %d, val at lb = %.3f\n",lb, val(1));
            fprintf("evalpoints(2) = %d, val(2) = %.3f\n",evalpoints(2),val(2));
            fprintf("evalpoints(3) = %d, val(3) = %.3f\n",evalpoints(3),val(3));
        end
        if abs(val(3)-val(2)) < ftol
            break
        end

    end
    [minVal,minId] = min(val);
    if verbose == true
        fprintf("Current minId = %d\n",minId);
        fprintf("Current min val = %.3f\n",minVal);
    end
end
[minVal,minId] = min(val);
sol = evalpoints(minId);
if verbose == true
    fprintf("%s\n","Optimum Found");
    fprintf("Optimum solution = %d\n",sol);
    fprintf("Objective value = %.3f\n",minVal);
end
out.sol = sol;
out.val = val(minId);
return
end


