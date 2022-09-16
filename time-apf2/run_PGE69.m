%% Initializations
SNAPSHOT = 'PGE69';

% For generating dispatch data sets 
NPROBS = 1260;
[PDXMIN, PDXMAX] = deal(0.9, 1.1);
[QDXMIN, QDXMAX] = deal(0.9, 1.1);
[RSXMIN, RSXMAX] = deal(0.95, 1.05);
[XSXMIN, XSXMAX] = deal(0.95, 1.05);
[BCXMIN, BCXMAX] = deal(0.95, 1.05);

% For anticipating supply injections
PU_REG = 1e-3;
QU_REG = 1e-3;
cvx_solver Mosek;

% For solving the APF equations
REFANG = 0;
EPSTOL = 1e-5;
NITERS = 100;
TRFALG = 'factorization';
LMDAMP = 1e-2;
LMDIAG = 'none';

% Go-to MATPOWER options
mpopt = getGoToMPOpts();

%% Snapshot data
load(sprintf("../snapshots/%s.mat", SNAPSHOT));

%% Solver configurations
NEVALS = max(1000, 500 * size(sdata.bus, 1));

%% Containers
areFeasible = zeros([NPROBS, 1], "logical");
phOks = zeros([NPROBS, 1], "logical");
phEvals = zeros([NPROBS, 1]);
phSteps = zeros([NPROBS, 1]);
phTimes = zeros([NPROBS, 1]);
phNorms = zeros([NPROBS, 1]);
lmOks = zeros([NPROBS, 1], "logical");
lmEvals = zeros([NPROBS, 1]);
lmSteps = zeros([NPROBS, 1]);
lmTimes = zeros([NPROBS, 1]);
lmNorms = zeros([NPROBS, 1]);

%% Warmup loop
parfor i = 1:17
    CONFIG = optimoptions("fsolve", ...
        SpecifyObjectiveGradient=true, ...
        FunctionTolerance=EPSTOL^2, ...
        SubproblemAlgorithm=TRFALG, ...
        InitDamping=LMDAMP, ...
        ScaleProblem=LMDIAG, ...
        FiniteDifferenceStepSize=sqrt(eps), ...
        MaxIterations=NITERS, ...
        MaxFunctionEvaluations=NEVALS, ...
        Display='off');

    % Dispatch data
    ddata = makeDispData(sdata, MPOptions=mpopt, ...
        PdScaleRange=[PDXMIN; PDXMAX], ...
        QdScaleRange=[QDXMIN; QDXMAX], ...
        RScaleRange=[RSXMIN; RSXMAX], ...
        XScaleRange=[XSXMIN; XSXMAX], ...
        BScaleRange=[BCXMIN; BCXMAX], ...
        assertOPFOk=false);
    if ~ddata.success, continue; end

    % Problem data
    pdata = buildProbData(sdata, ddata, RefAngle=REFANG, calcCu=true, calcCd=true, calcY=true);
    Vm_snap = pdata.snapshot.bus(:, 8);
    Va_snap = pdata.snapshot.bus(:, 9);
    [Ph, Qh] = calcPhQh(pdata.dispatch, Vm=Vm_snap);
    [Po, Qo] = calcPoQo(pdata.dispatch, Vm=Vm_snap, Va=Va_snap);
    Va_snap = deg2rad(Va_snap);

    % Anticipate supply injections
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=PU_REG, QuReg=QU_REG);
    PsDistribs = xed.Pu_max - xed.Pu_min;
    PsDistribs = PsDistribs / sum(PsDistribs);
    MinPs = sum(xed.Pu_min) - sum(xed.Pu);
    MaxPs = sum(xed.Pu_max) - sum(xed.Pu);
    InitPs = 0.5*(MinPs + MaxPs);

    % Solve APF equations via Powell hybrid method
    CONFIG.Algorithm = 'trust-region-dogleg';
    vol = solveAPFEqns(pdata, CONFIG, Pu=xed.Pu, Qu=xed.Qu, PsDistribs=PsDistribs, ...
        InitVm=Vm_snap, InitVa=Va_snap, InitPs=InitPs);

    % Solve APF equations via Levenberg-Marquardt algorithm
    CONFIG.Algorithm = 'levenberg-marquardt';
    vol = solveAPFEqns(pdata, CONFIG, Pu=xed.Pu, Qu=xed.Qu, PsDistribs=PsDistribs, ...
        InitVm=Vm_snap, InitVa=Va_snap, InitPs=InitPs);
end

%% Main loop
parfor i = 1:NPROBS
    CONFIG = optimoptions("fsolve", ...
        SpecifyObjectiveGradient=true, ...
        FunctionTolerance=EPSTOL^2, ...
        SubproblemAlgorithm=TRFALG, ...
        InitDamping=LMDAMP, ...
        ScaleProblem=LMDIAG, ...
        FiniteDifferenceStepSize=sqrt(eps), ...
        MaxIterations=NITERS, ...
        MaxFunctionEvaluations=NEVALS, ...
        Display='off');

    % Dispatch data
    ddata = makeDispData(sdata, "MPOptions", mpopt, ...
        "PdScaleRange", [PDXMIN; PDXMAX], ...
        "QdScaleRange", [QDXMIN; QDXMAX], ...
        "RScaleRange", [RSXMIN; RSXMAX], ...
        "XScaleRange", [XSXMIN; XSXMAX], ...
        "BScaleRange", [BCXMIN; BCXMAX], ...
        "assertOPFOk", false);
    if ddata.success, areFeasible(i) = true; else, continue; end

    % Problem data
    pdata = buildProbData(sdata, ddata, RefAngle=REFANG, calcCu=true, calcCd=true, calcY=true);
    Vm_snap = pdata.snapshot.bus(:, 8);
    Va_snap = pdata.snapshot.bus(:, 9);
    [Ph, Qh] = calcPhQh(pdata.dispatch, Vm=Vm_snap);
    [Po, Qo] = calcPoQo(pdata.dispatch, Vm=Vm_snap, Va=Va_snap);
    Va_snap = deg2rad(Va_snap);

    % Anticipate supply injections
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=PU_REG, QuReg=QU_REG);
    PsDistribs = xed.Pu_max - xed.Pu_min;
    PsDistribs = PsDistribs / sum(PsDistribs);
    MinPs = sum(xed.Pu_min) - sum(xed.Pu);
    MaxPs = sum(xed.Pu_max) - sum(xed.Pu);
    InitPs = 0.5*(MinPs + MaxPs);

    % Solve APF equation via Powell hybrid method
    CONFIG.Algorithm = 'trust-region-dogleg';
    vol = solveAPFEqns(pdata, CONFIG, Pu=xed.Pu, Qu=xed.Qu, PsDistribs=PsDistribs, ...
        InitVm=Vm_snap, InitVa=Va_snap, InitPs=InitPs);
    phOks(i) = vol.converged;
    phEvals(i) = vol.FuncEvals;
    phSteps(i) = vol.Steps;
    phTimes(i) = vol.Time;
    phNorms(i) = norm(vol.Phi, Inf);

    % Solve APF equations via Levenberg-Marquardt algorithm
    CONFIG.Algorithm = 'levenberg-marquardt';
    vol = solveAPFEqns(pdata, CONFIG, Pu=xed.Pu, Qu=xed.Qu, PsDistribs=PsDistribs, ...
        InitVm=Vm_snap, InitVa=Va_snap, InitPs=InitPs);
    lmOks(i) = vol.converged;
    lmEvals(i) = vol.FuncEvals;
    lmSteps(i) = vol.Steps;
    lmTimes(i) = vol.Time;
    lmNorms(i) = norm(vol.Phi, Inf);
end
clear sdata ddata pdata Vm_snap Va_snap xed vol;

%% Outro
fname = sprintf("%s_results-%f", SNAPSHOT, sum(areFeasible) / NPROBS);
save(sprintf("./%s.mat", fname), ...
    "SNAPSHOT", "NPROBS", ...
    "PDXMIN", "PDXMAX", ...
    "QDXMIN", "QDXMAX", ...
    "RSXMIN", "RSXMAX", ...
    "XSXMIN", "XSXMAX", ...
    "BCXMIN", "BCXMAX", ...
    "PU_REG", "QU_REG", ...
    "REFANG", "EPSTOL", "NITERS", "NEVALS", ...
    "TRFALG", "LMDAMP", "LMDIAG", ...
    "areFeasible", ...
    "phOks", "phEvals", "phSteps", "phTimes", "phNorms", ...
    "lmOks", "lmEvals", "lmSteps", "lmTimes", "lmNorms");
fprintf("Saved to './%s.mat'\n", fname);
clear;
