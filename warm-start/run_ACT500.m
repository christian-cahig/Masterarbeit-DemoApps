%% Initializations
SNAPSHOT_ID = 'ACT500';

% For generating dispatch data sets
NUM_INSTANCES = 490;
[PD_MUL_MIN, PD_MUL_MAX] = deal(0.9, 1.1);
[QD_MUL_MIN, QD_MUL_MAX] = deal(0.9, 1.1);
[RS_MUL_MIN, RS_MUL_MAX] = deal(0.9, 1.1);
[XS_MUL_MIN, XS_MUL_MAX] = deal(0.9, 1.1);
[BC_MUL_MIN, BC_MUL_MAX] = deal(0.9, 1.1);

% Go-to MATPOWER options
mpopt = getGoToMPOpts();

%% Snapshot data
load(sprintf("../snapshots/%s.mat", SNAPSHOT_ID));
BASEMVA = sdata.baseMVA;

% Deal with unbounded supply units
Pu_snap = sdata.gen(:, 2);
Qu_snap = sdata.gen(:, 3);

Pu_min = sdata.gen(:, 10);
nlbp_units = find(Pu_min == -Inf);
lbp_units = find(Pu_min ~= -Inf);
Pu_min(nlbp_units) = min(Pu_snap(nlbp_units), min(Pu_min(lbp_units)));
sdata.gen(:, 10) = Pu_min;
clear Pu_min nlbp_units lbp_units;

Pu_max = sdata.gen(:, 9);
nubp_units = find(Pu_max == Inf);
ubp_units = find(Pu_max ~= Inf);
Pu_max(nubp_units) = max(Pu_snap(nubp_units), max(Pu_max(ubp_units)));
sdata.gen(:, 9) = Pu_max;
clear Pu_max nubp_units ubp_units;

Qu_min = sdata.gen(:, 5);
nlbq_units = find(Qu_min == -Inf);
lbq_units = find(Qu_min ~= -Inf);
Qu_min(nlbq_units) = min(Qu_snap(nlbq_units), min(Qu_min(lbq_units)));
sdata.gen(:, 5) = Qu_min;
clear Qu_min nlbq_units lbq_units;

Qu_max = sdata.gen(:, 4);
nubq_units = find(Qu_max == Inf);
ubq_units = find(Qu_max ~= Inf);
Qu_max(nubq_units) = max(Qu_snap(nubq_units), max(Qu_max(ubq_units)));
sdata.gen(:, 4) = Qu_max;
clear Qu_max nubq_units ubq_units;

clear Pu_snap Qu_snap;

% Set supply costs
sdata.gencost(:, 5) = 1;
sdata.gencost(:, 6) = 1e-5;
sdata.gencost(:, 7) = 0;

%% For anticipating supply injections
SUPREG = 1.0;
cvx_solver Mosek;

%% For solving the APF equations
REFANG = 0;
EPSTOL = 1e-5;
NITERS = 100;
NEVALS = max(1000, 500 * size(sdata.bus, 1));
TRFALG = 'factorization';
LMDAMP = 1e-2;
LMDIAG = 'none';

%% Containers
areFeasible = zeros([NUM_INSTANCES, 1], 'logical');
snapObjs = zeros([NUM_INSTANCES, 1]);
uuObjs = zeros([NUM_INSTANCES, 1]);
uudiffCs = zeros([NUM_INSTANCES, 1]);
uudiffVms = zeros([NUM_INSTANCES, 1]);
uudiffVas = zeros([NUM_INSTANCES, 1]);
urObjs = zeros([NUM_INSTANCES, 1]);
urdiffCs = zeros([NUM_INSTANCES, 1]);
urdiffVms = zeros([NUM_INSTANCES, 1]);
urdiffVas = zeros([NUM_INSTANCES, 1]);
ruObjs = zeros([NUM_INSTANCES, 1]);
rudiffCs = zeros([NUM_INSTANCES, 1]);
rudiffVms = zeros([NUM_INSTANCES, 1]);
rudiffVas = zeros([NUM_INSTANCES, 1]);
rrObjs = zeros([NUM_INSTANCES, 1]);
rrdiffCs = zeros([NUM_INSTANCES, 1]);
rrdiffVms = zeros([NUM_INSTANCES, 1]);
rrdiffVas = zeros([NUM_INSTANCES, 1]);

%% Main loop
parfor i = 1:NUM_INSTANCES
    % Dispatch data
    ddata = makeDispData(sdata, MPOptions=mpopt, ...
        PdScaleRange=[PD_MUL_MIN; PD_MUL_MAX], ...
        QdScaleRange=[QD_MUL_MIN; QD_MUL_MAX], ...
        RScaleRange=[RS_MUL_MIN; RS_MUL_MAX], ...
        XScaleRange=[XS_MUL_MIN; XS_MUL_MAX], ...
        BScaleRange=[BC_MUL_MIN; BC_MUL_MAX], ...
        assertOPFOk=false);
    if ~ddata.success, continue; end
    areFeasible(i) = ddata.success;
    snapObjs(i) = ddata.f;
    snapIters(i) = ddata.raw.output.iterations;
    snapTimes(i) = ddata.et;
    [Cf, Ct] = makeCfCt(ddata);
    Cft = Cf - Ct;
    snapC = [ddata.var.val.Pg; ddata.var.val.Qg;];
    snapVm = ddata.var.val.Vm;
    snapVa = Cft * ddata.var.val.Va;

    % Problem data
    pdata = buildProbData(sdata, ddata, RefAngle=REFANG, calcCu=true, calcY=true);

    % APF essentials
    [Ph, Qh] = calcPhQh(pdata.dispatch, Vm=pdata.snapshot.bus(:, 8));
    [Po, Qo] = calcPoQo(pdata.dispatch, ...
        Vm=pdata.snapshot.bus(:, 8), ...
        Va=pdata.snapshot.bus(:, 9));
    Va_snap = deg2rad(pdata.snapshot.bus(:, 9));
    config = optimoptions("fsolve", ...
        Algorithm='trust-region-dogleg', ...
        SpecifyObjectiveGradient=true, ...
        FunctionTolerance=EPSTOL^2, ...
        SubproblemAlgorithm=TRFALG, ...
        InitDamping=LMDAMP, ...
        ScaleProblem=LMDIAG, ...
        FiniteDifferenceStepSize=sqrt(eps), ...
        MaxIterations=NITERS, ...
        MaxFunctionEvaluations=NEVALS, ...
        Display='off');

    % Warm-starting: UU
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=0, QuReg=0);
    InitPs = 0.5 * (sum(xed.Pu_min) - sum(xed.Pu) + sum(xed.Pu_max) - sum(xed.Pu));
    PsDistribs = xed.Pu_max - xed.Pu_min;
    PsDistribs = PsDistribs / sum(PsDistribs);
    config.Algorithm = 'trust-region-dogleg';
    vol = solveAPFEqns(pdata, config, Pu=xed.Pu, Qu=xed.Qu, ...
        InitVm=pdata.snapshot.bus(:, 8), InitVa=Va_snap, ...
        InitPs=InitPs, PsDistribs=PsDistribs);
    if ~vol.converged
        config.Algorithm = 'levenberg-marquardt';
        vol = solveAPFEqns(pdata, config, Pu=xed.Pu, Qu=xed.Qu, ...
            InitVm=pdata.snapshot.bus(:, 8), InitVa=Va_snap, ...
            InitPs=InitPs, PsDistribs=PsDistribs);
    end
    mpc = pdata.dispatch;
    mpc.bus(:, 8) = vol.Vm;
    mpc.bus(:, 9) = rad2deg(vol.Va);
    mpc.gen(:, 2) = pdata.sys_baseMVA * (xed.Pu + (vol.Ps * PsDistribs));
    mpc.gen(:, 3) = pdata.sys_baseMVA * xed.Qu;
    mpc = runopf(mpc, mpopt); mpc = ext2int(mpc);
    uuObjs(i) = mpc.f;
    uudiffCs(i) = norm([mpc.var.val.Pg; mpc.var.val.Qg] - snapC, Inf);
    uudiffVms(i) = norm(mpc.var.val.Vm - snapVm, Inf);
    uudiffVas(i) = norm((Cft * mpc.var.val.Va) - snapVa, Inf);

    % Warm-starting: UR
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=0, QuReg=SUPREG);
    InitPs = 0.5 * (sum(xed.Pu_min) - sum(xed.Pu) + sum(xed.Pu_max) - sum(xed.Pu));
    PsDistribs = xed.Pu_max - xed.Pu_min;
    PsDistribs = PsDistribs / sum(PsDistribs);
    config.Algorithm = 'trust-region-dogleg';
    vol = solveAPFEqns(pdata, config, Pu=xed.Pu, Qu=xed.Qu, ...
        InitVm=pdata.snapshot.bus(:, 8), InitVa=Va_snap, ...
        InitPs=InitPs, PsDistribs=PsDistribs);
    if ~vol.converged
        config.Algorithm = 'levenberg-marquardt';
        vol = solveAPFEqns(pdata, config, Pu=xed.Pu, Qu=xed.Qu, ...
            InitVm=pdata.snapshot.bus(:, 8), InitVa=Va_snap, ...
            InitPs=InitPs, PsDistribs=PsDistribs);
    end
    mpc = pdata.dispatch;
    mpc.bus(:, 8) = vol.Vm;
    mpc.bus(:, 9) = rad2deg(vol.Va);
    mpc.gen(:, 2) = pdata.sys_baseMVA * (xed.Pu + (vol.Ps * PsDistribs));
    mpc.gen(:, 3) = pdata.sys_baseMVA * xed.Qu;
    mpc = runopf(mpc, mpopt); mpc = ext2int(mpc);
    urObjs(i) = mpc.f;
    urdiffCs(i) = norm([mpc.var.val.Pg; mpc.var.val.Qg] - snapC, Inf);
    urdiffVms(i) = norm(mpc.var.val.Vm - snapVm, Inf);
    urdiffVas(i) = norm((Cft * mpc.var.val.Va) - snapVa, Inf);

    % Warm-starting: RU
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=SUPREG, QuReg=0);
    InitPs = 0.5 * (sum(xed.Pu_min) - sum(xed.Pu) + sum(xed.Pu_max) - sum(xed.Pu));
    PsDistribs = xed.Pu_max - xed.Pu_min;
    PsDistribs = PsDistribs / sum(PsDistribs);
    config.Algorithm = 'trust-region-dogleg';
    vol = solveAPFEqns(pdata, config, Pu=xed.Pu, Qu=xed.Qu, ...
        InitVm=pdata.snapshot.bus(:, 8), InitVa=Va_snap, ...
        InitPs=InitPs, PsDistribs=PsDistribs);
    if ~vol.converged
        config.Algorithm = 'levenberg-marquardt';
        vol = solveAPFEqns(pdata, config, Pu=xed.Pu, Qu=xed.Qu, ...
            InitVm=pdata.snapshot.bus(:, 8), InitVa=Va_snap, ...
            InitPs=InitPs, PsDistribs=PsDistribs);
    end
    mpc = pdata.dispatch;
    mpc.bus(:, 8) = vol.Vm;
    mpc.bus(:, 9) = rad2deg(vol.Va);
    mpc.gen(:, 2) = pdata.sys_baseMVA * (xed.Pu + (vol.Ps * PsDistribs));
    mpc.gen(:, 3) = pdata.sys_baseMVA * xed.Qu;
    mpc = runopf(mpc, mpopt); mpc = ext2int(mpc);
    ruObjs(i) = mpc.f;
    rudiffCs(i) = norm([mpc.var.val.Pg; mpc.var.val.Qg] - snapC, Inf);
    rudiffVms(i) = norm(mpc.var.val.Vm - snapVm, Inf);
    rudiffVas(i) = norm((Cft * mpc.var.val.Va) - snapVa, Inf);

    % Warm-starting: RR
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=SUPREG, QuReg=SUPREG);
    InitPs = 0.5 * (sum(xed.Pu_min) - sum(xed.Pu) + sum(xed.Pu_max) - sum(xed.Pu));
    PsDistribs = xed.Pu_max - xed.Pu_min;
    PsDistribs = PsDistribs / sum(PsDistribs);
    config.Algorithm = 'trust-region-dogleg';
    vol = solveAPFEqns(pdata, config, Pu=xed.Pu, Qu=xed.Qu, ...
        InitVm=pdata.snapshot.bus(:, 8), InitVa=Va_snap, ...
        InitPs=InitPs, PsDistribs=PsDistribs);
    if ~vol.converged
        config.Algorithm = 'levenberg-marquardt';
        vol = solveAPFEqns(pdata, config, Pu=xed.Pu, Qu=xed.Qu, ...
            InitVm=pdata.snapshot.bus(:, 8), InitVa=Va_snap, ...
            InitPs=InitPs, PsDistribs=PsDistribs);
    end
    mpc = pdata.dispatch;
    mpc.bus(:, 8) = vol.Vm;
    mpc.bus(:, 9) = rad2deg(vol.Va);
    mpc.gen(:, 2) = pdata.sys_baseMVA * (xed.Pu + (vol.Ps * PsDistribs));
    mpc.gen(:, 3) = pdata.sys_baseMVA * xed.Qu;
    mpc = runopf(mpc, mpopt); mpc = ext2int(mpc);
    rrObjs(i) = mpc.f;
    rrdiffCs(i) = norm([mpc.var.val.Pg; mpc.var.val.Qg] - snapC, Inf);
    rrdiffVms(i) = norm(mpc.var.val.Vm - snapVm, Inf);
    rrdiffVas(i) = norm((Cft * mpc.var.val.Va) - snapVa, Inf);
end
clear sdata ddata pdata snapC snapVm snapVa xed vol mpc;

%% Outro
fname = sprintf("%s_results", SNAPSHOT_ID);
save(sprintf("./%s.mat", fname), ...
    "SNAPSHOT_ID", "NUM_INSTANCES", ...
    "PD_MUL_MIN", "PD_MUL_MAX", ...
    "QD_MUL_MIN", "QD_MUL_MAX", ...
    "RS_MUL_MIN", "RS_MUL_MAX", ...
    "XS_MUL_MIN", "XS_MUL_MAX", ...
    "BC_MUL_MIN", "BC_MUL_MAX", ...
    "BASEMVA", ...
    "SUPREG", "REFANG", "EPSTOL", "NITERS", "NEVALS", ...
    "TRFALG", "LMDAMP", "LMDIAG", ...
    "areFeasible", ...
    "snapObjs", ...
    "uuObjs", "uudiffCs", "uudiffVms", "uudiffVas", ...
    "urObjs", "urdiffCs", "urdiffVms", "urdiffVas", ...
    "ruObjs", "rudiffCs", "rudiffVms", "rudiffVas", ...
    "rrObjs", "rrdiffCs", "rrdiffVms", "rrdiffVas");
fprintf("Saved to './%s.mat'\n", fname);
clear;
