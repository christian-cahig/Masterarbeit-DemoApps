%% Initializations
SNAPSHOT_ID = 'RTE1888';

% For generating dispatch data sets 
[PD_MUL_MIN, PD_MUL_MAX] = deal(0.9, 1.1);
[QD_MUL_MIN, QD_MUL_MAX] = deal(0.9, 1.1);
[RS_MUL_MIN, RS_MUL_MAX] = deal(0.95, 1.05);
[XS_MUL_MIN, XS_MUL_MAX] = deal(0.95, 1.05);
[BC_MUL_MIN, BC_MUL_MAX] = deal(0.95, 1.05);

% For anticipating supply injections
REGS = [0; 1e-4; 1];
[PUREGS, QUREGS] = ndgrid(REGS, REGS);
PUREGS = PUREGS(:);
QUREGS = QUREGS(:);
NUM_INSTANCES = numel(PUREGS);

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
load(sprintf("../snapshots/%s.mat", SNAPSHOT_ID));

%% Dispatch data
ddata = makeDispData(sdata, MPOptions=mpopt, ...
    PdScaleRange=[PD_MUL_MIN; PD_MUL_MAX], ...
    QdScaleRange=[QD_MUL_MIN; QD_MUL_MAX], ...
    RScaleRange=[RS_MUL_MIN; RS_MUL_MAX], ...
    XScaleRange=[XS_MUL_MIN; XS_MUL_MAX], ...
    BScaleRange=[BC_MUL_MIN; BC_MUL_MAX], ...
    assertOPFOk=false);
while ~ddata.success
    ddata = makeDispData(sdata, MPOptions=mpopt, ...
        PdScaleRange=[PD_MUL_MIN; PD_MUL_MAX], ...
        QdScaleRange=[QD_MUL_MIN; QD_MUL_MAX], ...
        RScaleRange=[RS_MUL_MIN; RS_MUL_MAX], ...
        XScaleRange=[XS_MUL_MIN; XS_MUL_MAX], ...
        BScaleRange=[BC_MUL_MIN; BC_MUL_MAX], ...
        assertOPFOk=false);
end

%% Problem data
pdata = buildProbData(sdata, ddata, RefAngle=0);

SYSMVA = pdata.sys_baseMVA;
Pu_snap = sdata.gen(:, 2) ./ SYSMVA;
Qu_snap = sdata.gen(:, 3) ./ SYSMVA;

Vm_snap = pdata.snapshot.bus(:, 8);
Va_snap = pdata.snapshot.bus(:, 9);
[Ph, Qh] = calcPhQh(pdata.dispatch, Vm=Vm_snap);
[Po, Qo] = calcPoQo(pdata.dispatch, Vm=Vm_snap, Va=Va_snap);
Va_snap = deg2rad(Va_snap);

Vm_min = pdata.dispatch.bus(:, 13);
Vm_max = pdata.dispatch.bus(:, 12);

%% Solver configuration and initialization
NEVALS = max(1000, 500 * size(sdata.bus, 1));
xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, useCQP=false);
PsDistribs = xed.Pu_max - xed.Pu_min;
PsDistribs = PsDistribs / sum(PsDistribs);
Pu_min = xed.Pu_min;
Pu_max = xed.Pu_max;
Qu_min = xed.Qu_min;
Qu_max = xed.Qu_max;

%% Containers
Pu = zeros([NUM_INSTANCES, pdata.N_u]);
Qu = zeros([NUM_INSTANCES, pdata.N_u]);
apfOks = zeros([NUM_INSTANCES, 1], 'logical');
Vm = zeros([NUM_INSTANCES, pdata.N_b]);
Va = zeros([NUM_INSTANCES, pdata.N_b]);
Ps_min = zeros([NUM_INSTANCES, 1]);
Ps_max = zeros([NUM_INSTANCES, 1]);
Ps = zeros([NUM_INSTANCES, 1]);

%% Main loop
cvx_solver Mosek;
parfor i = 1:NUM_INSTANCES
    preg = PUREGS(i);
    qreg = QUREGS(i);
    config = optimoptions("fsolve", ...
        SpecifyObjectiveGradient=true, ...
        FunctionTolerance=EPSTOL^2, ...
        SubproblemAlgorithm=TRFALG, ...
        InitDamping=LMDAMP, ...
        ScaleProblem=LMDIAG, ...
        FiniteDifferenceStepSize=sqrt(eps), ...
        MaxIterations=NITERS, ...
        MaxFunctionEvaluations=NEVALS, ...
        Display='off');

    % Anticipate supply injections
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=preg, QuReg=qreg);
    Pu(i, :) = xed.Pu;
    Qu(i, :) = xed.Qu;
    Ps_min(i) = sum(xed.Pu_min) - sum(xed.Pu);
    Ps_max(i) = sum(xed.Pu_max) - sum(xed.Pu);
    InitPs = 0.5 * (Ps_min(i) + Ps_max(i));

    % Solve APF equations
    config.Algorithm = 'trust-region-dogleg';
    vol = solveAPFEqns(pdata, config, Pu=xed.Pu, Qu=xed.Qu, PsDistribs=PsDistribs, ...
        InitVm=Vm_snap, InitVa=Va_snap, InitPs=InitPs);
    if ~vol.converged
        config.Algorithm = 'levenberg-marquardt';
        vol = solveAPFEqns(pdata, config, Pu=xed.Pu, Qu=xed.Qu, PsDistribs=PsDistribs, ...
            InitVm=Vm_snap, InitVa=Va_snap, InitPs=InitPs);
    end
    apfOks(i) = vol.converged;
    Vm(i, :) = vol.Vm;
    Va(i, :) = vol.Va;
    Ps(i) = vol.Ps;
end
clear sdata ddata pdata;

%% Outro
fname = sprintf("%s_results", SNAPSHOT_ID);
save(sprintf("./%s.mat", fname), ...
    "SNAPSHOT_ID", "REGS", "PUREGS", "QUREGS", "NUM_INSTANCES", ...
    "PD_MUL_MIN", "PD_MUL_MAX", ...
    "QD_MUL_MIN", "QD_MUL_MAX", ...
    "RS_MUL_MIN", "RS_MUL_MAX", ...
    "XS_MUL_MIN", "XS_MUL_MAX", ...
    "BC_MUL_MIN", "BC_MUL_MAX", ...
    "REFANG", "EPSTOL", "NITERS", "NEVALS", ...
    "TRFALG", "LMDAMP", "LMDIAG", ...
    "Pu_snap", "Pu_min", "Pu_max", "Pu", ...
    "Qu_snap", "Qu_min", "Qu_max", "Qu", ...
    "apfOks", ...
    "Vm_snap", "Vm_min", "Vm_max", "Vm", ...
    "Va_snap", "Va", ...
    "PsDistribs", "Ps_min", "Ps_max", "Ps");
fprintf("Saved to './%s.mat'\n", fname);
clear;
