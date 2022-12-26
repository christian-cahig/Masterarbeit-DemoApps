%% Initializations
SNAPSHOT_ID = 'AEP57';

% For generating dispatch data sets
[PD_MUL_MIN, PD_MUL_MAX] = deal(0.95, 1.05);
[QD_MUL_MIN, QD_MUL_MAX] = deal(0.95, 1.05);
[RS_MUL_MIN, RS_MUL_MAX] = deal(0.99, 1.01);
[XS_MUL_MIN, XS_MUL_MAX] = deal(0.99, 1.01);
[BC_MUL_MIN, BC_MUL_MAX] = deal(0.99, 1.01);
MUL_STD_DEV = 0.05;
NUM_SAMPLES = 1260;

%% Go-to MATPOWER options
mpopt = getGoToMPOpts();

%% Snapshot data
load(sprintf("../snapshots/%s.mat", SNAPSHOT_ID));

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

%% Base dispatch data
bdata = makeDispData(sdata, MPOptions=mpopt, ...
    PdScaleRange=[PD_MUL_MIN; PD_MUL_MAX], ...
    QdScaleRange=[QD_MUL_MIN; QD_MUL_MAX], ...
    RScaleRange=[RS_MUL_MIN; RS_MUL_MAX], ...
    XScaleRange=[XS_MUL_MIN; XS_MUL_MAX], ...
    BScaleRange=[BC_MUL_MIN; BC_MUL_MAX], ...
    assertOPFOk=false);
while ~bdata.success
    bdata = makeDispData(sdata, MPOptions=mpopt, ...
        PdScaleRange=[PD_MUL_MIN; PD_MUL_MAX], ...
        QdScaleRange=[QD_MUL_MIN; QD_MUL_MAX], ...
        RScaleRange=[RS_MUL_MIN; RS_MUL_MAX], ...
        XScaleRange=[XS_MUL_MIN; XS_MUL_MAX], ...
        BScaleRange=[BC_MUL_MIN; BC_MUL_MAX], ...
        assertOPFOk=false);
end
[FBUSES, TBUSES] = deal(bdata.branch(:, 1), bdata.branch(:, 2));

%% Base problem data
pdata = buildProbData(sdata, bdata, calcCft=true, calcCu=true, calcCd=true);
[Cf, Ct] = deal(pdata.Cf, pdata.Ct);
Cft = Cf - Ct;
clear Cf Ct;

%% For anticipating supply injections
SUPREG = 1.0;
cvx_solver Mosek;

%% For solving the APF equations
REFANG = 0;
EPSTOL = 1e-5;
NITERS = 100;
NEVALS = max(1000, 500*pdata.N_b);
TRFALG = 'factorization';
LMDAMP = 1e-2;
LMDIAG = 'none';

%% Containers
BASEMVA = pdata.sys_baseMVA;
meanPds = bdata.bus(pdata.dem_buses, 3) / pdata.sys_baseMVA;
meanQds = bdata.bus(pdata.dem_buses, 4) / pdata.sys_baseMVA;
uuFeasi = zeros([NUM_SAMPLES, 1], 'logical');
uuAPFOk = zeros([NUM_SAMPLES, 1], 'logical');
uuMuls = zeros([NUM_SAMPLES, 1]);
uuPus = zeros([NUM_SAMPLES, pdata.N_u]);
uuQus = zeros([NUM_SAMPLES, pdata.N_u]);
uuVms = zeros([NUM_SAMPLES, pdata.N_b]);
uuVas = zeros([NUM_SAMPLES, pdata.N_br]);
uuPs = zeros([NUM_SAMPLES, 1]);
urFeasi = zeros([NUM_SAMPLES, 1], 'logical');
urAPFOk = zeros([NUM_SAMPLES, 1], 'logical');
urMuls = zeros([NUM_SAMPLES, 1]);
urPus = zeros([NUM_SAMPLES, pdata.N_u]);
urQus = zeros([NUM_SAMPLES, pdata.N_u]);
urVms = zeros([NUM_SAMPLES, pdata.N_b]);
urVas = zeros([NUM_SAMPLES, pdata.N_br]);
urPs = zeros([NUM_SAMPLES, 1]);
ruFeasi = zeros([NUM_SAMPLES, 1], 'logical');
ruAPFOk = zeros([NUM_SAMPLES, 1], 'logical');
ruMuls = zeros([NUM_SAMPLES, 1]);
ruPus = zeros([NUM_SAMPLES, pdata.N_u]);
ruQus = zeros([NUM_SAMPLES, pdata.N_u]);
ruVms = zeros([NUM_SAMPLES, pdata.N_b]);
ruVas = zeros([NUM_SAMPLES, pdata.N_br]);
ruPs = zeros([NUM_SAMPLES, 1]);
rrFeasi = zeros([NUM_SAMPLES, 1], 'logical');
rrAPFOk = zeros([NUM_SAMPLES, 1], 'logical');
rrMuls = zeros([NUM_SAMPLES, 1]);
rrPus = zeros([NUM_SAMPLES, pdata.N_u]);
rrQus = zeros([NUM_SAMPLES, pdata.N_u]);
rrVms = zeros([NUM_SAMPLES, pdata.N_b]);
rrVas = zeros([NUM_SAMPLES, pdata.N_br]);
rrPs = zeros([NUM_SAMPLES, 1]);

%% Main loop: ED+ with unregularized P and unregularized Q
parfor i = 1:NUM_SAMPLES
    % Normally-distributed multipliers for the anticipated demand draws
    %   The demand draws would be centered at the "base" values from above.
    dmul = normrnd(1, MUL_STD_DEV);
    ddata = makeDispData(bdata, MPOptions=mpopt, ...
        PdScaleRange=[dmul, dmul], QdScaleRange=[dmul, dmul], ...
        RScaleRange=[1; 1], XScaleRange=[1; 1], BScaleRange=[1; 1], ...
        assertOPFOk=false);
    if ~ddata.success, continue; end

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

    % Anticipated supply injections
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=0, QuReg=0);
    InitPs = 0.5 * (sum(xed.Pu_min) - sum(xed.Pu) + sum(xed.Pu_max) - sum(xed.Pu));

    % Anticipated bus voltages
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

    if ~vol.converged, continue; end
    uuFeasi(i) = ddata.success;
    uuAPFOk(i) = vol.converged;
    uuMuls(i) = dmul;
    uuPus(i, :) = xed.Pu;
    uuQus(i, :) = xed.Qu;
    uuVms(i, :) = vol.Vm;
    uuVas(i, :) = Cft * vol.Va;
    uuPs(i, :) = vol.Ps;
end
clear ddata pdata;

%% Main loop: ED+ with unregularized P and regularized Q
parfor i = 1:NUM_SAMPLES
    % Normally-distributed multipliers for the anticipated demand draws
    %   The demand draws would be centered at the "base" values from above.
    dmul = normrnd(1, MUL_STD_DEV);
    ddata = makeDispData(bdata, MPOptions=mpopt, ...
        PdScaleRange=[dmul, dmul], QdScaleRange=[dmul, dmul], ...
        RScaleRange=[1; 1], XScaleRange=[1; 1], BScaleRange=[1; 1], ...
        assertOPFOk=false);
    if ~ddata.success, continue; end

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

    % Anticipated supply injections
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=0, QuReg=SUPREG);
    InitPs = 0.5 * (sum(xed.Pu_min) - sum(xed.Pu) + sum(xed.Pu_max) - sum(xed.Pu));

    % Anticipated bus voltages
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

    if ~vol.converged, continue; end
    urFeasi(i) = ddata.success;
    urAPFOk(i) = vol.converged;
    urMuls(i) = dmul;
    urPus(i, :) = xed.Pu;
    urQus(i, :) = xed.Qu;
    urVms(i, :) = vol.Vm;
    urVas(i, :) = Cft * vol.Va;
    urPs(i, :) = vol.Ps;
end
clear ddata pdata;

%% Main loop: ED+ with regularized P and unregularized Q
parfor i = 1:NUM_SAMPLES
    % Normally-distributed multipliers for the anticipated demand draws
    %   The demand draws would be centered at the "base" values from above.
    dmul = normrnd(1, MUL_STD_DEV);
    ddata = makeDispData(bdata, MPOptions=mpopt, ...
        PdScaleRange=[dmul, dmul], QdScaleRange=[dmul, dmul], ...
        RScaleRange=[1; 1], XScaleRange=[1; 1], BScaleRange=[1; 1], ...
        assertOPFOk=false);
    if ~ddata.success, continue; end

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

    % Anticipated supply injections
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=SUPREG, QuReg=0);
    InitPs = 0.5 * (sum(xed.Pu_min) - sum(xed.Pu) + sum(xed.Pu_max) - sum(xed.Pu));

    % Anticipated bus voltages
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

    if ~vol.converged, continue; end
    ruFeasi(i) = ddata.success;
    ruAPFOk(i) = vol.converged;
    ruMuls(i) = dmul;
    ruPus(i, :) = xed.Pu;
    ruQus(i, :) = xed.Qu;
    ruVms(i, :) = vol.Vm;
    ruVas(i, :) = Cft * vol.Va;
    ruPs(i, :) = vol.Ps;
end
clear ddata pdata;

%% Main loop: ED+ with regularized P and regularized Q
parfor i = 1:NUM_SAMPLES
    % Normally-distributed multipliers for the anticipated demand draws
    %   The demand draws would be centered at the "base" values from above.
    dmul = normrnd(1, MUL_STD_DEV);
    ddata = makeDispData(bdata, MPOptions=mpopt, ...
        PdScaleRange=[dmul, dmul], QdScaleRange=[dmul, dmul], ...
        RScaleRange=[1; 1], XScaleRange=[1; 1], BScaleRange=[1; 1], ...
        assertOPFOk=false);
    if ~ddata.success, continue; end

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

    % Anticipated supply injections
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=SUPREG, QuReg=SUPREG);
    InitPs = 0.5 * (sum(xed.Pu_min) - sum(xed.Pu) + sum(xed.Pu_max) - sum(xed.Pu));

    % Anticipated bus voltages
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

    if ~vol.converged, continue; end
    rrFeasi(i) = ddata.success;
    rrAPFOk(i) = vol.converged;
    rrMuls(i) = dmul;
    rrPus(i, :) = xed.Pu;
    rrQus(i, :) = xed.Qu;
    rrVms(i, :) = vol.Vm;
    rrVas(i, :) = Cft * vol.Va;
    rrPs(i, :) = vol.Ps;
end
clear ddata pdata;

%% Outro
fname = sprintf("%s_results", SNAPSHOT_ID);
save(sprintf("./%s.mat", fname), ...
    "SNAPSHOT_ID", "NUM_SAMPLES", "MUL_STD_DEV", ...
    "PD_MUL_MIN", "PD_MUL_MAX", ...
    "QD_MUL_MIN", "QD_MUL_MAX", ...
    "RS_MUL_MIN", "RS_MUL_MAX", ...
    "XS_MUL_MIN", "XS_MUL_MAX", ...
    "BC_MUL_MIN", "BC_MUL_MAX", ...
    "FBUSES", "TBUSES", "BASEMVA", ...
    "SUPREG", "REFANG", "EPSTOL", "NITERS", "NEVALS", ...
    "TRFALG", "LMDAMP", "LMDIAG", ...
    "meanPds", "meanQds", ...
    "uuFeasi", "uuAPFOk", ...
    "uuMuls", "uuPus", "uuQus", "uuVms", "uuVas", "uuPs", ...
    "urFeasi", "urAPFOk", ...
    "urMuls", "urPus", "urQus", "urVms", "urVas", "urPs", ...
    "ruFeasi", "ruAPFOk", ...
    "ruMuls", "ruPus", "ruQus", "ruVms", "ruVas", "ruPs", ...
    "rrFeasi", "rrAPFOk", ...
    "rrMuls", "rrPus", "rrQus", "rrVms", "rrVas", "rrPs");
fprintf("Saved to './%s.mat'\n", fname);
clear;
