%% Initializations
SNAPSHOT_ID = 'ACT500';

% For generating dispatch data sets
NUM_INSTANCES = 490;
[PD_MUL_MIN, PD_MUL_MAX] = deal(0.9, 1.1);
[QD_MUL_MIN, QD_MUL_MAX] = deal(0.9, 1.1);
[RS_MUL_MIN, RS_MUL_MAX] = deal(0.95, 1.05);
[XS_MUL_MIN, XS_MUL_MAX] = deal(0.95, 1.05);
[BC_MUL_MIN, BC_MUL_MAX] = deal(0.95, 1.05);

% Go-to MATPOWER options
mpopt = getGoToMPOpts();

%% Snapshot data
load(sprintf("../snapshots/%s.mat", SNAPSHOT_ID));
[Cf, Ct] = makeCfCt(sdata);
Cft = Cf - Ct;
clear Cf Ct;

% Set supply costs
sdata.gencost(:, 5) = 1;
sdata.gencost(:, 6) = 1e-5;
sdata.gencost(:, 7) = 0;

%% For solving the extended DCOPF
cvx_solver Mosek;

%% For solving the APF equations
REFANG = 0;
EPSTOL = 1e-3;
NITERS = 100;
NEVALS = max(1000, 500 * size(sdata.bus, 1));
TRFALG = 'factorization';
LMDAMP = 1e-2;
LMDIAG = 'none';
Pu_min = sdata.gen(:, 10) ./ sdata.baseMVA;
Pu_max = sdata.gen(:, 9) ./ sdata.baseMVA;
Pu_min_sum = sum(Pu_min);
Pu_max_sum = sum(Pu_max);
PsDistrib = Pu_max - Pu_min;
PsDistrib = PsDistrib / sum(PsDistrib);
clear Pu_min Pu_max;

%% Containers
apxPus = zeros([NUM_INSTANCES, size(sdata.gen, 1)]);
apxQus = zeros([NUM_INSTANCES, size(sdata.gen, 1)]);
apxVms = zeros([NUM_INSTANCES, size(sdata.bus, 1)]);
apxVas = zeros([NUM_INSTANCES, size(sdata.bus, 1)]);
apxVbs = zeros([NUM_INSTANCES, 1]);
apfOks = zeros([NUM_INSTANCES, 1], "logical");
apfVms = zeros([NUM_INSTANCES, size(sdata.bus, 1)]);
apfVas = zeros([NUM_INSTANCES, size(sdata.bus, 1)]);
apfSks = zeros([NUM_INSTANCES, 1]);
apfVbs = zeros([NUM_INSTANCES, 1]);
difVms = zeros([NUM_INSTANCES, 1]);
difVas = zeros([NUM_INSTANCES, 1]);

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
    while ~ddata.success
        ddata = makeDispData(sdata, MPOptions=mpopt, ...
            PdScaleRange=[PD_MUL_MIN; PD_MUL_MAX], ...
            QdScaleRange=[QD_MUL_MIN; QD_MUL_MAX], ...
            RScaleRange=[RS_MUL_MIN; RS_MUL_MAX], ...
            XScaleRange=[XS_MUL_MIN; XS_MUL_MAX], ...
            BScaleRange=[BC_MUL_MIN; BC_MUL_MAX], ...
            assertOPFOk=false);
    end

    % Problem data
    pdata = buildProbData(sdata, ddata, RefAngle=REFANG, calcCu=true, calcCd=true, calcY=true);

    % APF essentials
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

    % Approximated supply injections and bus voltages
    apx = runAugmDCOPF(pdata, MPOptions=mpopt);
    apxPus(i, :) = apx.Pu;
    apxQus(i, :) = apx.Qu;
    apxVms(i, :) = apx.Vm;
    apxVas(i, :) = apx.Va;
    apxVbs(i) = norm(Cft * apx.Va, 1);

    % Nearest power flow feasible point via APF equations
    config.Algorithm = 'trust-region-dogleg';
    vol = solveAPFEqns(pdata, config, ...
        Pu=apx.Pu, Qu=apx.Qu, InitVm=apx.Vm, InitVa=apx.Va, ...
        InitPs=0.5 * (Pu_min_sum + Pu_max_sum - (2 * sum(apx.Pu))), ...
        PsDistribs=PsDistrib);
    if ~vol.converged
        config.Algorithm = 'levenberg-marquardt';
        vol = solveAPFEqns(pdata, config, ...
            Pu=apx.Pu, Qu=apx.Qu, InitVm=apx.Vm, InitVa=apx.Va, ...
            InitPs=0.5 * (Pu_min_sum + Pu_max_sum - (2 * sum(apx.Pu))), ...
            PsDistribs=PsDistrib);
    end
    if ~vol.converged, continue; end
    apfOks(i) = vol.converged;
    apfVms(i, :) = vol.Vm;
    apfVas(i, :) = vol.Va;
    apfSks(i) = vol.Ps;
    apfVbs(i) = norm(Cft * vol.Va, 1);

    % Metrics
    difVms(i) = norm(vol.Vm - apx.Vm, 2);
    difVas(i) = norm(Cft * (vol.Va - apx.Va), 2);
end
clear sdata ddata pdata apx vol;

%% Outro
fname = sprintf("%s_results-%.5f", SNAPSHOT_ID, sum(apfOks) / NUM_INSTANCES);
save(sprintf("./%s.mat", fname), ...
    "SNAPSHOT_ID", "NUM_INSTANCES", ...
    "PD_MUL_MIN", "PD_MUL_MAX", ...
    "QD_MUL_MIN", "QD_MUL_MAX", ...
    "RS_MUL_MIN", "RS_MUL_MAX", ...
    "XS_MUL_MIN", "XS_MUL_MAX", ...
    "BC_MUL_MIN", "BC_MUL_MAX", ...
    "REFANG", "EPSTOL", "NITERS", "NEVALS", ...
    "TRFALG", "LMDAMP", "LMDIAG", ...
    "Cft", "PsDistrib", ...
    "apxVms", "apxVas", "apxPus", "apxQus", "apxVas", "apxVbs", ...
    "apfOks", "apfVms", "apfVas", "apfSks", "apfVbs", ...
    "difVms", "difVas");
fprintf("Saved to './%s.mat'\n", fname);
clear;
