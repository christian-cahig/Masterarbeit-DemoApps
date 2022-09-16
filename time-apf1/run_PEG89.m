%% Initializations
SNAPSHOT_ID = 'PEG89';

% For generating dispatch data sets 
NUM_INSTANCES = 1260;
[PD_MUL_MIN, PD_MUL_MAX] = deal(0.9, 1.1);
[QD_MUL_MIN, QD_MUL_MAX] = deal(0.9, 1.1);
[RS_MUL_MIN, RS_MUL_MAX] = deal(0.95, 1.05);
[XS_MUL_MIN, XS_MUL_MAX] = deal(0.95, 1.05);
[BC_MUL_MIN, BC_MUL_MAX] = deal(0.95, 1.05);

% For anticipating supply injections
PU_REG = 1e-3;
QU_REG = 1e-3;

% Go-to MATPOWER options
mpopt = getGoToMPOpts();

%% Snapshot data
load(sprintf("../snapshots/%s.mat", SNAPSHOT_ID));

%% Containers
areFeasible = zeros([NUM_INSTANCES, 1], "logical");
sdpSteps = zeros([NUM_INSTANCES, 1]);
sdpTimes = zeros([NUM_INSTANCES, 1]);
sdmSteps = zeros([NUM_INSTANCES, 1]);
sdmTimes = zeros([NUM_INSTANCES, 1]);

%% Warmup loop
parfor i = 1:17
    % Dispatch data
    ddata = makeDispData(sdata, MPOptions=mpopt, ...
        PdScaleRange=[PD_MUL_MIN; PD_MUL_MAX], ...
        QdScaleRange=[QD_MUL_MIN; QD_MUL_MAX], ...
        RScaleRange=[RS_MUL_MIN; RS_MUL_MAX], ...
        XScaleRange=[XS_MUL_MIN; XS_MUL_MAX], ...
        BScaleRange=[BC_MUL_MIN; BC_MUL_MAX], ...
        assertOPFOk=false);
    if ~ddata.success, continue; end

    % Problem data
    pdata = buildProbData(sdata, ddata, RefAngle=0);
    [Ph, Qh] = calcPhQh(pdata.dispatch, Vm=pdata.snapshot.bus(:, 8));
    [Po, Qo] = calcPoQo(pdata.dispatch, Vm=pdata.snapshot.bus(:, 8), Va=pdata.snapshot.bus(:, 9));

    % Anticipate supply injections via Extended Economic Dispatch
    cvx_solver SDPT3
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=PU_REG, QuReg=QU_REG);
    cvx_solver SeDuMi
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=PU_REG, QuReg=QU_REG);
end

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
    if ddata.success, areFeasible(i) = true; end
    if ~ddata.success
        [sdpSteps(i), sdpTimes(i)] = deal(NaN);
        [sdmSteps(i), sdmTimes(i)] = deal(NaN);
        continue;
    end

    % Problem data
    pdata = buildProbData(sdata, ddata, RefAngle=0);
    [Ph, Qh] = calcPhQh(pdata.dispatch, Vm=pdata.snapshot.bus(:, 8));
    [Po, Qo] = calcPoQo(pdata.dispatch, Vm=pdata.snapshot.bus(:, 8), Va=pdata.snapshot.bus(:, 9));

    % Anticipate supply injections via Extended Economic Dispatch
    cvx_solver SDPT3
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=PU_REG, QuReg=QU_REG);
    sdpSteps(i) = xed.Steps;
    sdpTimes(i) = xed.Time;
    cvx_solver SeDuMi
    xed = solveExtdEconDisp(pdata, Ph=Ph, Qh=Qh, Po=Po, Qo=Qo, PuReg=PU_REG, QuReg=QU_REG);
    sdmSteps(i) = xed.Steps;
    sdmTimes(i) = xed.Time;
end
clear sdata ddata pdata Ph Qh Po Qo xed;

%% Outro
fname = sprintf("%s_results-%f", SNAPSHOT_ID, sum(areFeasible) / NUM_INSTANCES);
save(sprintf("./%s.mat", fname), ...
    "SNAPSHOT_ID", "NUM_INSTANCES", ...
    "PD_MUL_MIN", "PD_MUL_MAX", ...
    "QD_MUL_MIN", "QD_MUL_MAX", ...
    "RS_MUL_MIN", "RS_MUL_MAX", ...
    "XS_MUL_MIN", "XS_MUL_MAX", ...
    "BC_MUL_MIN", "BC_MUL_MAX", ...
    "PU_REG", "QU_REG", ...
    "areFeasible", ...
    "sdpSteps", "sdpTimes", ...
    "sdmSteps", "sdmTimes");
fprintf("Saved to './%s.mat'\n", fname);
clear;
