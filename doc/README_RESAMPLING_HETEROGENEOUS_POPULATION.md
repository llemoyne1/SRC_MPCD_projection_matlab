# Weighted resampling: distributed heterogeneous population initial state

This step adds a designed Taylor--Green test with a spatially distributed initial
population heterogeneity.  It is less artificial than the two-pocket test, but it
still deliberately exercises the closed-loop population-management method.

The initial edit conserves the total number of active particles.  A target cell
population field is generated around `gamma` with a prescribed standard deviation
in particles per cell.  Particles are then moved from donor cells to receiver
cells so that the initial population follows this target field.

The intended closed-loop sequence is then:

```text
weighted SRC/Q6
save pre-edit weighted velocity
extract particles from cells with N > NMax
insert particles into cells with N < NMin
local conservative mass/momentum remap
optional weighted thermostat after remap
```

The edit is selected with:

```matlab
'initialDepletion','heterogeneous_population'
```

Key options:

```matlab
'PopulationHeterogeneityStd',4   % target std in particles/cell
'PopulationHeterogeneityMin',0   % lower target clamp
'PopulationHeterogeneityMax',40  % upper target clamp
'NMin',14
'NTarget',20
'NMax',26
'ExtractEvery',1
'InsertEvery',1
```

A visual single case can be launched with:

```matlab
outH = run_resamp_population_heterogeneity_demo( ...
    'PopulationHeterogeneityStd',4, ...
    'steps',300, ...
    'visualEvery',5, ...
    'summaryEvery',25);
```

A compact non-visual sweep can be launched with:

```matlab
outSweep = run_resamp_population_heterogeneity_sweep( ...
    'PopulationHeterogeneityStdList',[0 2 4 6 8], ...
    'steps',300, ...
    'summaryEvery',25);
```

The sweep writes:

```text
runs/resamp_population_heterogeneity_sweep/heterogeneity_sweep_summary.csv
```

The important fields are:

```text
initialNMin / initialNMax / initialNStd
extractedCumulative / insertedCumulative
NMinFinal / NMaxFinal / NStdFinal
MRelRmsFinal
mParticleRelStdFinal
kBTFinal
rmsDivAfterFinal
```

This is not a final physical validation case.  Its role is to check whether the
weighted recycling loop converges from a controlled distributed heterogeneity and
whether particle masses remain bounded as the initial heterogeneity level is
increased.
