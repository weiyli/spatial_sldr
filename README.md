The spatial range exponent in human mobility behaviour
================
2025-09-18

# Code

All code used in this study is available at
<https://github.com/weiyli/spatial_sldr>.  
The repository contains scripts to reproduce all results, including
network construction, metric computation, simulation, and figure
generation. The code is organized into the following modules:

## Global configuration

- `sldr_global_vars_funs.R`: Provides global parameters and shared
  settings required across the entire SLDR analysis workflow. This
  script is sourced at the beginning of all modules to ensure consistent
  definitions of study regions, disaster events, and analysis windows.

  Specifically, it:

  - Loads required packages for **data processing** (`data.table`,
    `plyr`, `dplyr`) and **visualization** (`ggplot2`, `patchwork`,
    `scales`, `RColorBrewer`, `ggthemes`, `latex2exp`, etc.).
  - Defines the set of **metropolitan statistical areas (MSAs)** under
    study, including their names, numeric identifiers, states, counties,
    and FIPS codes.
  - Defines the set of **disaster-affected regions**, with corresponding
    identifiers and metadata.
  - Provides **event labels** (e.g., `COVID-19`, `Hurricane_Harvey`,
    `Storm_Texas`) aligned with MSAs or disaster areas.
  - Sets up **analysis time windows**, including baseline durations
    (`startdate.msa`, `durdate.msa`, `enddate.msa`).

These variables act as the backbone of the project, ensuring that
downstream scripts for spatial range exponent fitting, spatial
regression, and simulation all refer to a unified set of identifiers and
dates.

## Data processing

- `od_flow.R`: Processes SafeGraph mobility data to construct
  origin-destination (OD) flows at the MSA and county levels. For each
  Census Block Group (CBG), the script computes daily flow volumes and
  degree-based connectivity measures.

  **Outputs:**

  - `Intra_Flow_<date>.csv`: daily CBG-level OD flows within each MSA or
    county  
  - `Intra_Degree_<date>.csv`: daily CBG-level in-/out-degree within
    each MSA or county

- `select_region.R`: Defines study regions by mapping counties to each
  MSA and extracting all CBGs within their boundaries. ACS 2015-2019
  demographic attributes and geometric features are merged with CBG
  polygons to produce standardized inputs for spatial modeling.

  **Outputs:**

  - `<MSA>.geojson`: per-MSA CBG layers with population, income, and
    area attributes  
  - `selected_msa_county.csv`: county roster associated with each MSA

## SLDR model fitting

- `sldr_fit_msa.R`: Implements the Spatial Lag Decay Autoregressive
  (SLDR) model fitting procedure at the **MSA level**. This script
  carries out the following steps:

  - **Setup**: Defines working directories, loads required packages
    (`sf`, `spdep`, `igraph`), and sources `sldr_global_vars_funs.R` for
    global variables and region/event definitions.  
  - **Region specification**: Loads the set of MSAs and
    disaster-affected regions. County–MSA mappings are imported from
    `selected_msa_county.csv`, with standardized FIPS codes.  
  - **Adjacency construction**: For each region, reads CBG-level
    polygons (`.geojson`), extracts CBG identifiers and attributes
    (population, area), and constructs queen-contiguity adjacency
    matrices. The maximum spatial lag (`lag.max`) is determined from the
    network diameter, and multi-lag adjacency matrices are generated
    (`queen.adj`).  
  - **Parameter estimation**: For each day and each CBG, constructs
    mobility matrices from flow data. Estimates optimal spatial range
    exponents ($\rho$) under three functional forms:
    - no-decay process  
    - power-law decay  
    - exponential decay  
  - **Model fitting**: Using the estimated $\rho$, generates fitted
    mobility flows for each CBG under the three functional forms.  
  - **Model evaluation**: Computes model performance at the MSA level,
    including:
    - $R^2$ for each functional form  
    - Correlations between empirical and simulated flows

  **Outputs:**

  - `*_SLDR_params_*.csv`: estimated parameters ($\rho$)  
  - `*_SLDR_fit_*.csv`: fitted vs. empirical flows  
  - `SLDR_r2_*.csv`: goodness-of-fit statistics  
  - `SLDR_rank_cor_*.csv`: correlation metrics

- `sldr_fit_county.R`: Extends the above SLDR fitting procedure to the
  **county level**, applying the same workflow (adjacency construction,
  parameter estimation, model fitting, and evaluation) but at a finer
  geographic scale.

Together, these scripts provide the core SLDR estimation pipeline at
multiple spatial resolutions. Their outputs (parameters, fitted flows,
and evaluation metrics) serve as the foundation for subsequent
visualization and comparative analysis.

## Spatially correlated SIR simulations

- `spreading_model.R`: Simulates the **spatially correlated SIR
  (Susceptible-Infected-Recovered) model** across MSAs while varying the
  spatial range exponent ($\rho$). The script evaluates how spatial
  dependence modulates disease dynamics under different epidemiological
  parameters. This script carries out the following steps:
  - **Setup**: Defines paths, loads packages (`igraph`, `doParallel`,
    `foreach`), and imports region definitions from
    `sldr_global_vars_funs.R`.  
  - **Functions**:
    - `row_normalize()`: Row-normalizes spatial weight matrices.  
    - `update_SIR()`: Updates SIR states based on weighted adjacency and
      epidemiological parameters.  
    - `run_SIR_simulation()`: Runs repeated simulations over time,
      incorporating spatial lag decay weights parameterized by $\rho$.  
  - **Simulation design**:
    - Constructs adjacency and distance matrices from CBG polygons.  
    - Applies exponential decay with $\rho$ to generate spatial weight
      matrices.  
    - Initializes infections with $N_0$ nodes and infection probability
      $I_0$.  
    - Runs multiple stochastic simulations and averages results.  
  - **Parameter variation**:
    - Varies $\rho$ around empirical maxima (e.g., ±25% or ±50%).  
    - Uses selected infection ($\beta$) and recovery ($\mu$) rates.  
    - Repeats simulations for all parameter combinations across MSAs.  
    - Saves simulated epidemic trajectories with time series of
      susceptible, infected, and recovered populations.

  **Outputs:**
  - `SIR_model_*.csv`: Averaged epidemic trajectories across time for
    each MSA under varied $\rho,\beta,\mu$.  
  - `SIR_model_rho_*.csv`: Simulation results highlighting sensitivity
    to perturbations in $\rho$.

These outputs are used to generate **Fig. 2b-c** and **Fig. S9**,
showing infection curves under varying spatial range exponents and
epidemiological parameters, and to support the robustness analysis of
spatial dependence in epidemic spread.

## Figure generation

- `fig1_a_c.R`: Produces adjacency and decay illustrations for **Fig.
  1a-c**. It visualizes CBG-level queen neighborhoods and
  first-/second-order spatial lags within a selected MSA, and plots the
  exponential spatial decay curve used as an example.

  **Outputs:**

  - `lag1.pdf`: edges for 1-lag (queen) adjacency among CBG centroids  
  - `lag2.pdf`: edges for 2-lag adjacency among CBG centroids  
  - `queen.pdf`: 5×5 schematic of queen neighborhood (lag = 0/1/2)  
  - `lag.h.pdf`: panel combining 1-lag and 2-lag edge maps  
  - `brief_lag_exp.pdf`: exemplar exponential decay curve

- `fig1_d_i.R`: Generates **Fig. 1d-i**, visualizing the performance of
  the SLDR model at both MSA and county scales.

  - **Parameter distributions:** Ranks MSAs by the spatial range
    exponent ($\rho$), with error bars indicating min/median/max values
    across days.  
  - **Population dependence:** Plots $\rho$ against population size,
    combining MSA-level estimates (large points) and county-level
    estimates (small points).  
  - **Flow maps:** Compares empirical and fitted mobility flows side by
    side for selected MSAs, with one-lag queen adjacency networks as the
    geographic background.

  **Outputs:**

  - `brief_rho_fit_*.pdf`: panel combining parameter distributions,
    population dependence, and empirical vs. fitted flow maps.

- `fig2.R`: Produces **Fig. 2**. The script computes the
  **before→during** percentage change in $\rho$, and visualizes
  epidemiological implications with a spatially correlated SIR model
  (peak sensitivity and time paths under $\Delta\rho$).

  **Outputs:**

  - `us-msas-daily-cumulative-case.csv`: daily cumulative and new cases
    aggregated to MSAs  
  - `brief_rho_case_*.pdf`: scatter of % change in $\rho$ versus
    cumulative cases per population (used alongside panel **b**)  
  - `brief_rho_pop_mob_sir_*.pdf`: panel combining **(a)** ranked %
    change in $\rho$ across MSAs (with IQR bars and event labels),
    **(b)** peak-height sensitivity to $\Delta\rho$, and **(c)** SIR
    infection trajectories for exemplar MSAs under selected $\Delta\rho$

- `fig_S1.R`: Generates **Fig. S1** showing how neighborhood scope
  expands with spatial lag under queen contiguity across MSAs.

  - **Fig. S1a:** Number of neighboring CBGs by lag.  
  - **Fig. S1b:** Cumulative population of neighbors by lag
    (millions).  
  - **Fig. S1c:** Cumulative land area of neighbors by lag (km²).  
    Here, each lag corresponds to one additional queen step outward from
    each CBG (lag 1 = first-order neighbors, lag 2 = second-order,
    etc.).

  **Outputs:**

  - `SI_queen_lag.pdf`: panel combining (a) neighbors, (b) cumulative
    population, and (c) cumulative area by lag for all MSAs.

- `fig_S2_S3.R`: Generates **Fig. S2-S3**, comparing SLDR model accuracy
  across MSAs using **SMAPE** (Symmetric Mean Absolute Percentage Error)
  by day. The script reads daily fitted vs. empirical flows from
  `sldr_fit_msa.R`, computes SMAPE for three functional
  forms-**no-decay**, **power-law decay**, and **exponential decay**-and
  facetes results by MSA. Plots are shown on a log-scaled y-axis,
  separately for **before** and **during** periods.

  **Outputs:**

  - `SMAPE_*.csv`: per-day SMAPE for each model and MSA
  - `SI_fit_model_before_*.pdf`: Fig. S2 - SMAPE time series (before
    period) across MSAs
  - `SI_fit_model_during_*.pdf`: Fig. S3 - SMAPE time series (during
    period) across MSAs

- `fig_S4_S5.R`: Generates **Fig. S4-S5**, contrasting empirical
  mobility with SLDR outcomes and checking residual stability by day.

  - **Fig. S4 (maps):** For each MSA, selects the representative day
    (most frequent maximum $R^2$), overlays the **one-lag queen**
    adjacency, and shows **Empirical** vs. **Model (exponential decay)**
    node maps with size/color proportional to normalized flows.  
  - **Fig. S5 (residuals):** Computes CBG-level residuals
    $R=empirical-model$ and, per MSA, plots **violin + box**
    distributions by day with pairwise $t$-test annotations to assess
    temporal stability of the spatial range exponent $\rho$.

  **Outputs:**

  - `SI_fit_*.pdf`: Fig. S4 - empirical/model maps for all MSAs  
  - `SI_fit_res_*.pdf`: Fig. S5 - residual violin/box plots with
    significance marks across days

- `fig_S6_S9.R`: Generates **Fig. S6–S9**, combining COVID-19 case
  alignment, OD distance distributions, commuting correlations, and SIR
  simulations. The script performs the following tasks:

  - **COVID-19 case alignment**: Aggregates county-level confirmed cases
    to MSAs, aligns them with daily SLDR estimates of $\rho$, and plots
    percentage change in $\rho$ vs. cumulative confirmed cases per
    capita (log scale). (**Fig. S6**)  
  - **OD distance distributions**: Computes weighted average OD travel
    distances across MSAs, using visit counts as weights, and visualizes
    their distributions. (**Fig. S7**)  
  - **Commuting correlations**: Calculates Spearman’s rank correlation
    between commuting mode shares and the spatial range exponent $\rho$,
    highlighting sustainable vs. car-based modes. (**Fig. S8**)  
  - **SIR simulations**: Runs spatially correlated SIR models under
    varying epidemiological parameters ($\beta,\mu$) across selected
    MSAs, showing infection trajectories for different $\rho$ values.
    (**Fig. S9**)

  **Outputs:**

  - `us-msas-daily-cumulative-case.csv`: Daily cumulative and new cases
    aggregated to MSAs  
  - `brief_rho_case_*.pdf`: **Fig. S6** - Scatter of % change in $\rho$
    vs. cumulative cases per population  
  - `SI_od_distance_*.pdf`: **Fig. S7** - Distribution of weighted OD
    travel distances across MSAs  
  - `SI_commute_rho_corr_*.pdf`: **Fig. S8** - Correlation between
    commuting mode shares and $\rho$  
  - `SI_sir_beta_mu_*.pdf`: **Fig. S9** - SIR trajectories under
    different $\beta,\mu,\rho$ combinations

**Note on data access**:  
The raw mobility data used in this study are derived from SafeGraph,
which are subject to access restrictions. We therefore provide only the
processed, results-level outputs (e.g.,
`*_SLDR_fit_*.csv`,`SIR_model_rho_*.csv`), which are sufficient to
reproduce all figures and analyses presented in this study.
