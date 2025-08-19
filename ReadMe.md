# Scenario-Based Projection Tool for Pandemic Spread

## Overview

This is the repository for the scenario-based projection tool concerning the impact of pandemic spread. This tool simulates COVID-19 infections in Germany (2020–2021) using age- and sex-specific weekly incidence from RKI SurvStat and a simple MicSim state model (`Susceptible → Infected → Recovered`). The tool is developed within the Leibniz Lab for Pandemic Preparedness by the team at DIW (Alex Lepe and Sabine Zinn).

## Repository contents

-   app.R - Shiny tool for the simulation
-   00_prepping_rki_data.qmd - Code to load and clean the RKI data
-   00_prepping_rki_data.html - Rendered HTML of the quarto file
-   data/
    -   SurvStat_complete_data - All the processed data in both csv and rds formats
    -   SurvStat_female_data - Processed data for females in both csv and rds format
    -   SurvStat_male_data - Processed data for males in both csv and rds format
    -   weekly - Folder with all the raw weekly data files from SurvStat\@RKI 2.0
-   shinycovid-simulation.png – Example output (screenshot)
-   README.md - Description of the project
-   ProjectionPandemicIncidence.Rproj - R project

## Running the tool

### Prerequisites

Ensure the following packages are installed: `shiny`, `shinyBS`, `MicSim`, `data.table`, and `ggplot2`

### Running the tool

From the project root run:

``` r
shiny::runApp("app.R")
```

Adjust the sliders in the left-hand panel and click **Run simulation**. **Note:** larger values of **N** (cohort size) increase runtime.

If you run the tool with **all default settings**, the output should match the screenshot ![Simulation output](shinycovid-simulation.png)

## How the projection tool works

1.  Import infection rates: Weekly SARS‑CoV‑2 infection rates for 2020–2021 are read from `data/weekly/` (source: Robert Koch Institute – SurvStat\@RKI 2.0).
2.  Generate a synthetic cohort: `N` individuals are sampled with uniform ages between the user‑chosen age range and the sample is evenly split by sex.
3.  Microsimulation (using MicSim): Each individual can move through the states `Susceptible → Infected → Recovered`
    -   Transition rates are scaled by `1 – (intervention effectiveness / 100)`, which is defined by the user.
    -   The amount of time spent being ill is also defined by the user.
    -   For more details on how MicSim works, see the [CRAN documentation](https://cran.r-project.org/web/packages/MicSim/index.html).
4.  Post‑processing: Getting daily counts for
    -   active infections (14‑day rolling sum),
    -   total number of people recovered from illness
    -   remaining number of people susceptible to infection
5.  Visualisation: The Shiny UI plots the percentage of people actively infected and recovered per day.

## Inputs (left panel)

-   **Sample size (N):** 1,000–10,000
-   **Minimum/Maximum age at start:** age bounds for the initial cohort
-   **Duration of illness (days):** used to define recovery timing
-   **Intervention effectiveness (%):** multiplies incidence by `1 − x/100`
-   **Run simulation:** executes the model with current settings
