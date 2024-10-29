# LAAKE

        =============================================================
             __           ____       ____     __   ___   ________
            |  |         /    \     /    \   |  | /  /  |   _____|
            |  |        /  /\  \   /  /\  \  |  |/  /   |  |_____
            |  |       |  /__\  | |  /__\  | |      \   |   _____|
            |  |_____  |   __   | |   __   | |  |\   \  |  |_____
            |________| |__|  |__| |__|  |__| |__| \___\ |________|
            Lightcurve  Assembly Architecture for KAMP   Extracts
        =============================================================


Lightcurve Assembly Architecture for KAMP Extracts (LAAKE) is a python and C based data processing pipeline specifically designed for KMTNet AGN Monitoring Program (KAMP), which utilizes Korea Microlensing Telescope Network (KMTNet) to generate light curves of 500+ AGNs (Active Galactic Nuclei) in the Southern Hemisphere. 

## Installation

Install the latest development version:
```
git clone https://github.com/ericjsh/LAAKE.git
```

Resolve environment:
```
conda env create --file environment.yaml --name LAAKE
conda activate LAAKE
```

Visit Notion webpage (https://congruous-value-af8.notion.site/LAAKE-12c48de1964b80e09902ce791bc1d833) for tutorial on running data process and/or reading output files