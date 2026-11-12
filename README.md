<p align="left" >
<a href='https://carbonplan.org'>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://carbonplan-assets.s3.amazonaws.com/monogram/light-small.png">
  <img alt="CarbonPlan monogram." height="48" src="https://carbonplan-assets.s3.amazonaws.com/monogram/dark-small.png">
</picture>
</a>
</p>

# iam-rtm-compare

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

This repository contains the analysis for "Integrated assessment models overestimate carbon removal from enhanced weathering" by Tyler Kukla, Jay Fuhrman, Zach Perzan, Connor Mack, and Freya Chay. We use the rock application schedule of an integrated assessment model, [GCAM](https://github.com/JGCRI/gcam-core), to drive a reactive transport model of soil geochemistry, [SCEPTER](https://github.com/carbonplan/SCEPTER), and compare the carbon removal each model projects.

## Quick start

This repository uses functions from [ew-workflows](https://github.com/carbonplan/ew-workflows) and expects it to be cloned alongside:

```bash
git clone https://github.com/carbonplan/ew-workflows.git
git clone https://github.com/carbonplan/iam-rtm-compare.git
cd iam-rtm-compare
uv sync
```

Most scripts read from and write to CarbonPlan's S3 storage. Batch SCEPTER simulations run on [Coiled](https://docs.coiled.io/user_guide/setup/index.html) via `run/batch-driver.py`.

## License

This work is licensed under the MIT license.

## About us

CarbonPlan is a nonprofit organization that uses data and science for climate action. We aim to improve the transparency and scientific integrity of climate solutions with open data and tools. Find out more at [carbonplan.org](https://carbonplan.org/).
