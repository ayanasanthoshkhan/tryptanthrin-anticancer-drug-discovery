# Tryptanthrin Analog Library Generation

## Overview
This module focuses on generating and organising a computational library of tryptanthrin analogs for downstream drug-likeness filtering, ADMET assessment, bioactivity prediction and molecular docking.

## Workflow
1. Parent tryptanthrin structure verification.
2. Reference analog generation.
3. Exploratory analog library generation.
4. Expanded library generation.
5. Broad scaffold enumeration.
6. Synthetic accessibility (SA) scoring.

## Generated Library
- 11,393 unique analog structures.
- Raw and cleaned analog datasets.
- Molecular structures and supporting metadata.
- Synthetic accessibility scores.

## Tools
- Python
- RDKit
- Pandas

## Folder Structure
- `notebooks/` — generation and scoring notebooks.
- `data/` — generated analog datasets and supporting archives.

## Downstream Integration
The generated analog library is provided for drug-likeness and ADMET filtering, followed by target-specific bioactivity prediction and molecular docking.

## Notes
This work represents computational analog generation and prioritisation. Predicted properties and scores do not establish experimental anticancer activity.
