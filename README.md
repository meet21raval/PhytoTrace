
# Plant Phytoconstituent Retriever — Version 1

## What it does

Enter a scientific plant name, for example:

    Rheum palmatum

Version 1 retrieves plant–phytochemical associations from the current public IMPPAT website and exports an Excel workbook.

The workbook contains:

- `Summary`
- `IMPPAT`
- `NPASS`
- `Combined`

The IMPPAT sheet includes, where available:

- Plant
- Plant part
- IMPPAT ID
- Phytochemical name
- Synonyms
- External identifiers
- PubChem CID
- SMILES
- InChI
- InChIKey
- DeepSMILES
- Functional groups
- Chemical classification
- NP-Likeness
- RDKit-computed molecular formula
- Molecular weight
- Exact mass
- XlogP
- TPSA
- HBD/HBA
- Rotatable bonds
- Heavy atoms
- Carbon/nitrogen/oxygen/sulfur/halogen counts
- Ring count
- Aromatic ring count
- Fraction Csp3
- QED

## Important database-version note

This application is designed for the current IMPPAT 3.x web structure. IMPPAT's current website reports version 3.0, released September 5, 2026.

Do not hard-code an expected number such as 132. The number of records can change between database releases, and the plant page can contain multiple plant-part associations and duplicate chemical names with different IDs. Version 1 preserves the database IDs and plant-part associations and separately reports unique IDs.

## NPASS

NPASS currently provides downloadable 2026 datasets. The application attempts to discover the official download links automatically and cache the required General Information, Structure, and Species Source files under:

    NPASS_cache/

Because NPASS periodically changes download endpoints/file formats, the program does not fabricate NPASS records. If automatic NPASS download fails, download the official NPASS 2026 files manually from its download page and put the relevant files in `NPASS_cache`, then rerun.

## Installation — Windows

Install Python 3.11 or newer.

Open Command Prompt in this folder:

    python -m pip install -r requirements.txt

If `rdkit` cannot be installed by pip on your machine, use conda/mamba:

    conda create -n plantphyto python=3.11
    conda activate plantphyto
    conda install -c conda-forge rdkit pandas openpyxl requests beautifulsoup4 lxml

Then run:

    python main.py

## Output

By default:

    Rheum_palmatum_phytoconstituents.xlsx

The program uses live database pages and requires an internet connection.

## Data integrity

The application does not generate phytoconstituent names, IDs, or SMILES when the source database does not provide them. RDKit-derived descriptors are explicitly labelled as RDKit values.

## Official sources

IMPPAT:
https://cb.imsc.res.in/imppat/

IMPPAT download:
https://cb.imsc.res.in/imppat/download

NPASS:
https://bidd.group/NPASS/index.php

NPASS download:
https://bidd.group/NPASS/downloadnpass.html
