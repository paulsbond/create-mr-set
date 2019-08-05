# paper3

Data for PhD paper 3 test set

## Requirements

CCP4 7.0.076
Python 3.6.8 with biopython-1.74

## 1. Choose structures

Ran `choose_structures /data/pdb/validation_reports` with default parameters.
This chose 2000 targets and wrote them to `chosen_targets`.
The targets don't share chains with 50% sequence idenity.
They are distributed evenly in 10 bins between 1.0 and 3.5 A resolution.
The PDB validation sliders were checked
relative to similar resolution structures.
R-free had to be at least the 50th percentile
and the other 4 sliders had to be in at least the 45th percentile.

## 2. Refine deposited structure

Run the commands in `refine_commands` to prepare a standardised MTZ file
and refine the deposited PDB file against it.
This should create the following files:

- `joined.mtz` with FREE, FP, SIGFP and reference ABCD
- `metadata.json` with cell info, reference R-factors, resolution and semet
- `reference.pdb` the refined deposited PDB without UNL residues

After running, run `refine_check`
to get a list of structures that failed in `refine_failed`
and a list of structures that passed in `refine_passed`.
Check any structures that errored for an unknown reason,
and update the reasons in `refine_failed` for future reference.
Finally, run `refine_clean` to delete files
that are no longer needed for the next steps.

## 3. Extract unique protein sequences

Run the script `extract_protein_sequences`,
which creates `protein.fasta` files in each directory listed in `refine_passed`.
The file contains unique annotated protein sequences with 20+ amino acids.
The ID of each chain is a comma separated list of IDs with that sequence.

## 4. Find structural homologues

## 5. Molecular replacement (or just gesamt depending on output from 4)
