# Paper: Test Set

Data for PhD test set paper.

## Requirements

CCP4 7.0.076
Python 3.6.8 with biopython-1.74

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

The script `find_homologues` (rename from run_gesamt)
will run GESAMT to search for structural homologues.
It needs the path to a GESAMT archive, which can be created using:

```bash
gesamt --make-archive
```

This script loops through the structures in `refine_passed` one at a time
but will parallelise each search over the threads available on the machine.
It creates files called `chain?_gesamt.txt`
where `?` is the chain ID for each representative chain in `protein.fasta`.

## 5. Choose models

The `choose_models` script looks at the hits in the `chain?_gesamt.txt` files
and chooses which ones should be made into models.
It creates a folder called `models` in each directory
and a folder within that for each model with the format
`{query chain}_{hit pdb}{hit chain}`,
e.g. `A_3ecrB` means that chain B of 3ECR
has been chosen as an MR model for chain A.

There is a `choose_models_check` script that can be run after
and will list PDBs for which no models were chosen.

## 7. Molecular replacement

Script `mr` needs to be run once for each model.
It runs the following steps:

1. gesamt to superpose hit chain over query chain and write out an alignment.  
2. phaser.sculptor to trim the hit chain using the gesamt alignment.  
3. phaser to do MR with the sculpted model.
4. refmac to refine the top solution from phaser.

## 8. Extract metadata
