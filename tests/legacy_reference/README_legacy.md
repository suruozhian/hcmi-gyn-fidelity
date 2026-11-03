# Legacy reference file (test use only; not an input of any analysis)

`LEGACY_uterine_units_plus_nextgen.csv` is the intermediate file formerly shipped as `inputs/intermediate/uterine_units_plus_nextgen.csv`
and published as `supp/source_data/static/uterine_units_plus_nextgen.csv` (until round 5b, 2026-10-08). It was written by an
exploratory script that is not part of the package. Since round 5c it is replaced by `uterine_units_classes.csv`, which
`05_derive_endometrial_classes.py` derives from open data with the final unit definitions (`inputs/units_215.csv`).

It is kept only so that `tests/compare_molecular.py` can show what changed:
- the exploratory script merged the E-MTAB-10664 (Berg) samples per patient (units `Berg-pt<NN>`; unit histology taken from
  the first sample; unit MLH1 = mean over the patient's samples). For patient 04 its row 92 (`Berg-pt04`, `UT_endometrioid`,
  from sample OEC-04-E) was therefore inconsistent with Supplementary Tables S2 and S10, which use the final unit OEC-04-CC
  (`UT_clear_cell`);
- the surrogate class of every unit (102/102), and hence every published number, is unchanged.
