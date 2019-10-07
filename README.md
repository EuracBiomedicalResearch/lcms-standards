# Measurement of standards in different LC/MS configurations

This repository contains measurement of sets of standards on LC-MS systems with
specific configurations. The aim is to provide approximate retention times for
the standard compounds in the respective setups.

Measured retention time for each standards (in potentially different setups) are
supposed to be stored in the table
[standards_rtime.txt](data/standards_rtime.txt) while information on the
measurement run (LC-MS setup) is supposed to be put into
[measurement_run.txt](data/measurement_run.txt).


## Data tables

- [standards_rtime.txt](data/standards_rtime.txt): table containing retention
  times for each standard on each setup.
- [measurement_run.txt](data/measurement_run.txt): table providing additional
  information about a certain measurement run (eventually also LC-MS setup) in
  which standards were measured.
- [internal_standards.txt](data/internal_standards.txt): table containing the
  (exact) masses for the internal standards used in our LC-MS setup.

## Measurement runs

As up to now retention times for standards were determined by measuring them in
water (sets of standards spiked into water and measured with LC-MS). This
information was then used to define their approximate retention time in serum
samples by manually evaluating the peaks close to the retention time measured in
water. This was done by Giuseppe Paglia. Standards/compounds are supposed to
elute later in a more *rich* sample like serum. Thus, the retention times for
water and serum might differ.
