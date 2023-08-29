# Use of malaria Pf7
1. Start from the VCF file directly
2. Filter by VQSR (already in the file) and Fws (monoclonal only)
3. Phasing using dominant allele (deploid method might have issue of non-overlapping sites)
4. Calling IBD (hapibd) for subsampling unrelated samples
5. Subseting dataset to build structured population
6. Subsetting by region and year to build a homogeneous population for Ne



# Pf7 metadata fields
[source](https://malariagen.github.io/parasite-data/pf7/Data_access.html#variant-calls)

- The `Sample` column gives the unique sample identifier used throughout all Pf7
analyses.
- The `Study` refers to the partner study which collected the sample.
- The `Country` & `Admin level 1` describe the location where the sample was
collected from.
- The `Country latitude`, `Country longitude`, `Admin level 1 latitude` and
`Admin 1 longitude` contain the GADM coordinates for each country &
administrative level 1.
- The `Year` column gives the time of sample collection.
- The `ENA` column gives the run accession(s) for the sequencing read data for
each sample.
- The `All samples same case` column identifies samples set collected from the
same individual.
- The `Population` column gives the population to which the sample has been
assigned. The possible values are: Africa - West (AF-W), Africa-Central (AF-C),
Africa - East (AF-E), Africa - Northeast (AF-NE), Asia - South - East (AS-S-E),
Asia - South – Far East (AS-S-FE), Asia - Southeast - West (AS-SE-W), Asia -
Southeast - East (AS-SE-E), Oceania - New Guinea (OC-NG), South America (SA).
- The `% callable` column refers to the % of the genome with coverage of at
least 5 reads and less than 10% of reads with mapping quality 0.
- The `QC pass` column defines whether the sample passed (True) or failed
(False) QC.
- The `Exclusion` reason describes the reason why the particular sample was
excluded from the main analysis.
- The `Sample type` column gives details on the DNA preparation method used
- The `Sample was in Pf6` column defines whether the sample was included in the
previous version of the data release (Pf6) or if it is new to Pf7.
