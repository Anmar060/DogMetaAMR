# Kraken2 input used by CUPID

The CUPID host-analysis stage reused existing standard per-read Kraken2 output rather than rerunning Kraken2.

Recorded barcode36 input:

`/exports/cmvm/eddie/eb/groups/mazeri_grp/Ana/scripts_methyl/krakens/kraken2_standard_c05_barcode36.txt`

Expected naming:

`kraken2_standard_c05_barcodeXX.txt`

The standard per-read output is required because the report file does not retain read IDs.

Kraken2 version recorded: 2.1.3.
