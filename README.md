# mapping_toolkit

This repository contains python scripts to help process mapping files (BED, GFF, BAM).

## bed_to_gff.py 

This script takes a BED file and converts it in a GFF formatted file. Works on BED12 and BED6.

```
usage: bed_to_gff.py [-h] --bed_file BED_FILE --source SOURCE --mol_type
                     MOL_TYPE [--is_bed12] [--make_gff3]

Creates GFF file from a given BED file. Note that the features of the GFF are
created based on the ID of the BED file.

optional arguments:
  -h, --help            show this help message and exit
  --bed_file BED_FILE, -f BED_FILE
                        Name of the BED file to be converted
  --source SOURCE, -s SOURCE
                        Name of the source
  --mol_type MOL_TYPE, -m MOL_TYPE
                        Name of the molecular type of elements from BED
  --is_bed12            Specify this argument if bed file is bed12 formated
                        and contain blocks
  --make_gff3           Specify if you want to make the output a proper gff3
  ```

  