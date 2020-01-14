# mapping_toolkit

This repository contains python scripts to help process mapping files (BED, GFF, BAM).

## gff_to_bed.py 

Reverse complement of the bed_to_gff script. This script takes a GFF file and converts it into a BED formatted file (BED12 or BED6). Selection on the molecular type and feature type to extract in the arguments.

```
usage: gff_to_bed.py [-h] --gff_file GFF_FILE [--bed12] [--no_bed6]
                     [--mol_type MOL_TYPE] [--feature_type FEATURE_TYPE]
                     [--id_as_features] [--path PATH] [--name NAME]
                     [--verbose] [--discard] [--skip_exon_number]

example run : ./gff_to_bed.py -f test_files/tiny_dmel_sample_r5-57.genes.gff -b12 -v

Creates BED files from a given GFF file with specific filters. In BED12,
groups all the elements of a selected molecular type according to their
feature type.

optional arguments:
  -h, --help            show this help message and exit
  --gff_file GFF_FILE, -f GFF_FILE
                        Name of the GFF file to be converted
  --bed12, -b12         Creates the corresponding BED12 file
  --no_bed6, -nb6       Prevents from creating the corresponding BED6 file
  --mol_type MOL_TYPE, -mt MOL_TYPE
                        The molecular type (column 3 of the GFF file) selected
                        for the BED files, default is exon
  --feature_type FEATURE_TYPE, -ft FEATURE_TYPE
                        The feature type (column 9 of the GFF file) selected
                        for the BED files, default is Parent
  --id_as_features, -id
                        Will set the ID of each element as a string containing
                        all its features
  --path PATH, -p PATH  The location where BED files will be created, default
                        is current working directory
  --name NAME, -n NAME  The name of the BED files, default is the GFF file
                        name
  --verbose, -v         Will outpout in stdout the command arguments and the
                        name of each element raising a warning in consistency
                        check
  --discard, -d         Will discard the element raising a warning in strand
                        consistency and overlapping check
  --skip_exon_number, -s
                        If set, the program will skip addiing _# for exon
                        number.
```


## bed_to_gff.py 

Reverse complement of the gff_to_bed script. This script takes a BED file and converts it into a GFF formatted file. Works on BED12 and BED6.

```
usage: bed_to_gff.py [-h] --bed_file BED_FILE --source SOURCE --mol_type
                     MOL_TYPE [--is_bed12] [--make_gff3]

example run : ./bed_to_gff.py -f tiny_dmel_sample_r5-57.genes.bed6 -m exon -s dmel --make_gff3
		 	  ./bed_to_gff.py -f tiny_dmel_sample_r5-57.genes.bed12 --is_bed12 -m exon -s dmel --make_gff3

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
