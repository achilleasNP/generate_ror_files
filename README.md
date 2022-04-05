# Get dosages and counts from bcf files for specific variants

## Installation:
1. Clone the repository `git clone https://github.com/achilleasNP/generate_ror_files`

Running
-------
1. Make sure you have poetry available. On our cluster it's availabe as part of the python 3 module e.g: (module load python3/3.8.10)
2. Enter the cloned directory, e.g: `cd generate_ror_files`
3. Install dependencies 'poetry install'
4. Modify the config file `sample_generate_ror.config` file and either copy to "generate_ror.config" or
   provide its location as an optional argument of the command.
5. Run the program using `poetry run generate_ror_files annotated_popseq_files.csv output_dir/prefix`

The annotated_popseq_files.csv files have the format

```
  Variant,Amino acid,spid,chrom,pos,ref,alt,bcftoolsQuery
  LMNA 06 949G>A,p.Glu317Lys,NC_000001.11:156135912:G:A,1,156135913,G,A,chr1:156135913-156135913
  TNNT2 08 236T>A,p.Ile79Asn,NC_000001.11:201365637:A:T,1,201365638,A,T,chr1:201365638-201365638
  MYBPC3 25 2429G>A,p.Arg810His,NC_000011.10:47337563:C:T,11,47337564,C,T,chr11:47337564-47337564
  MYBPC3 24 2374T>C,p.Trp792Arg,NC_000011.10:47337728:A:G,11,47337729,A,G,chr11:47337729-47337729
```
