from pysam import VariantFile
import argparse
import os.path as path
import pandas as pd
import configparser

def locations_from_file(filename):
    with open(filename) as fin:
        header = next(fin).strip().split(",")
        enum_header = {j: i for (i, j) in enumerate(header)}
        cols = ["chrom", "pos", "ref", "alt", "Variant", "Amino acid"]
        col_indices = [enum_header[i] for i in cols]
        results = list()
        for line in fin:
            split_line = line.strip().split(",")
            results.append(tuple([split_line[i] for i in col_indices]))

    return results

def cdna_from_var_string(s):
    return s.split(" ")[2]

def get_variant(chrom, pos, ref, alt, variant, amino, bcf_filename_template):
    contig = f"chr{chrom}"
    bcf_filename = bcf_filename_template.format(chrom=chrom)
    bcf_in = VariantFile(bcf_filename)
    results = list()
    header = bcf_in.header
    for rec in bcf_in.fetch(contig, int(pos)-1, int(pos)):
        if (rec.pos == int(pos) and
                rec.ref == ref and
                rec.alts[0] == alt):
            results.append(rec)
    bcf_in.close()
    return (header, results)

def get_dosages(location_filename, output_filename, bcf_filename_template): 
    dat = locations_from_file(location_filename)
    results = [get_variant(*info,bcf_filename_template) for info in dat]

    samples_in_first_vcf = list(results[0][0].samples)
    sample_orders_are_identical = all([list(x[0].samples)
                                       == samples_in_first_vcf for x in results])
    if sample_orders_are_identical:
        with open(output_filename, "wt") as fout:
            samples_str  =",".join(samples_in_first_vcf)
            header=f"variant,{samples_str}\n"
            fout.write(header)
            for ((_, res), info) in zip(results, dat):
                if len(res) == 1:
                    rec = res[0]
                    if len(rec.alts) > 1:
                        raise NameError("Multi-allelic")
                    name =":".join([rec.chrom, str(rec.pos), rec.ref, rec.alts[0]])
                    gts = [rec.samples[s]["GT"] for s in samples_in_first_vcf]
                    dosages = [str(sum(x)) if (-1 not in x) else "" for x in gts]
                    dosages_str = ",".join(dosages)
                    out = f"{name},{dosages_str}\n"
                    fout.write(out)
                elif len(res) > 1:
                    print(f"Something is off {info}")
                else:
                    print(f"No results were found for {info}")
    else:
        raise NameError("Sample order/ samples are not the same between all the files used can't do merge")


def get_counts(dosages_filename, locations_filename, dictionary_csv, counts_output_filename):
    dosages = pd.read_csv(dosages_filename).set_index("variant")
    nwd_ids = dosages.apply(lambda x: list(dosages.columns[x != 0]), axis=1)
    df = pd.DataFrame( { "nwd_ids": nwd_ids, "counts": nwd_ids.apply(len)})

    dict_df = pd.read_csv(dictionary_csv)[["framid", "NWD_ID"]]
    dict_df = dict_df.loc[~dict_df.NWD_ID.isna()]

    my_dict = dict(zip(dict_df["NWD_ID"],dict_df["framid"]))

    df["framids"] = df.nwd_ids.apply(lambda x: [str(my_dict[y]) for y in x])
    df["nwd_ids"] = df.nwd_ids.apply(";".join)
    df["framids"] = df.framids.apply(";".join)

    popseq = pd.read_csv(locations_filename)
    popseq["variant"] = popseq.apply(lambda x: "chr" + ":".join(map(str,x[["chrom","pos","ref","alt"]])),axis=1)
    result_df = popseq.set_index("variant").join(df, how="left")
    result_df.to_csv(counts_output_filename)



def main():
    parser = argparse.ArgumentParser(description="Get bcf file with variants defined in input file")
    parser.add_argument('input', metavar='FILE', 
                        help="File with variants csv with cols variant,Amino acid,spid,chrom,pos,ref,alt,bcftoolsQuery")
    parser.add_argument('--output_prefix', '-o', metavar='OUTPUT_PREFIX', 
                        help="Output directory")

    parser.add_argument('--config', '-c', metavar='CONFIG_FILE', 
                        help="configuration file", default="generate_ror.config")

 
    # Get args
    args = parser.parse_args()
    locations_filename = args.input
    # Read config file
    config = configparser.ConfigParser()
    config.read(args.config)
    dictionary_csv = config["config"]["dictionary"]
    bcf_filename_template = config["config"]["bcf_template"]

    if not path.exists(dictionary_csv):
        raise NameError("Dictionary file does not exist")
    test_bcf = bcf_filename_template.format(chrom=1)
    if not path.exists(test_bcf):
        raise NameError(f"BCF file for chromosome 1 does not exist: {test_bcf}")


    # Construct output filenames
    locations_filename = args.input
    output_dosages_filename = f"{args.output_prefix}_dosages.csv"
    counts_filename = f"{args.output_prefix}_counts.csv"

    # Get dosage file
    get_dosages(locations_filename, output_dosages_filename,bcf_filename_template)

    # Get counts file 
    get_counts(output_dosages_filename, locations_filename, dictionary_csv, counts_filename)   

