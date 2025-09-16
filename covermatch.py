# using a report from DIA-NN or MSGF+ combined with a proteome FASTA file,
# calculate the approximate coverage for all proteins

import argparse
import polars as pl
import fasta
import logging
logger = logging.getLogger(__name__)

def compute_coverage(entry :fasta.Entry, 
                     table :pl.DataFrame,
                     protname_colname :str = 'Protein.Names',
                     s_colname :str = 'Stripped.Sequence') -> float:
    "Compute coverage in a conservative manner"
    subtab = table.filter(pl.col(protname_colname).str.contains(entry.e_name, literal=True))
    if subtab.shape[0] > 0:
        precurs = subtab.select(pl.col(s_colname))[:,0]
        covered_spans = [] # in [start, end) format
        for precursor in precurs:
            try:
                align = entry.seq.index(precursor)
            except ValueError as ve:
                msg = "Warning: Failed alignment of {} in gene {}, hit: {}".format(precursor, entry.gene_name, ve)
                logger.warning(msg)
                #print(msg)
                continue
            covering_span = (align, align + len(precursor))
            covered_spans.append(covering_span)
        covered_spans.sort(key=lambda pair : pair[0])
        if len(covered_spans) == 1:
            begin, end = covered_spans[0]
            frac = (end - begin) / len(entry.seq)
            return frac
        else:
            for i in range(len(covered_spans)-1):
                left_span = covered_spans[i]
                right_span = covered_spans[i+1]
                # previous span completely envelopes next span
                if left_span[1] >= right_span[1]:
                    # right span adds no additional coverage so it should be zero length
                    covered_spans[i+1] = (right_span[1], right_span[1])
                # previous span partially envelops next span
                elif left_span[1] > right_span[0]: # only care about > because intervals are [include, exclude)
                    # only the non-overlapping part of right span adds additional coverage
                    covered_spans[i+1] = (left_span[1], right_span[1])
                # in all other cases, nothing needs doing
            coverlen = 0
            for begin, end in covered_spans:
                coverlen += end - begin
            frac = coverlen / len(entry.seq)
            return frac
    else:
        return 0

def main():
    
    logging.basicConfig(filename="covermatch.log", level=logging.INFO)
    
    parser = argparse.ArgumentParser(prog="Coverage Matcher",
                                     description="Computes per-protein coverage for a proteome FASTA and a mass spec report")
    parser.add_argument("--fasta", help="specify fasta file to use for sequence information", required=True)
    parser.add_argument("--report", help="specify mass-spec report file", required=True,)
    parser.add_argument("--report_format", help="""Specify format of mass-spec report. 
                                                   Supported formats are DIA-NN parquet reports, 
                                                   and MGSF+ .TSV reports""", required=True, choices=["DIANN", "MSGF+"])
    parser.add_argument("--output", 
                        help="specify output file name. Format will be csv regardless of chosen extension",
                        required=True)
    args = parser.parse_args()
    
    if args.report_format == "DIANN":
        tab = pl.read_parquet(args.report)
    elif args.report_format == "MSGF+":
        repname = args.report
        if repname.endswith(".xlsx"):
            tab = pl.read_excel(repname)
        elif repname.endswith(".tsv"):
            tab = pl.read_csv(repname, separator="\t")
            tab = tab.filter(pl.col("SpecEValue").lt(1e-10))
            tab = tab.filter(pl.col("Peptide").str.contains(r"\d").not_())
            tab = tab.unique(pl.col('Peptide', 'Charge'))
        tab = tab.with_columns(pl.col("Peptide").str.extract(r"\.([A-Z]+)\.").alias("Strip.Peptide"))
    #breakpoint()
        
    
    cover_table = {}
    with open(args.fasta) as fhandle:
        fasta_iter = fasta.iter_entries(fhandle)
        for entry in fasta_iter:
            if args.report_format == "DIANN":
                coverage = compute_coverage(entry, tab)
                cover_table[entry.e_name] = coverage
                print(entry.gene_name, coverage)
            elif args.report_format == "MSGF+":
                coverage = compute_coverage(entry,tab,"Protein","Strip.Peptide")
                cover_table[entry.e_name] = coverage
                print(entry.e_name, coverage)
    
    ctabT = pl.from_dict(cover_table)
    coverframe = ctabT.transpose(include_header=True, header_name="ENTRY.NAME", column_names=["Coverage.Fraction"])
    outname = args.output + ".csv"
    coverframe.write_csv(outname)


if __name__ == '__main__':
    main()
    
    