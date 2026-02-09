# using a report from DIA-NN or MSGF+ combined with a proteome FASTA file,
# calculate the approximate coverage for all proteins

import argparse
import polars as pl
import fasta
import logging
logger = logging.getLogger(__name__)


def compute_union_coverage(entry :fasta.Entry, dda_tables :list[pl.DataFrame], dia_tables :list[pl.DataFrame]):
    "compute union coverage of DIA and DIA data"
    subtables = []
    all_zerocount = True
    for dda_table in dda_tables:
        subtab = dda_table.filter(pl.col("Protein").str.contains(entry.e_name, literal=True))
        if len(subtab) > 0:
            subtab = subtab.rename({"Protein":"Protein.Names", "Strip.Peptide":"Stripped.Sequence"})
            subtables.append(subtab)
            all_zerocount = False
    for dia_table in dia_tables:
        subtab = dia_table.filter(pl.col("Protein.Names").str.contains(entry.e_name, literal=True))
        if len(subtab) > 0:
            subtables.append(subtab)
            all_zerocount = False
    if all_zerocount:
        return 0
    precursor_lists = [s_tab.select(pl.col("Stripped.Sequence"))[:,0] for s_tab in subtables]
    precursors = []
    for p_list in precursor_lists:
        precursors.extend(p_list)
    
    covered_spans = []
    for precursor in precursors:
        try:
            align = entry.seq.index(precursor)
        except ValueError as ve:
            msg = "Warning: Failed alignment of {} in gene {}, hit: {}".format(precursor, entry.gene_name, ve)
            logger.warning(msg)
            continue
        covering_span = (align, align+len(precursor))
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
            if left_span[1] >= right_span[1]:
                covered_spans[i+1] = (right_span[1],right_span[1])
            elif left_span[1] > right_span[0]:
                covered_spans[i+1] = (left_span[1], right_span[1])
        coverlen = 0
        for begin, end in covered_spans:
            coverlen += end - begin
        frac = coverlen / len(entry.seq)
        return frac


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
                                     description="Computes per-protein coverage for a proteome FASTA and a mass spec report",
                                     epilog="Note that format should be specified for each report file")
    parser.add_argument("--fasta", help="specify fasta file to use for sequence information", required=True)
    
    parser.add_argument("--reports", help="specify mass-spec report files", nargs="+", required=True, action="extend")
    
    parser.add_argument("--report_formats",
                        help="""specify formats of mass-spec reports in corresponding order.
                        First argument to this controls expected format for the first report file, and so on""",
                        nargs="+",
                        required=True,
                        action="extend",
                        choices=["DIANN", "MSGF+"]
                        )
    parser.add_argument("--output", 
                        help="specify output file name. Format will be csv regardless of chosen extension",
                        required=True)
    args = parser.parse_args()
    
    if len(args.reports) != len(args.report_formats):
        raise Exception("discrepancy exists: {} reports but {} formats specified. Quitting".format(len(args.reports), len(args.report_formats)))
    
    
    diann_tables = []
    msgf_plus_tables = []
    
    for i, fmt in enumerate(args.report_formats):
        report_name = args.reports[i]
        if fmt == 'DIANN':
            tab = pl.read_parquet(report_name)
            tab = tab.with_columns(pl.lit(report_name).alias('report_name'))
            diann_tables.append(tab)
        elif fmt == 'MSGF+':
            if report_name.endswith(".xlsx"):
                tab = pl.read_excel(report_name)
            elif report_name.endswith(".tsv"):
                tab = pl.read_csv(report_name, separator='\t')
                tab = tab.filter(pl.col("SpecEValue").lt(1e-10))
                tab = tab.filter(pl.col("Peptide").str.contains(r"\d").not_())
                tab = tab.unique(pl.col('Peptide', 'Charge'))
            #breakpoint()
            tab = tab.with_columns(
                pl.col('Peptide')
                # extract content between periods
                .str.extract(r"\.(.*)\.")
                # remove mass-added modifications
                .str.replace_all(r"\+\d+\.\d+", "")
                .alias("Strip.Peptide"),
                pl.lit(report_name).alias('report_name')
            )
            #breakpoint()
            #tab = tab.with_columns(pl.col('Peptide').str.extract("r\.([A-Z]+)\.").alias("Strip.Peptide"))
            msgf_plus_tables.append(tab)
    
    #cover_table = {}
    cover_table = {}
    cover_table["entry_name"] = []
    cover_table["gene_name"] = []
    cover_table["coverage"] = []
    cover_table["library_kind"] = []
    cover_table['report_name'] = []
    
    with open(args.fasta) as fhandle:
        fasta_iter = fasta.iter_entries(fhandle)
        for entry in fasta_iter:
            for dia_table in diann_tables:
                coverage = compute_coverage(entry, dia_table)
                cover_table["entry_name"].append(entry.e_name)
                cover_table["gene_name"].append(entry.gene_name)
                cover_table["coverage"].append(coverage)
                cover_table["library_kind"].append("DIA")
                cover_table['report_name'].append(dia_table[0, 'report_name'])
            for dda_table in msgf_plus_tables:
                #breakpoint()
                coverage = compute_coverage(entry, dda_table, "Protein", "Strip.Peptide")
                cover_table["entry_name"].append(entry.e_name)
                cover_table["gene_name"].append(entry.gene_name)
                cover_table["coverage"].append(coverage)
                cover_table["library_kind"].append("DDA")
                cover_table['report_name'].append(dda_table[0, 'report_name'])
                pass
            union_cov = compute_union_coverage(entry, msgf_plus_tables, diann_tables)
            cover_table['entry_name'].append(entry.e_name)
            cover_table['gene_name'].append(entry.gene_name)
            cover_table['coverage'].append(union_cov)
            cover_table['library_kind'].append("Union DIA,DDA")
            cover_table['report_name'].append('N/A')
    cover_tab = pl.from_dict(cover_table, 
                         schema={'entry_name':pl.String,
                                 'gene_name':pl.String,
                                 'coverage':pl.Float64,
                                 'library_kind':pl.String,
                                 'report_name':pl.String,
                                }
                        )
    outname = args.output + ".csv"
    print("writing output to: {}".format(outname))
    cover_tab.write_csv(outname)


if __name__ == '__main__':
    main()
    
    