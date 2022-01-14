import gzip
import csv
import re
import logging
import sys

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')

def split_metadata_fields(value: str) -> {}:
    metadata_dict = {}
    for entry in value.split(';'):
        metadata_dict.update({entry.split('=')[0]: entry.split('=')[1]})
    return(metadata_dict)

def convert_bed_to_tsv(bedfile_path: str, tsvfile_path: str):
    header_written = False

    logging.debug("Commencing conversion")

    with gzip.open(tsvfile_path, mode='wt') as tsvfile:
        with gzip.open(bedfile_path, mode='rt') as bedfile:
            bedreader = csv.reader(bedfile, delimiter="\t")

            # Skip the first line
            next(bedreader)

            for row in bedreader:
                if not header_written:
                    tsvfile.write("{chrom}\t{start}\t{end}\t{strand}\t{metadata}\n".format(
                        chrom='CHROM',
                        start='START',
                        end='END',
                        strand='STRAND',
                        metadata="\t".join(split_metadata_fields(row[3]).keys())
                    ))
                    header_written = True

                entry = "{chrom}\t{start}\t{end}\t{strand}\t{metadata}\n".format(
                    chrom=row[0],
                    start=int(row[1]) + 1,
                    end=row[2],
                    strand=row[5],
                    metadata="\t".join(split_metadata_fields(row[3]).values())
                )

                tsvfile.write(entry)

                if bedreader.line_num % 10000 == 0:
                    logging.debug("Written {} lines".format(bedreader.line_num))

    logging.debug("Finishing writing")

if __name__ == "__main__":
    convert_bed_to_tsv(sys.argv[1], sys.argv[2])
