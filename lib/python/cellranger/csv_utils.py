import csv
from collections import Counter
import itertools
import os

class CSVParseException(Exception):
    pass

# Unicode BOM
BOM ="\xef\xbb\xbf"

def load_csv_filter_comments(filename, descriptive_name, required_cols):
    """Returns each non-comment line as an element in a list, also verifies ASCII encoding and no
    duplicate columns."""
    if not os.path.isfile(filename):
        raise CSVParseException("Could not find the {} csv file {}".format(descriptive_name, filename))

    if not os.access(filename, os.R_OK):
        msg = "The {} csv is not readable, please check file permissions: {}".format(descriptive_name, filename)
        raise CSVParseException(msg)

    with open(filename, 'rU') as f:
        rows = f.readlines()

    if len(rows) == 0:
        raise CSVParseException("The {} csv file {} has no data.".format(descriptive_name, filename))

    # Remove a BOM if present
    rows[0] = rows[0].replace(BOM, "")
    try:
        rows = map(lambda string: string.decode('ascii'), rows)
    except UnicodeDecodeError:
        raise CSVParseException("The {} csv file {} contains non-ascii characters."
                                  .format(descriptive_name, filename))
    rows = list(itertools.ifilter(lambda x: not x.startswith('#'), rows))

    # Remove whitespace around column headers
    rows[0] = ",".join([x.strip() for x in rows[0].split(",")])
    reader = csv.DictReader(rows)
    col_names = reader.fieldnames
    if not set(col_names).issuperset(set(required_cols)):
        raise CSVParseException(
            'The {} file header must contain the following comma-delimited fields: "{}".'.format(
                descriptive_name,
                ', '.join(required_cols)))
    col_counts = Counter(col_names)
    for (k, v) in col_counts.items():
        if v > 1:
            msg = "{} csv has a duplicated column: {}".format(descriptive_name, k)
            raise CSVParseException(msg)
    return reader

def write_filtered_barcodes(out_csv, bcs_per_genome):
    """ Args:
        bcs_per_genome (dict of str to list): Map each genome to its cell-associated barcodes
    """
    with open(out_csv, 'w') as f:
        writer = csv.writer(f, lineterminator='\n')
        for (genome, bcs) in bcs_per_genome.iteritems():
            for bc in bcs:
                writer.writerow([genome, bc])
