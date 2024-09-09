#!/usr/bin/env python3
import argparse
import genemap_utils as gmu

#!/usr/bin/env python3
import argparse
import genemap_utils as gmu

def main():
    parser = argparse.ArgumentParser(description="Map gene symbols to Ensembl gene IDs")
    parser.add_argument('-g', '--genesymbols', nargs='+', required=True, help="Gene symbols, space or comma separated")
    args = parser.parse_args()

    # Handle both comma and space-separated gene symbols
    list_of_genesymbols = []
    for item in args.genesymbols:
        list_of_genesymbols.extend(item.split(','))

    # Strip any whitespace from gene symbols
    list_of_genesymbols = [gene.strip() for gene in list_of_genesymbols]

    # Map gene symbols
    list_of_genesymbols_mapped = gmu.mapGenesymbol('genesymbol_to_ensgeneid', list_of_genesymbols)

    # Print results
    for i in range(len(list_of_genesymbols)):
        print(list_of_genesymbols[i], list_of_genesymbols_mapped[i], sep='\t')

if __name__ == '__main__':
    main()
