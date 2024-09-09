#!/usr/bin/env python3
import argparse
import genemap_utils as gmu

def main():
    # Get available data handles
    available_data_handles = gmu.listAvailableData()

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Map gene symbols using available data handles")
    parser.add_argument('-g', '--genesymbols', nargs='+', required=True, help="Gene symbols, space or comma separated")
    parser.add_argument('-dh', '--data_handle', default='genesymbol_to_ensgeneid', help="Data handle for mapping (default: genesymbol_to_ensgeneid)")
    parser.add_argument('-H', '--handles', action='store_true', help="Print available data handles and exit")
    parser.add_argument('-d', '--debug', action='store_true', help="Enable debug output")

    args = parser.parse_args()

    # Debugging output for available data handles
    if args.debug:
        if available_data_handles is None:
            print("Debug: gmu.listAvailableData() returned None.")
        else:
            print("Debug: gmu.listAvailableData() returned:", available_data_handles)

    # If --handles is provided, print available data handles and exit
    if args.handles:
        print("Available data handles:", ", ".join(available_data_handles))
        return

    # Check if provided data_handle is valid
    if args.data_handle not in available_data_handles:
        print(f"Error: '{args.data_handle}' is not a valid data handle.")
        print("Available data handles:", ", ".join(available_data_handles))
        return

    # Handle both comma and space-separated gene symbols
    list_of_genesymbols = []
    for item in args.genesymbols:
        list_of_genesymbols.extend(item.split(','))

    # Strip any whitespace from gene symbols
    list_of_genesymbols = [gene.strip() for gene in list_of_genesymbols]

    # Debugging output for gene symbols
    if args.debug:
        print("Debug: list_of_genesymbols:", list_of_genesymbols)

    # Map gene symbols
    list_of_genesymbols_mapped = gmu.mapGenesymbol(args.data_handle, list_of_genesymbols)

    # Debugging output for gene symbol mappings
    if args.debug:
        if list_of_genesymbols_mapped is None:
            print("Debug: gmu.mapGenesymbol() returned None.")
        else:
            print("Debug: gmu.mapGenesymbol() returned:", list_of_genesymbols_mapped)

    # Replace None or empty mappings with "NA"
    list_of_genesymbols_mapped = ['NA' if not mapping else mapping for mapping in list_of_genesymbols_mapped]

    # Print results
    for i in range(len(list_of_genesymbols)):
        print(list_of_genesymbols[i], list_of_genesymbols_mapped[i], sep='\t')

if __name__ == '__main__':
    main()
