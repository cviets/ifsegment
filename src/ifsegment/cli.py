import argparse
from .run_cyto_mask import main as cyto_main

def main():
    parser = argparse.ArgumentParser(
        prog="ifsegment", 
        description="Cell segmentation from IF images", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    subparsers = parser.add_subparsers(dest="command", required=False)

    cyto_parser = subparsers.add_parser("cyto-mask", help="Generate cytoplasmic masks", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cyto_parser.add_argument("-i", "--input", type=str, required=True, help="Path to directory with .czi images")
    cyto_parser.add_argument("-o", "--output", type=str, required=True, help="Path to directory to save masks")
    cyto_parser.add_argument("-c", "--channel", type=int, required=True, help="Channel number of cytoplasmic marker (0-indexed)")

    args = parser.parse_args()

    if args.command == "cyto-mask":
        cyto_main(args.input, args.output, args.channel)