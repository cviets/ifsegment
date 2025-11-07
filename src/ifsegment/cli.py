import argparse
from .run_cyto_mask import main as cyto_main
from .run_mask import main as mask_main
from .run_quant import main as quant_main

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

    mask_parser = subparsers.add_parser("mask", help="Generate nuclear and cytoplasmic masks", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    mask_parser.add_argument("-i", "--input", type=str, required=True, help="Path to directory with .czi images")
    mask_parser.add_argument("-o", "--output", type=str, required=True, help="Path to directory to save masks")
    mask_parser.add_argument("-n", "--nuclear", type=int, default=0, help="Nuclear channel number (0-indexed)")
    mask_parser.add_argument("-c", "--cytoplasmic", type=int, default=3, help="Cytoplasmic channel number (0-indexed)")
    mask_parser.add_argument("-m", "--mode", type=str, default="max", help="z-projection type (choose from 'max' or 'mean')")

    quant_parser = subparsers.add_parser("quantify", help="Quantify protein fluorescence in nucleus and cytoplasm", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    quant_parser.add_argument("-i", "--input", type=str, required=True, help="Path to CZI images")
    quant_parser.add_argument("-m", "--masks", type=str, required=True, help="Path to mask tiff files")
    quant_parser.add_argument("-o", "--output", type=str, required=True, help="Path to output CSV file to store data")
    quant_parser.add_argument("-c", "--channels", type=int, nargs="+", default=[1, 2], help="Channel numbers to quantify (0-indexed)")
    quant_parser.add_argument("-md", "--mode", type=str, default="max", help="z-projection type (choose from 'max' or 'mean')")

    args = parser.parse_args()

    if args.command == "cyto-mask": 
        cyto_main(args.input, args.output, args.channel)
    elif args.command == "mask":
        mask_main(args.input, args.output, args.nuclear, args.cytoplasmic, args.mode)
    elif args.command == "quantify":
        quant_main(args.input, args.masks, args.output, args.channels, args.mode)