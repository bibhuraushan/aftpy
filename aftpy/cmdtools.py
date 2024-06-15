__all__ = ['get_arg_file', 'get_arg_dir', 'aftview', 'aftconvertor']

import argparse
import os
from .aftmap import AFTmap, AFTmaps


def get_arg_file():
    """
    Parse command line arguments for processing a single AFTmap file.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="Process a single AFTmap file.")
    parser.add_argument("filename",
                        help="Input filename")
    parser.add_argument("-s", "--save",
                        help="Save image as png.",
                        action="store_true")
    parser.add_argument("-sm", "--show-mask",
                        default=True,
                        help="Show mask image.",
                        action="store_true")
    return parser.parse_args()


def get_arg_dir():
    """
    Parse command line arguments for processing multiple AFTmap files in a directory.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="Process multiple AFTmap files in a directory.")
    parser.add_argument("path",
                        help="Input directory path.",
                        type=str,
                        default="./")
    parser.add_argument("-o", "--outdir",
                        type=str,
                        help="Output directory path for saving images.",
                        default="./")
    parser.add_argument("-ft", "--file-type",
                        type=str,
                        help="Input File Type",
                        default="aftmap",
                        choices=["aftmap", "aftold", "hipft"])
    parser.add_argument("-df", "--date-format",
                        type=str,
                        help="Date format for input files",
                        default="AFTmap_%Y%m%d_%H%M.h5")
    parser.add_argument("-cto", "--convert-to",
                        type=str,
                        help="Convert images to specified format",
                        default="png",
                        choices=["png", "fits"])
    parser.add_argument("-v", "--verbose",
                        help="Show conversion progress.",
                        action="store_true")
    parser.add_argument("-sm", "--show-mask",
                        default=True,
                        help="Show mask image.",
                        action="store_true")
    return parser.parse_args()


def aftview():
    """
    View or save an AFTmap image based on command line arguments.

    This function reads a filename from command line arguments, creates an AFTmap object,
    and either displays or saves the map based on the provided options.
    """
    args = get_arg_file()
    aftmap = AFTmap(args.filename)
    if args.save:
        aftmap.plot(save=True, show_mask=args.show_mask)
    else:
        aftmap.plot(save=False, show_mask=args.show_mask)


def aftconvertor():
    """
    Convert multiple AFTmap files in a directory to a specified format.

    This function reads directory and conversion options from command line arguments,
    creates an AFTmaps object, and converts all files in the directory to the specified format.
    """
    try:
        args = get_arg_dir()
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir, exist_ok=True)
        aftmaps = AFTmaps(path=args.path, filetype=args.file_type, date_fmt=args.date_format, verbose=args.verbose)
        if args.outdir == "./":
            os.makedirs(f"./aft_{args.convert_to}", exist_ok=True)
            outdir = os.path.abspath(f"./aft_{args.convert_to}")
        else:
            outdir = os.path.abspath(args.outdir)
        print(f"Output saved in: {os.path.abspath(outdir)}")
        aftmaps.convert_all(convert_to=args.convert_to, outpath=outdir, show_mask=args.show_mask)
    except OSError:
        print("Command returned with Exit Status 1")


