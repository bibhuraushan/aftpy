from .aftmap import AFTmap
import argparse


def get_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="input filename")
    parser.add_argument("-s", "--save", help="Save image as png.", action="store_true")

    return parser.parse_args()


def aftview():
    args = get_arg()
    aftmap = AFTmap(args.filename)
    if args.save:
        aftmap.plot(save=True)
    else:
        aftmap.plot(save=False)
