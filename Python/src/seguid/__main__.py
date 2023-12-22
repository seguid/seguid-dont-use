from argparse import ArgumentParser
from . import __version__

from seguid.chksum import seguid
from seguid.chksum import slseguid
from seguid.chksum import scseguid
from seguid.chksum import dlseguid
from seguid.chksum import dcseguid

parser = ArgumentParser(prog="python -m seguid", description="seguid: Sequence Globally Unique Identifier (SEGUID) for Nucleotide and Amino-Acid Sequences")
parser.add_argument("--version", action="store_true", help="Show version")
parser.add_argument("--type", type=str, nargs="?", help="Type of checksum to calculate")

args = vars(parser.parse_args())

if args.pop("version"):
    print(__version__)
else:
    type=args.pop("type")
    if type == None:
        type="seguid"
        
    ## Read sequence data from the standard input
    seq=input()
    
    if type == "seguid":
        res=seguid(seq)
    elif type == "slseguid":
        res=slseguid(seq)
    elif type == "scseguid":
        res=scseguid(seq)
    else:
        raise ValueError("Unknown --type='" + type + "'")
    
    print(res)
