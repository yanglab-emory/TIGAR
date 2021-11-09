import os, sys, argparse

# path = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
# TIGAR_dir = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
# path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
# print(path)



# TIGAR_dir = '/home/rparrish/github/TIGAR/'

TIGAR_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
sys.path.append(TIGAR_dir)
import TIGARutils as tg

parser = argparse.ArgumentParser(description='scrap_description')
parser.add_argument('--bleh', type=str, choices=['yes', 'no', 'ugh'], 
	help='the help text', metavar='METAVAR')
args = parser.parse_args()

# tg.print_tigarutils_path()

args.__dict__
ls(args)

dir(parser)

parser.print_help()


bleh={'y': 0, 'n': 1}