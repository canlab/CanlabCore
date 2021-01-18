try:
    import sys
    import argparse
    import os
    import zipfile
    import csv
    import time
    import pandas
    import numpy
    import numba
    import umap
    print('All imports essential to doUmap.py are found')
    exit(0)
except ImportError as exc:
    print(exc)
    exit(100)
