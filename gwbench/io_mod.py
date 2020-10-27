import os
import sys

#-----Block and unblock printing-----
def block_print(active=1):
    if active:
        sys.stdout = open(os.devnull, 'w')
        return
    else:
        return

def unblock_print(active=1):
    if active:
        sys.stdout = sys.__stdout__
        return
    else:
        return
