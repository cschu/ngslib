#!/usr/bin/env python
'''
Created on Apr 16, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys

from matplotlib_venn import *
from matplotlib import pyplot as plt
import numpy as np




def main(argv):
    
    fig = plt.figure(figsize=(4, 4))
    v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('A', 'B', 'C'))
    c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
    plt.title("Sample Venn Diagram")
    plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]),
                 xytext=(-70, -70), ha='center', textcoords='offset points',
                 bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', color='gray'))
    fig.savefig('test_venn.png')
    
    pass

if __name__ == '__main__': main(sys.argv[1:])
