import sys,os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.image as mpimg
import cPickle as pickle
import tty, termios

def getchar():
    # https://gist.github.com/jasonrdsouza/1901709
    # Returns a single character from standard input, does not support special keys like arrows
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch

#Allow user to run from terminal and parse variables:
if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    
    o.add_option('-f', '--filedir', dest = 'filedir', default='.', help = 'master directory of the .pkl and .png files, no / at the end - this will be the directory of the file ABOVE the cluster of folders labelled by month - output2 in alfaburst case')
    o.add_option('-q', '--queue', dest = 'queue', default='.', help = 'path of the priority queue dataframe file')
    o.add_option('-n', '--number', dest = 'num', default = 1, help = 'how many timeseries you would like to label')
    #o.add_option('-o', '--output', dest = 'output', default = './labelleddf', help = 'path of the output (partially) labelled #dataframe')
    opts, args = o.parse_args(sys.argv[1:])
    
    #Print possible class allocations:
    print 'Classes:'
    print '0 : Interesting, follow up'
    print '1 : Unclipped RFI/Noise'
    print '2 : Wide-band, long-duration in time clipped RFI (2016+)'
    print '3 : Wide-band, short-duration in time clipped RFI (2016+)'
    print '4 : Wide-band, short duration clipped RFI (2015)'
    print '5 : Sharp bandpass transition'
    print '6 : Wide-band, bursty clipped RFI (2015)'
    print '7 : Error in spectra captue or replacement'
    print '8 : Systematic int/float overflow'
    print '9 : Known Pulsar'
               
                 
    #Read in sorted priority queue dataframe:
    priorityq = pd.read_pickle(opts.queue)
    labdf = priorityq
                                  
    for i in range(int(opts.num)):
        #Find datname of highest priority datapoint, read in it's image path, must be one exactly one sub-directory below master dir:
        datname = priorityq.ix[i, 'datfile']
        datname = datname[:-4]
        datname1 = datname[0:6]
        datname2 = datname[8:]
        datname = datname1 + 'fb' + datname2
        #print datname
        
        #Get binfactor and buffer numbers, and add to filestring:
        buffer = priorityq.ix[i,'Buffer']
        binfactor = priorityq.ix[i, 'BinFactor']                      
        
        path = opts.filedir + '/*/' + datname + '.buffer%d.d%d.png'%(buffer, binfactor)
        #path = opts.filedir + '/*/' + datname + '.buffer%d.meta.pkl'%(buffer)


        print path
        imgpath = glob.glob(path)
        
        if len(imgpath) == 1:
            imgpath = imgpath[0]
            print 'file found in ' + imgpath

            #Load up data point image:
            print 'Loading timeseries image...'
            plt.ion()
            fig = plt.figure(frameon=False, figsize=(12,8))
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)

            img = mpimg.imread(imgpath)
            plt.imshow(img)
            plt.pause(0.01)

            #Ask for classification
            print 'Input class number:'
            ch = getchar()

            #If valid class is given, move row to end of priority queue with label in place
            if ch > -0.0001 and ch < 9.0001:
                idx = labdf.index[df['index'] == priorityq[i, 'index']]
                labdf.ix[idx, 'Label'] = ch
                labrow = labdf.iloc[idx]
                labdf = pd.concat([labdf.drop(idx, axis=0), labrow], axis=0, ignore_index=True)            

        elif len(imgpath) == 0:
            print 'file not found'
        else:
            print 'multiple files found with matching dat name'
            print imgpath

    #Save new priority queue and new labelled dataframe:
    labdf = labdf.reset_index(drop=True)
    labdf.to_pickle(opts.queue)
    