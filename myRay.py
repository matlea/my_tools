
"""

Import data from Ray-UI csv-files

Updates:
May, 2020:  -   updated fileInfo(), load1p(), load2p(), rawrays(), and mapRawrays()
                to make the file name argument meet current... standards.

Sept 19, 2022   -   Updated and simplified fileInfo() and load1p().
                    They are no longer using the csv module.
                    They are not sensitive to delimeter and the annoying sep=.
                    The previous load1p() is still here, now called load1p_().
                    Also, the data is returned as read, i.e. not transposed
                    as it used to be.

"""

__version__ = "22.09.19"
__author__  = "Mats Leandersson"

import csv
import numpy as np






# ---------------------------------------------------------

def fileInfo(fn = ''):
    """
    Enumerate and print column names.
    """

    try:
        f = open(fn, 'r')
    except:
        print('fileInfo(): could not find or open file "{0}".'.format(fn))
        return
    
    header = ''
    # data = []
    for i, row in enumerate(f):
        if row.startswith('sep='):
            pass
        elif row.startswith('P0'):
            header = row.split('\t')
            if len(header) == 1:
                header = row.split(sep = ';')
                break
        else:
            pass
    f.close()

    if len(header) == 0:
        print('fileInfo(): could not find column descriptions in file "{0}".'.format(fn))
        return
    else:
        print('Columns in "{0}":'.format(fn))
        for i, h in enumerate(header):
            print('{0:<3}\t{1}'.format(i,h))





# ---------------------------------------------------------

def load1p(fn = '', col = [], shup = True):
    """
    Load Ray loop data with one scan parameter.
    col = [] loads all data, while e.g. col = [1, 3, 9] only loads column 1 (usually the scan parameter) and columns 3 and 9.
    shup is not used.
    """
    try:
        f = open(fn, 'r')
    except:
        print('load1p(): could not find or open file "{0}".'.format(fn))
        return
    
    header = ''
    data = []
    for i, row in enumerate(f):
        if row.startswith('sep='):
            pass
        elif row.startswith('P0'):
            header = row.split('\t')
            if len(header) == 1:
                header = row.split(sep = ';')
        else:
            r = row.split(sep = '\t')
            if len(r) == 1:
                r = row.split(sep = ';')
            r[-1] = r[-1].replace('\n','')
            for i, rr in enumerate(r):
                r[i] = np.float64(rr)
            data.append(r)
    f.close()


    if len(col) > 0:

        data_ = []
        header_ = []
        for c in col:
            if c < len(header):
                data_.append(np.transpose(data)[c])
                header_.append(header[c])
            else:
                print('load1p(): column {0} does not exist.'.format(c))
        data = np.transpose(data_)
        header = list(header_)

    return data, header
    


# ---------------------------------------------------------


def load2p(fn = '', numP = [10, 13], xyCol = [2, 3], zCol = [32, 33, 34, 35], shup = False, **kwargs):
    """
    For loading data from Ray-UI loop-data csv 
    files with two scan parameters.

    numpP is an array with the number of points
    for each of the scan parameters. xyCol is an 
    array with the column numbers for the scan 
    parameters. zCol is an array with the column 
    numbers of the columns to load.

    Returns an array for scan parameter 1, an 
    array for scan parameter 2, and an array with 
    2d-arrays for the outputs.

    """

    #f len(kwargs) > 0: 
    #    for kw in ['filename', 'fileName', 'file_name']:
    #        _kw = kwargs.get(kw, '')
    #        if type(_kw) is str:
    #            if not _kw == '':
    #                fn = _kw
    #                break

    fields = []
    rows = []
    try:
        with open(fn, 'r') as csvfile: 
            try:
                csvreader = csv.reader(csvfile) 
                tmp = next(csvreader)
                fields = next(csvreader)
                fields = list(fields[0].split("\t"))
            except:
                if not shup: 
                    print("\n\tError. Could not open or read '{0}'".format(fn))
                return False
            try:
                for row in csvreader: 
                    rows.append(row)
            except:
                if not shup: 
                    print("\n\tAn unexpected error when reading file '{0}'".format(fn))
                return False
        csvfile.close
    except:
        if not shup: 
            print("\n\tError. Could not open or read '{0}'".format(fn))
        return False
        
    size1 = 0; size2 = 0
    try:
        size1, size2 = np.shape(rows)
    except:
        if not shup: 
            print("\n\tAn unexpected error when reading file '{0}'".format(fn))
        return False
    if size1<1 or size2<1:
        if not shup:
            print("\n\tAn unexpected error when reading file '{0}'".format(fn))
        return False
    allData = []
    try:
        with open(fn, 'r') as csvfile: 
            csvreader = csv.reader(csvfile) 
            tmp = next(csvreader)
            tmp = next(csvreader)
            for row in csvreader: 
                r = list(row[0].split("\t"))
                allData.append(r)
            csvfile.close
    except:
        if not shup: 
            print("\n\tAn unexpected error when reading file '{0}'".format(fn))
        return False
    allData = np.array(allData)
    data = []
    for i in zCol:
        data.append(allData[:,i])
    data = np.array(data)
    yData = allData[0:numP[1], xyCol[1]]
    yData = [float(numeric_string) for numeric_string in yData]
    yData = np.array(yData)
    xData = []
    for i in range(0, numP[0]):
        xData.append(float(allData[i*numP[1], xyCol[0]]))
    xData = np.array(xData)
    xData = np.asfarray(xData, float)
    yData = np.asfarray(yData, float)
    data = np.asfarray(data, float)
    dataArray_all = []
    for h in range(len(zCol)):
        dataArray = np.zeros(((numP[0],numP[1])))
        for i in range(numP[0]):
            for j in range(numP[1]):
                dataArray[i][j] = data[ h ][ i*numP[1] + j]
        dataArray_all.append(dataArray)
    return xData, yData, dataArray_all
    





# ---------------------------------------------------------

def rawrays(rayfilename = '', columns = [3,4], **kwargs):
    """
    For Ray-UI raw data and footprint csv files.

    For raw files use columns = [3,4], [3,5], or [4,5] depending on orientation.

    For footprint files use columns = [2,3]

    Columns in a Ray raw file:
    3: ox   4: oy   5: oz   6: dx   7: dy   8: dz

    Columns in Ray footprint file:
    2: ox   3: oy

    """

    #if len(kwargs) > 0:
    #    for kw in ['filename', 'fileName', 'file_name']:
    #        _kw = kwargs.get(kw, '')
    #        if type(_kw) is str:
    #            if not _kw == '':
    #                fn = _kw
    #                break

    cols = columns
    fields = []
    rows = []
    Y = []
    try:
        with open(rayfilename, 'r') as csvfile: 
            csvreader = csv.reader(csvfile) 
            tmp = next(csvreader)
            fields = next(csvreader)
            fields = list(fields[0].split("\t"))
            for row in csvreader: 
                rows.append(row)
            csvfile.close
    except:
        print('\n\tError when trying to open or read the file {0}'.format(rayfilename))
        return Y
    try:
        for i, row in enumerate(rows):
            R = list(row[0].split("\t")) 
            y = []
            for c in cols:
                y.append(float(R[c]))
            Y.append(y)
        Y = np.array(Y)
    except:
        print('Unexpected error when trying to extract data from file {0}'.format(rayfilename))
    return Y







# ---------------------------------------------------------

def mapRawrays(rawrays = [], xbins = 10, ybins = 10, xlims = [None, None], ylims=[None, None], transpose = True):
    """
    Make a map out of data loaded with rawrays().
    Returns xaxis (1d), yaxis (1d), and intensity map (2d).
    """
    rays = np.array(rawrays)
    imap = np.zeros([xbins+1,ybins+1])
    if type(xlims[0]) is type(None): minX = np.min(rays[:,0])
    else: minX = xlims[0]
    if type(xlims[1]) is type(None): maxX = np.max(rays[:,0])
    else: maxX = xlims[1]
    if type(ylims[0]) is type(None): minY = np.min(rays[:,1])
    else: minY = ylims[0]
    if type(ylims[1]) is type(None): maxY = np.max(rays[:,1])
    else: maxY = ylims[1]

    X = np.linspace(minX,maxX,xbins+1)
    Y = np.linspace(minY,maxY,ybins+1)
    
    dx = X[1]-X[0]; dy = Y[1]-Y[0]
    for ix, x in enumerate(X):
        for iy, y in enumerate(Y): 
            tmp = rays[ np.all( [ rays[:,0] > x, rays[:,0] <= x+dx, rays[:,1] > y, rays[:,1] <= y+dy  ], axis = 0)]
            imap[ix,iy] += len(tmp)
    
    if transpose: imap = imap.transpose()
    return X, Y, imap

def mapRawRays(rawrays = [], xbins = 10, ybins = 10, xlims = [None,None], ylims=[None,None], transpose = True):
    """ A dummy for mapRawrays() as I don't remember how I name methods..."""
    return mapRawrays(rawrays = rawrays, xbins = xbins, ybins = ybins, xlims = xlims, ylims = ylims, transpose = transpose)


# ---------------------------------------------------------










# --------------------------------------------------------- Old

def load1p_(fn = '', col = [1, 30, 31, 32, 33], shup = False, **kwargs):
    """

    Old load1p().

    >> For Ray-UI loop data csv files <<

    For loading data from Ray-UI loop-data csv files
    with one scan parameter.

    Returns the data from the specified columns as
    an array. The first column should contain values
    for the scan parameter, the next columns should 
    contain Ray outputs.

    """

    fields = []
    rows = []
    try:
        with open(fn, 'r') as csvfile: 
            csvreader = csv.reader(csvfile) 
            tmp = next(csvreader)
            fields = next(csvreader)
            fields = list(fields[0].split("\t"))
            for row in csvreader: 
                rows.append(row)
            csvfile.close
    except:
        print("\n\tError. Could not open or read '{0}'".format(fn))
        return None

    size1 = 0; size2 = 0
    try:
        size1, size2 = np.shape(rows)
    except:
        if not shup: print("\n\tUnexpected error with file '{0}'".format(fn))
        return None
    if size1<1 or size2<1:
        if not shup: print("\n\tUnexpected error with file '{0}'".format(fn))
        return None

    try:
        Y = []
        for row in rows:
            R = list(row[0].split("\t")) 
            y = []
            for c in col:
                y.append(float(R[c]))
            Y.append(y)
        Y = np.transpose(Y)
    except:
        if not shup: print("\n\tAn unexpected error when reading file '{0}'".format(fn))
        return None
    return Y