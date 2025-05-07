"""
Created 2021-12-17, version 1

Tools to view and/or load data from HDF5 files created by the motorscan
GUI. This file is a smaller and vamped up version of myH5tools.py

Method(s):

H5()    H5 is a general version
h5()    h5 is more specified for the motorscan gui

"""

import h5py as h5
import numpy as np
from datetime import datetime as DT
from matplotlib import pyplot as plt



# =======================================================================================


def H5(file='',item1='', item2='', item3='', item4=''):
    """

    A simple h5 viewer / browser / what ever

    Example:

        H5('filename.h5')

        prints a list of all entries in filename.h5

        H5('filename.h5', 'entry1')

        prints a list of all items in entry1 and shows if they are datasets, groups, or what ever

        H5('filename.h5', 'entry1', 'measurement')

        prints the datasets, groups and what ever in measurement

        H5('filename.h5', 'entry1', 'measurement', 'mono_energy')

        returns the data in the dataset mono_energy, eg. energy = H5('filename.h5', 'entry1', 'measurement', 'mono_energy')

        H5('filename.h5', 'entry1', 'measurement', 'pre_scan_snapshot')

        prints the content of that group.

        etc., etc., etc.

    """

    if not( type(file) is str and  
            type(item1) is str and 
            type(item2) is str and 
            type(item3) is str and 
            type(item4) is str ):
                print("\nH5(): all arguments passed to H5() must be strings.")
                if type(file) is str:
                    if file == '':
                        print("\nH5(): a good start is to pass a file name. pass (at least) attribute 'file'.)")
                        return None

    #  Check if the file can be opened
    try:
        h5f = h5.File(file, 'r')
    except:
        print("\nH5(): could not find or open the file {0}".format(file))
        return None
    
    if not (h5f.__class__.__name__ == 'File'): 
        print("\nH5(): not sure this is a h5 file but we'll continue anyway...")

    # Check if the first key exists (argument item1)
    key1exists = False
    if item1 == '':
        try:
            for k in h5f.items():
                print(k[0])
            return
        except:
            print("\nH5(): could not find any items in the file {0}".format(file))
            return
    else:
        try:
            for k in h5f.items():
                if k[0] == item1:
                    key1exists = True
                    break
            if not key1exists:
                print('\nH5(): could not find {0}'.format(item1) )
                return
        except:
            print('\nH5(): could not find {0}'.format(item1) )
            return

    Key1 = h5f.get(item1)

    key2exists = False; key2type = ''
    if item2 == '':
        print(item1)
        try:
            for k in Key1.items():
                Name = k[0]
                TypeN = k[1].__class__.__name__
                print("\t{0: <12} {1}".format(TypeN ,Name))
            return
        except:
            print('\nH5(): could not find anything in {0}'.format(item2) )
            return
    else:
        try:
            for k in Key1.items():
                if k[0] == item2:
                    key2exists = True
                    key2type = k[1].__class__.__name__
                    break
            if not key2exists:
                print('\nH5(): could not find ' + item2 + ' in ' + item1)
                return None
        except:
            print('\nH5(): could not find ' + item2 + ' in ' + item1)
            return

    Key2 = Key1.get(item2)

    if key2type == 'Dataset':
        data = h5f[item1+'/'+item2]
        return data[()]

    key3exists = False; key3type = ''
    if item3 == '':
        print(item1)
        print('\t'+item2)
        try:
            for k in Key2.items():
                Name = k[0]
                TypeN = k[1].__class__.__name__
                print("\t\t{0: <12} {1}".format(TypeN ,Name))
            return
        except:
            print('\nH5(): could not find anything in {0} in {1}'.format(item2,item1) )
            return
    else:
        try:
            for k in Key2.items():
                if k[0] == item3:
                    key3exists = True
                    key3type = k[1].__class__.__name__
                    break
            if not key3exists:
                print('\nH5(): could not find ' + item3 + ' in ' + item2 + ' in ' + item1)
                return 
        except:
            return
    Key3 = Key2.get(item3)

    if key3type == 'Dataset':
        data = h5f[item1+'/'+item2+'/'+item3]
        return data[()]
            

    key4exists = False; key4type = ''
    if item4 == '':
        print(item1)
        print('\t'+item2)
        print('\t\t'+item3)
        try:
            for k in Key3.items():
                Name = k[0]
                TypeN = k[1].__class__.__name__
                print("\t\t\t{0: <12} {1}".format(TypeN ,Name))
            return
        except:
            print('\nH5(): could not find anything in ' + item3 )
            return
    else:
        try:
            for k in Key3.items():
                if k[0] == item4:
                    key4exists = True
                    key4type = k[1].__class__.__name__
                    break
            if not key4exists:
                print('\nh5browser():\tcould not find ' + item4 + ' in ' + item3 + ' in ' + item2 + ' in ' + item1)
                return None
        except:
            return
    Key4 = Key3.get(item4)

    if key4type == 'Dataset':
        data = h5f[item1+'/'+item2+'/'+item3+'/'+item4]
        return data[()]
    else:
        print('\nh5browser():')
        print('sorry, but this is as far as this script goes.')
        print('for datasets further down, consult a propper')
        print('h5-file viewer / browser or write your own script.')
        return





# =======================================================================================



def H5ms(file_name='', entry = -1, measurements = [], prescans = []):
    """


    """

    if not( type(file_name) is str and type(entry) is int and type(measurements) is list and type(prescans) is list):  
        print("\nArguments:")
        print("\tfilename       (str)")
        print("\tentry          (int)")
        print("\tmeasurements   (list of strings)")
        print("\tprescans       (list of strings)")
        return

    #  Check if the file can be opened / found
    try:
        h5f = h5.File(file_name, 'r')
    except:
        print("\nCould not find or open the file '{0}'".format(file_name))
        return
    
    if not (h5f.__class__.__name__ == 'File'): 
        print("\nNot sure this is a h5 file but let's try anyway...")
    
    if len(h5f.items()) < 1:
        print("\nCould not find any entries in the file.")
        return
    
    # Make sure that the entry is okay
    key1name = 'entry{0}'.format(entry)
    Key1 = h5f.get(key1name)
    if type(Key1) is type(None):
        print('\nAvailable entries are:')
        tmp = 'entry'
        for item in h5f.items(): tmp += item[0][5:] + ', '
        print(tmp[0:-1])
        return
    
    measurementKeyExists = False
    try:
        for k in Key1.items():
            if k[0] == 'measurement':
                measurementKeyExists = True
    except:
        print("\nCould not find a 'measurement' group in {0}".format(key1name))
        return
    
    if not measurementKeyExists:
        print("\nCould not find a 'measurement' group in {0}".format(key1name))
        return
    
    Key2 = Key1.get('measurement')

    DataSets = []
    Groups = []
    try:
        for k in Key2.items():
            Name = k[0]
            TypeOfItem = k[1].__class__.__name__
            if TypeOfItem == 'Dataset':
                DataSets.append(Name)
            if TypeOfItem == 'Group':
                Groups.append(Name)
    except:
        print("\nCould not find anything in 'measurement'")
        return
    
    if len(DataSets) == 0 and len(Groups) == 0:
        print("\nCould not find anything in 'measurement'")
        return

    if measurements == [] and prescans == []:
        if len(DataSets) > 0:
            print('Datasets:')
            for s in DataSets:
                print('\t{0}'.format(s))
        if len(Groups) > 0:
            print('Groups:')
            for s in Groups:
                print('\t{0}'.format(s))
                if s == 'pre_scan_snapshot':
                    try:
                        tmpk = Key2.get('pre_scan_snapshot')
                        if len(tmpk.items())>0:
                            for item in tmpk.items():
                                print('\t\t{0} ({1})'.format(item[0], item[1].__class__.__name__ ))
                    except:
                        pass    

    if len(measurements) > 0:
        dsok = True
        for dsR in measurements:
            if not dsR in DataSets:
                print("\ndataset '{0}' does not exist in the 'measurement' group".format(dsR))
                dsok = False
    else:
        dsok = False
    
    # deal with datasets in the subgroup here
    Key3 = None
    prsok = True
    if len(prescans) > 0:
        try:
            Key3 = Key2.get('pre_scan_snapshot')
        except:
            print("\nthere is no group 'pre_scan_snapshot' in the 'measurement' group in entry {0}".format(key1name))
            return
        try:
            if len(Key3.items()) == 0:
                print("\nthe group 'pre_scan_snapshot' in the 'measurement' group in entry {0} is empty.".format(key1name))
                return
        except:
            print("\n_unexpected error connected to the group 'pre_scan_snapshot' in the 'measurement' group in entry {0}.".format(key1name))
            return
        try:
            tmp = []
            for item in Key3.items():
                tmp.append(item[0])
            for p in prescans:
                if not p in tmp:
                    print('p',p,'tmp',tmp)
                    prsok = False  
            del tmp  
        except:
            print("\nunexpected error connected to the group 'pre_scan_snapshot' in the 'measurement' group in entry {0}.".format(key1name))
            return

        if not prsok:
            print("available dataset(s) in 'pre_scan_snapshot' is/are:")
            for item in Key3.items():
                if item[1].__class__.__name__ == 'Dataset':
                    print(item[0])


    returnData = []
    if dsok:
        for dsR in measurements:
            for dsE in DataSets:
                if dsR == dsE:
                    returnData.append(Key2.get(dsR)[()])
    if prsok:
        for pssR in prescans:
            for pssE in Key3.items():
                if pssR == pssE[0]:
                    returnData.append(Key3.get(pssR)[()])

    
    return returnData
 
    


    




    
    
    



    



    

    

    

    

    """
    try:
        for k in h5f.items():
            print(k[0])
        return
    except:
        print("\nCould not find any items (groups, datasets,...) in the file {0}".format(file_name))
        return
    
    
    else:
        try:
            for k in h5f.items():
                if k[0] == item1:
                    key1exists = True
                    break
            if not key1exists:
                print('\nH5(): could not find {0}'.format(item1) )
                return
        except:
            print('\nH5(): could not find {0}'.format(item1) )
            return

    Key1 = h5f.get(item1)

    key2exists = False; key2type = ''
    if item2 == '':
        print(item1)
        try:
            for k in Key1.items():
                Name = k[0]
                TypeN = k[1].__class__.__name__
                print("\t{0: <12} {1}".format(TypeN ,Name))
            return
        except:
            print('\nH5(): could not find anything in {0}'.format(item2) )
            return
    else:
        try:
            for k in Key1.items():
                if k[0] == item2:
                    key2exists = True
                    key2type = k[1].__class__.__name__
                    break
            if not key2exists:
                print('\nH5(): could not find ' + item2 + ' in ' + item1)
                return None
        except:
            print('\nH5(): could not find ' + item2 + ' in ' + item1)
            return

    Key2 = Key1.get(item2)

    if key2type == 'Dataset':
        data = h5f[item1+'/'+item2]
        return data[()]

    key3exists = False; key3type = ''
    if item3 == '':
        print(item1)
        print('\t'+item2)
        try:
            for k in Key2.items():
                Name = k[0]
                TypeN = k[1].__class__.__name__
                print("\t\t{0: <12} {1}".format(TypeN ,Name))
            return
        except:
            print('\nH5(): could not find anything in {0} in {1}'.format(item2,item1) )
            return
    else:
        try:
            for k in Key2.items():
                if k[0] == item3:
                    key3exists = True
                    key3type = k[1].__class__.__name__
                    break
            if not key3exists:
                print('\nH5(): could not find ' + item3 + ' in ' + item2 + ' in ' + item1)
                return 
        except:
            return
    Key3 = Key2.get(item3)

    if key3type == 'Dataset':
        data = h5f[item1+'/'+item2+'/'+item3]
        return data[()]
            

    key4exists = False; key4type = ''
    if item4 == '':
        print(item1)
        print('\t'+item2)
        print('\t\t'+item3)
        try:
            for k in Key3.items():
                Name = k[0]
                TypeN = k[1].__class__.__name__
                print("\t\t\t{0: <12} {1}".format(TypeN ,Name))
            return
        except:
            print('\nH5(): could not find anything in ' + item3 )
            return
    else:
        try:
            for k in Key3.items():
                if k[0] == item4:
                    key4exists = True
                    key4type = k[1].__class__.__name__
                    break
            if not key4exists:
                print('\nh5browser():\tcould not find ' + item4 + ' in ' + item3 + ' in ' + item2 + ' in ' + item1)
                return None
        except:
            return
    Key4 = Key3.get(item4)

    if key4type == 'Dataset':
        data = h5f[item1+'/'+item2+'/'+item3+'/'+item4]
        return data[()]
    else:
        print('\nh5browser():')
        print('sorry, but this is as far as this script goes.')
        print('for datasets further down, consult a propper')
        print('h5-file viewer / browser or write your own script.')
        return
"""

