import sys
import argparse
import os
script = os.path.basename(sys.argv[0])
nArgs=len(sys.argv)
if nArgs < 2:
    aid=script+' takes YourFile.csv (or YourFile.zip) of rowsXcols numbers\n\t ' + \
        'and produces a UMAP result in YourFile.umap.csv'
else:
    print('script arguments:  '+ " ".join(sys.argv[1:]))
    file, ext = os.path.splitext(sys.argv[1])
    aid=script + ' takes ' + sys.argv[1] + ' and produces a\n\t' + \
        'UMAP result in ' + file + 'umap.csv'
    print(aid)

p=argparse.ArgumentParser(aid + '\n')
p.add_argument("--verbose", help="Increase output verbosity", action="store_true")
p.add_argument('csvOrZipFile', help='The csv file to process or zip file which contains')
p.add_argument('--metric', default='euclidean', help='The UMAP type of distance ',
               choices=['euclidean', 'manhattan', 'cosine', 'cityblock', 'correlation'])
p.add_argument('--firstRow', type=int, default=-2, help='The starting row in csv file')
p.add_argument('--rowCount', type=int, default=-2, help='The # of rows in the csv file')
p.add_argument('--n_neighbors', type=int, default=30, help='The # of rows in the csv file')
p.add_argument('--output_dimensions', type=int, default=2, help='UMAP''s # of output dimensions',
               choices=[2,3])
p.add_argument('--plot2D', type=int, default=0, help='Draw 2D plot')
p.add_argument('--min_dist',default=.2, type=float, help='The UMAP min_dist parameter')
p.add_argument('--random_state',default=None, type=int, help='The seed for the random state')
p.add_argument('--labels',default=None, help='A column of categorical labels used to perform supervised UMAP')
p.add_argument('--useTemplate',default=None, help='The UMAP object to use to transform new data')
p.add_argument('--saveTemplate',default=None, help='The name to which to save the UMAP object after dimension reduction completes')
argv = p.parse_args()


csvFile = argv.csvOrZipFile
if not os.path.exists(csvFile):
    print(csvFile + ' does not exist')
    sys.exit(11)

path=os.path.dirname(csvFile);
baseFile=os.path.basename(csvFile)
file, ext = os.path.splitext(baseFile)
usingZip=False;

if ext.upper() == '.ZIP':
    import zipfile

    usingZip=True
    ext = '.csv'
    with zipfile.ZipFile(csvFile) as zf:
        zf.extractall(path)

elif not ext.upper() == '.CSV':
    print('The file ' + csvFile + ' must end in ".csv" and be comma delimited data')
    sys.exit(12)

argsFile=os.path.join(path, file+'.args'+ext)
if  nArgs < 3 and os.path.exists(argsFile):
    import csv
    print('Reading arguments from '+argsFile)
    with open(argsFile, newline='') as aFile:
        argReader = csv.reader(aFile, skipinitialspace=True)
        args = next(argReader)
        argv = p.parse_args(args)

csvFile=os.path.join(path, file+ext)
if not os.path.exists(csvFile):
    print(csvFile + ' does not exist')
    sys.exit(11)

if argv.saveTemplate is not None:
    newTemplate = argv.saveTemplate
    saving = True
else:
    saving = False

if argv.useTemplate is not None:
    templateFile = argv.useTemplate
    if not os.path.exists(templateFile):
        print(templateFile + ' does not exist')
        sys.exit(11)

    import pickle

    with open(templateFile, 'rb') as f:
        trans = pickle.load(f)

    transforming = True
    print('Guiding with template')
else:
    transforming = False

import time
t1_start = time.perf_counter()

import pandas as pd
D=pd.read_csv(csvFile)
rows, cols = D.shape
print(file + '.csv contains ' + "{:,}".format(rows) +
      ' rows X ' + str(cols) + ' columns....')

if argv.labels is not None:
    target = pd.read_csv(argv.labels)
    supervising = True
else:
    supervising = False

starting=argv.firstRow
if starting<0:
    starting=int(rows*.01)

count=argv.rowCount
if count<2:
    ending=rows-1
else:
    ending=starting + count -1
    if ending >= rows:
        ending=rows-1

if ending <= starting:
    print('Ending row ' + str(ending) + ' is less that starting row ' +str(starting))
    sys.exit(13)


labelsFile=os.path.join(path, file+'.labels'+ext)
dataSlice=D.loc[starting:ending, :]

if supervising:
    target_slice = target.loc[starting:ending, :]
    rows, cols = target_slice.shape
    strRows="{:,}".format(rows)
    strCols=str(cols)
    print('Supervising labels ' +strRows+ ' rows X ' + strCols + ' columns....')


rows, cols = dataSlice.shape
strRows="{:,}".format(rows)
strCols=str(cols)
print('Computing UMAP for ' +strRows+ ' rows X ' + strCols + ' columns....')
if argv.random_state == None:
    strRandom='None'
else:
    strRandom=str(argv.random_state)
print('\tmetric='+argv.metric+', min_dist='+str(argv.min_dist)+
      ', n_neighbors='+str(argv.n_neighbors) +',\n\t output_dimensions='+
      str(argv.output_dimensions)+', random='+strRandom + '\n')

import umap
if transforming:
    embedding = trans.transform(dataSlice)
else:
    reducer = umap.UMAP(metric=argv.metric, n_components=argv.output_dimensions,
                    n_neighbors=argv.n_neighbors, min_dist=argv.min_dist,
                    verbose=argv.verbose, random_state=argv.random_state)

    if supervising:
        umap = reducer.fit(dataSlice, y=target_slice)
        embedding = umap.embedding_
    else:
        umap = reducer.fit(dataSlice)
        embedding = umap.embedding_

    if saving:
        newTemplateFile = os.path.join(newTemplate + '.umap')
        import pickle
        with open(newTemplateFile, 'wb') as f:
            pickle.dump(umap, f)

embedRows, embedCols = embedding.shape
if embedCols == 2:
    XX=pd.DataFrame(data=embedding, columns=['Y umap', 'X umap'])
else:
    XX = pd.DataFrame(data=embedding, columns=['Y umap', 'X umap', 'Z umap'])

outFile=os.path.join(path, file+'.umap.csv')
XX.to_csv(outFile, index=False)
t1_stop = time.perf_counter()
print("Elapsed time: %.1f [min]" % ((t1_stop-t1_start)/60))
if embedCols == 2:
    print('Finished... UMAP X/Y of ' + strRows + ' rows is saved\r\n\t to  ' + outFile)
else:
    print('Finished... UMAP X/Y/Z of ' + strRows + ' rows is saved\r\n\t to  ' + outFile)
if usingZip:
    outZipFile = os.path.join(path, file + '.umap.zip')
    zf = zipfile.ZipFile(outZipFile, mode='w')
    try:
        zf.write(outFile)
    finally:
        print('Storing in zip file ' + outZipFile)
        zf.close()
if argv.plot2D == 1 and embedCols == 2:
    import matplotlib as mpl
    mpl.use('TkAgg')
    import matplotlib.pyplot as plt
    plt.scatter(embedding[:, 0], embedding[:, 1])
    plt.gca().set_aspect('equal', 'datalim')
    plt.title('UMAP 2D projection of ' + strCols + ' dimensions', fontsize=24)
    plt.show()

sys.exit(0)
