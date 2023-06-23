import os
import re
import sys
import glob
import json
import shutil
import tarfile
import subprocess
import numpy as np


try:
    import pandas as pd
except:
    print('pandas module could not be loaded.')

try:
    from astroquery.alma import Alma
except:
    print('astroquery module could not be loaded')

try:
    from casatools import image as iatool
except:
    print('ia tool could not be loaded')

try:
    from casatools import quanta as qatool
except:
    print('qa tool could not be loaded')

try:
    from astropy.io import fits
except:
    print('astropy module could not be loaded')

try:
    import pyvo
except:
    print('pyvo module could not be loaded')

try:
    from shapely.geometry import Point
    from shapely.affinity import scale, rotate
except:
    print('shapely module could not be loaded')


def queryArchive(uidList=[], uidListFile='', forceDownload=False, sciOnly=True, mfsOnly=True, PLweblog=False):

    topdirname = os.getcwd()

    if uidListFile != '':
        f = open(uidListFile)
        uidList = f.read().splitlines()
        f.close()

    df = pd.DataFrame()

    for uid in uidList:

        os.chdir(topdirname)

        dirname = uid.replace('/', '_').replace(':', '_')

        if os.path.isdir(dirname) == False:
            os.mkdir(dirname)

        os.chdir(dirname)

        datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(f"https://almascience.eso.org/datalink/sync?ID={uid}")
        data_info = datalink.to_table().to_pandas()

        for i in data_info.index:

            if re.search('^DataLink', data_info['service_def'][i], re.I) != None:

                uid = data_info['service_def'][i].replace('DataLink.', '')
                datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(f"https://almascience.eso.org/datalink/sync?ID={uid}")
                data_info1 = datalink.to_table().to_pandas()

                data_info = pd.concat([data_info, data_info1], ignore_index=True)

        df_dict = {}
        df_dict['image file url'] = []
        df_dict['image filename'] = []
        df_dict['PB file url'] = []
        df_dict['PB filename'] = []

        for url in data_info['access_url']:

            filename = os.path.basename(url)

            filenamePattern = '.pbcor.fits$'
            if mfsOnly == True: filenamePattern = '.mfs.*' + filenamePattern
            if sciOnly == True: filenamePattern = '_sci.*' + filenamePattern

            if re.search(filenamePattern, filename) != None:
                df_dict['image file url'].append(url)
                df_dict['image filename'].append(filename)

        for filename in df_dict['image filename']:

            filename = filename.replace('.pbcor.fits', '.pb.fits.gz')

            url_list = []

            for url in data_info['access_url']:

                if re.search(filename, os.path.basename(url)) != None:
                    url_list.append(url)

            if len(url_list) == 1:
                df_dict['PB file url'].append(url_list[0])
                df_dict['PB filename'].append(os.path.basename(url_list[0]).replace('.gz', ''))
            else:
                df_dict['PB file url'].append('')
                df_dict['PB filename'].append('')

        df1 = pd.DataFrame(df_dict)
        df = pd.concat([df, df1], ignore_index=True)

        for i in range(len(df_dict['image file url'])):

            if os.path.exists(df_dict['image filename'][i]) == False or forceDownload == True:
                os.system('curl -O '+df_dict['image file url'][i])

            if os.path.exists(df_dict['PB filename'][i]) == False or forceDownload == True:

                os.system('curl -O '+df_dict['PB file url'][i])

                filename = os.path.basename(df_dict['PB file url'][i])
                os.system('gunzip '+filename)

        if PLweblog == True:

            url_list = []

            for url in data_info['access_url']:

                if re.search('weblog', os.path.basename(url)) != None:
                    url_list.append(url)

            print(url_list)

            if len(url_list) == 1:

                filename = os.path.basename(url_list[0])

                if os.path.exists(filename) == False or forceDownload == True:
                    os.system('curl -O '+url_list[0])

                tar = tarfile.open(filename)
                weblog_path = tar.getnames()[0].split('/')[0]
                tar.extractall()
                tar.close()

            else:

                print('ERROR: multiple PL weblogs found in archive. Skipping.')

    os.chdir(topdirname)

    return(df)

def readPLweblog(path=''):

    with open(path + '/html/t1-1.html', "r") as f:
        df = pd.read_html(f.read(), thousands='')
    f.close()

    df = df[1]

    for i in df.index:

        if re.search('Pipeline Version', df.iat[i, 0], re.I) != None:
            PLversion = df.iat[i, 1].replace('(documentation)', '').strip()

    with open(path + '/html/t1-4.html', "r") as f:
        df1 = pd.read_html(f.read(), thousands='')
    f.close()

    df1 = df1[0]

    df = pd.DataFrame()
    dff = pd.DataFrame()

    for i in df1.index:

        if re.search('hif_makeimlist', df1.iat[i, 0], re.I) != None:
            df2 = readPLweblogPage(path=path, PLversion=PLversion, stage=df1.iat[i, 0].split('.')[0], template='hif_makeimlist')
            df = pd.concat([df, df2], ignore_index=True)

        if re.search('hif_makeimages', df1.iat[i, 0], re.I) != None:
            df2 = readPLweblogPage(path=path, PLversion=PLversion, stage=df1.iat[i, 0].split('.')[0], template='hif_makeimages')
            dff = pd.concat([dff, df2], ignore_index=True)

    df['image file'] = df['imagename']

    for i in df.index:

        ij = np.where(dff['image file'].str.contains(df['imagename'][i].replace('+', '\+').replace('STAGENUMBER', '[0-9_]+')) == True)[0]
        if len(ij) != 1: continue

        df.at[i, 'image file'] = dff.at[ij[0], 'image file']

    df = df.merge(dff, left_on='image file', right_on='image file')

    df['PL image filename'] = df['image file']

    for i in df.index:

        myregexp = '(' + df['imagename'][i].replace('+', '\+').replace('.sSTAGENUMBER', ')(?:\.s[0-9_]+)?(') + ')(.*)(?:iter[0-9]+)(?:.image)'

        filename = re.findall(myregexp, df['PL image filename'][i])
        if len(filename) != 1: sys.exit('ERROR')

        df.at[i, 'PL image filename'] = ''.join(filename[0]) + 'pbcor.fits'

    return(df)

def readPLweblogPage(path='', stage='', template='', PLversion=''):

    with open(path + '/html/stage'+stage+'/t2-4m_details.html', "r") as f:
        df0 = pd.read_html(f.read(), thousands='')
    f.close()

    df = pd.DataFrame()

    if template == 'hif_makeimlist':

        for i in range(len(df0)):

            if df0[i].empty == True: continue

            if PLversion in ['2021.2.0.128', '42866M (Pipeline-CASA56-P1-B)', '40896 (Pipeline-CASA51-P2-B)', '40896M (Pipeline-CASA51-P2-B)', '42254M (Pipeline-CASA54-P1-B)']:

                if 'imagename' in df0[i].columns:
                    df = df0[i]
                    break

    if template == 'hif_makeimages':

        for i in range(len(df0)):

            if df0[i].empty == True: continue

            if PLversion == '2021.2.0.128':

                if 'image file' in df0[i].iloc[:, 0].values.tolist():

                    df1 = df0[i]

                    df1 = np.split(df1, df1[df1.isnull().all(1)].index)

                    for i in range(len(df1)):

                        df2 = df1[i]
                        if len(df2.index) == 0: continue

                        df2 = df2.dropna(how='all').reset_index(drop=True).drop(1).transpose()
                        df2.iat[0, 0] = 'spw'
                        df2 = df2.rename(columns=df2.iloc[0]).reset_index(drop=True).drop(0).dropna(how='all')

                        df = pd.concat([df, df2], ignore_index=True)

                    df.rename(columns = {'spw':'spwname'}, inplace = True)

                    break

            if PLversion in ['42866M (Pipeline-CASA56-P1-B)', '40896 (Pipeline-CASA51-P2-B)', '40896M (Pipeline-CASA51-P2-B)', '42254M (Pipeline-CASA54-P1-B)']:

                if 'Field' in df0[i].columns:

                    df1 = df0[i]

                    grouped1 = df1.groupby('Field')

                    for fieldname in df1['Field'].unique():

                        df2 = grouped1.get_group(fieldname)
                        grouped2 = df2.groupby('Spw')

                        for spwname in df2['Spw'].unique():

                            df3 = grouped2.get_group(spwname)

                            df4 = df3.iloc[:, 3:5]

                            df4 = df4.transpose()
                            df4 = df4.rename(columns=df4.iloc[0]).reset_index(drop=True).drop(0)

                            df4['Field'] = fieldname
                            df4['Spw'] = spwname

                            df = pd.concat([df, df4], ignore_index=True)

                    df.rename(columns = {'Spw':'spwname'}, inplace = True)

                    break

    return(df)

def extractSources(uidList=[], uidListFile=''):

    if uidListFile != '':
        f = open(uidListFile)
        uidList = f.read().splitlines()
        f.close()

    image_fname = []

    for uid in uidList:

        dirname = uid.replace('/', '_').replace(':', '_')
        fname = glob.glob(dirname+'/*.pbcor.fits')

        for i in fname:
            image_fname.append(i)

    if len(image_fname) != 0:

        f = open('temp_imagelist.txt', 'w')
        for i in image_fname:
            f.write(i+'\n')
        f.close()

        with open('temp_casascript.py', "w") as f:
            f.write("import os\n")
            f.write("import sys\n")
            f.write("import json\n")
            f.write("sys.path.append('"+os.getcwd()+"')\n")
            f.write("import ADPutils as au\n")
            f.write("\n")
            f.write("f = open('temp_imagelist.txt')\n")
            f.write("image_fname = f.read().splitlines()\n")
            f.write("f.close()\n")
            f.write("\n")
            f.write("for i in range(len(image_fname)):\n")
            f.write("    src_dict = au.extractSourcesFromImage(image_fname[i])\n")

        f.close()

        subprocess.run('/usr/local/bin/casa --agg --nologger -c "temp_casascript.py"', shell=True)

def extractSourcesFromImage(img_sky_pbcor, img_pb='', min_snr=6, fit_width=10, nmax_src=1000, zerooff_snr=3, write_to_file=True):

    myia = iatool()
    myqa = qatool()

    if img_pb == '':
        img_pb = img_sky_pbcor.replace('.pbcor.fits', '.pb.fits')
        if os.path.exists(img_pb) == False:
            print('ERROR: no PB file found.')
            return({})

    img_sky_pbuncor = 'temp_pbuncor.image'

    myia.open(img_sky_pbcor)
    myia.pbcor(pbimage=img_pb, outfile=img_sky_pbuncor, mode='multiply', overwrite=True)
    myia.close()

    myia.open(img_sky_pbuncor)

    img_stats = myia.statistics()
    img_peak = img_stats['max'][0]

    img_stats = myia.statistics(algorithm='chauvenet')
    img_rms = img_stats['rms'][0]

    flux_cutoff = img_rms * min_snr / img_peak
    fit_width1 = round(fit_width/2.)

    src_dict = {}

    if flux_cutoff < 1:
        src_dict = myia.findsources(nmax=nmax_src, cutoff=flux_cutoff, width=fit_width1, point=False)

    myia.close()

    src_dict2 = {}

    if len(src_dict) != 0:

        myia.open(img_sky_pbuncor)

        img_shape = myia.shape()

        for i in range(src_dict['nelements']):

            src_ra = src_dict['component'+str(i)]['shape']['direction']['m0']['value']
            src_dec = src_dict['component'+str(i)]['shape']['direction']['m1']['value']

            src_pixelpos = myia.topixel([src_ra, src_dec])

            src_dict['component'+str(i)]['pixelcoords'] = src_pixelpos['numeric']

        myia.close()

        for i in range(src_dict['nelements']):

            if 'type' in src_dict['component'+str(i)]['shape'].keys():
                if src_dict['component'+str(i)]['shape']['type'] == 'Point': continue

            peakX = src_dict['component'+str(i)]['pixelcoords'][0]
            peakY = src_dict['component'+str(i)]['pixelcoords'][1]

            if peakX < 0 or peakX > img_shape[0]-1: continue
            if peakY < 0 or peakY > img_shape[1]-1: continue

            peakI = src_dict['component'+str(i)]['flux']['value'][0]
            bMaj = src_dict['component'+str(i)]['shape']['majoraxis']['value']
            bMajUnit = src_dict['component'+str(i)]['shape']['majoraxis']['unit']
            bMin = src_dict['component'+str(i)]['shape']['minoraxis']['value']
            bMinUnit = src_dict['component'+str(i)]['shape']['minoraxis']['unit']
            bPa = src_dict['component'+str(i)]['shape']['positionangle']['value']
            bPaUnit = src_dict['component'+str(i)]['shape']['positionangle']['unit']

            with open('temp_estimates.txt', 'w') as f:
                f.write(f"{peakI}, {peakX}, {peakY}, {bMaj}{bMajUnit}, {bMin}{bMinUnit}, {bPa}{bPaUnit}\n")
            f.close()

            myia.open(img_sky_pbuncor)
            
            img_shape = myia.shape()

            found = False

            fit_width1 = round(fit_width/2.)
            
            if peakX-fit_width1 > 0 and peakY-fit_width1 > 0 and peakX+fit_width1 < img_shape[0]-1 and peakY+fit_width1 < img_shape[1]-1:
            
                src_box = str(peakX-fit_width1)+','+str(peakY-fit_width1)+','+str(peakX+fit_width1)+','+str(peakY+fit_width1)
                src_dict1 = myia.fitcomponents(box=src_box, rms=img_rms, dooff=True, estimates='temp_estimates.txt', logfile='temp_logfile.txt')

                if src_dict1['converged'][0] == True:

                    fit_width1 = round( 3 * src_dict1['results']['component0']['shape']['majoraxis']['value'] * max(src_dict1['pixelsperarcsec']) / 2.) # convert value of 3 into an input parameter

                    if peakX-fit_width1 > 0 and peakY-fit_width1 > 0 and peakX+fit_width1 < img_shape[0]-1 and peakY+fit_width1 < img_shape[1]-1:

                        src_box = str(peakX-fit_width1)+','+str(peakY-fit_width1)+','+str(peakX+fit_width1)+','+str(peakY+fit_width1)
                        src_dict1 = myia.fitcomponents(box=src_box, rms=img_rms, dooff=True, estimates='temp_estimates.txt', logfile='temp_logfile.txt')
                    
                        found = True

            myia.close()

            if found == False: continue

            myia.open(img_sky_pbcor)
            if abs(src_dict1['zerooff']['value']/img_rms) < zerooff_snr:
                src_dict['component'+str(i)]['fitcomp'] = myia.fitcomponents(box=src_box, rms=img_rms, dooff=False, estimates='temp_estimates.txt', logfile='temp_logfile.txt')
            else:
                src_dict['component'+str(i)]['fitcomp'] = myia.fitcomponents(box=src_box, rms=img_rms, dooff=True, estimates='temp_estimates.txt', logfile='temp_logfile.txt')
            myia.close()

        for i in range(src_dict['nelements']):

            if 'fitcomp' not in src_dict['component'+str(i)].keys(): continue

            src_dict1 = src_dict['component'+str(i)]['fitcomp']

            if src_dict1['results']['nelements'] != 1: continue
            if src_dict1['converged'][0] != True: continue

            src_dict1 = src_dict1['results']['component0']

            peakX0 = src_dict['component'+str(i)]['pixelcoords'][0]
            peakY0 = src_dict['component'+str(i)]['pixelcoords'][1]
            peakX1 = src_dict1['pixelcoords'][0]
            peakY1 = src_dict1['pixelcoords'][1]
            peakDist = np.sqrt((peakX0-peakX1)**2 + (peakY0-peakY1)**2)

            if peakDist > 1: continue

            if 'zerooff' in src_dict['component'+str(i)]['fitcomp']:
                src_dict1['zerooff'] = src_dict['component'+str(i)]['fitcomp']['zerooff']

            src_dict1['pixelsperarcsec'] = src_dict['component'+str(i)]['fitcomp']['pixelsperarcsec']

            src_dict1['image filename'] = img_sky_pbcor
            src_dict1['shape']['direction']['RA'] = myqa.angle(src_dict1['shape']['direction']['m0'], form=['time'], prec=9)
            src_dict1['shape']['direction']['DEC'] = myqa.angle(src_dict1['shape']['direction']['m1'], form=['time'], prec=9)

            src_dict2[i] = src_dict1

    if write_to_file == True:

        if src_dict2 != {}:
            for j in src_dict2.keys():
                src_dict2[j]['pixelcoords'] = src_dict2[j]['pixelcoords'].tolist()
                src_dict2[j]['pixelsperarcsec'] = src_dict2[j]['pixelsperarcsec'].tolist()
                src_dict2[j]['flux']['error'] = src_dict2[j]['flux']['error'].tolist()
                src_dict2[j]['flux']['value'] = src_dict2[j]['flux']['value'].tolist()
                if 'zerooff' in src_dict2[j].keys(): src_dict2[j]['zerooff']['value'] = src_dict2[j]['zerooff']['value'].tolist()

        fname = img_sky_pbcor.replace('.fits', '.sourcelist.json')
        with open(fname, 'w') as outfile:
            json.dump(src_dict2, outfile)
        outfile.close()

    return(src_dict2)

def create_ellipse(center, axes, inclination):

    p = Point(*center)
    c = p.buffer(1)
    ellipse = scale(c, *axes)
    ellipse = rotate(ellipse, inclination)

    return ellipse

def filterSources(src_dict, pointSourcesOnly=True, noBackgroundEmission=True, isolatedSourcesOnly=True):

    src_dict1 = {}

    for i in src_dict.keys():

        if pointSourcesOnly == True and src_dict[i]['ispoint'] != True: continue

        if noBackgroundEmission == True and 'zerooff' in src_dict[i].keys(): continue

        if isolatedSourcesOnly == True:

            peakX0 = src_dict[i]['pixelcoords'][0]
            peakY0 = src_dict[i]['pixelcoords'][1]
            bMaj = src_dict[i]['shape']['majoraxis']['value']
            bMin = src_dict[i]['shape']['minoraxis']['value']
            bPa = src_dict[i]['shape']['positionangle']['value']
            pixelsperarcsec = max(src_dict[i]['pixelsperarcsec'])

            ellipse0 = create_ellipse((peakX0, peakY0), (3*bMaj/pixelsperarcsec, 3*bMin/pixelsperarcsec), 90+bPa)

            found = False

            for j in src_dict.keys():

                if j == i: continue

                peakX0 = src_dict[j]['pixelcoords'][0]
                peakY0 = src_dict[j]['pixelcoords'][1]
                bMaj = src_dict[j]['shape']['majoraxis']['value']
                bMin = src_dict[j]['shape']['minoraxis']['value']
                bPa = src_dict[j]['shape']['positionangle']['value']
                pixelsperarcsec = max(src_dict[j]['pixelsperarcsec'])

                ellipse1 = create_ellipse((peakX0, peakY0), (3*bMaj/pixelsperarcsec, 3*bMin/pixelsperarcsec), 90+bPa)

                if ellipse1.intersects(ellipse0):
                    found = True
                    break

            if found == True: continue

        src_dict1[i] = src_dict[i]

    return(src_dict1)

def generateCatalog(uidList=[], uidListFile='', runFilterSources=True):

    if uidListFile != '':
        f = open(uidListFile)
        uidList = f.read().splitlines()
        f.close()

    sourcelist_fname = []

    for uid in uidList:

        dirname = uid.replace('/', '_').replace(':', '_')
        fname = glob.glob(dirname+'/*.sourcelist.json')

        for i in fname:
            sourcelist_fname.append(i)

    src_dict = {}
    src_dict['RA'] = []
    src_dict['DEC'] = []
    src_dict['flux'] = []
    src_dict['flux error'] = []
    src_dict['image filename'] = []

    for i in sourcelist_fname:
        
        with open(i, 'r') as f:
            src_dict1 = json.load(f)
        f.close()

        if runFilterSources == True:
            src_dict1 = filterSources(src_dict1)

        for i in src_dict1.keys():

            src_dict['RA'].append(src_dict1[i]['shape']['direction']['RA'][0])
            src_dict['DEC'].append(src_dict1[i]['shape']['direction']['DEC'][0])
            src_dict['flux'].append(src_dict1[i]['flux']['value'][0])
            src_dict['flux error'].append(src_dict1[i]['flux']['error'][0])
            src_dict['image filename'].append(os.path.basename(src_dict1[i]['image filename']))

    src_df = pd.DataFrame(src_dict)

    src_csv_data = src_df.to_csv('catalog.csv', index = False)

#     hdr = fits.Header()
#     hdr['OBSERVER'] = 'ALMA'
#     empty_primary = fits.PrimaryHDU(header=hdr)
#
#     df_rec = df.filter(items=['ra', 'dec', 'flux', 'flux_error', 'rms']).to_records(index=False)
#     table_hdu = fits.BinTableHDU(data=df_rec)
#
#     hdul = fits.HDUList([empty_primary, table_hdu])
#
#     hdul.writeto('results.fits')

    return(src_df)

def extractMetadata(uidList=[], uidListFile=''):

    topdirname = os.getcwd()

    if uidListFile != '':
        f = open(uidListFile)
        uidList = f.read().splitlines()
        f.close()

    df = pd.DataFrame()

    for uid in uidList:

        os.chdir(topdirname)

        dirname = uid.replace('/', '_').replace(':', '_')

        if os.path.isdir(dirname) == False:
            print('ERROR: no data found. Skipping.')
            continue

        os.chdir(dirname)

        df_PLweblog = pd.DataFrame()
        PLweblog = glob.glob('pipeline-*')
        if len(PLweblog) == 1:
            df_PLweblog = readPLweblog(path=PLweblog[0])
        else:
            print(PLweblog)

        filename = glob.glob('*.pbcor.fits')

        for i in range(len(filename)):

            with fits.open(filename[i]) as hdul:
                hdr_info = hdul[0].header.cards

            df_dict = {}
            df_dict['image filename'] = filename[i]

            for j in range(len(hdr_info)):
                keyw = hdr_info[j][0]
                if keyw in ['BMAJ', 'BMIN', 'ROBUST']:
                    df_dict[keyw] = [hdr_info[j][1]]

            if 'PL image filename' in df_PLweblog.columns:

                ij = np.where(df_PLweblog['PL image filename'].str.contains(filename[i].replace('member.', '')) == True)[0]
                if len(ij) == 1:

# Index(['field', 'intent', 'spw', 'phasecenter', 'cell', 'imsize', 'imagename',
#        'specmode', 'start', 'width', 'nbin', 'nchan', 'robust', 'nterms',
#        'uvrange', 'image file', 'centre frequency of image', 'beam',
#        'beam p.a.', 'final theoretical sensitivity', 'cleaning threshold',
#        'clean residual peak / scaled MAD', 'non-pbcor image RMS',
#        'pbcor image max / min', 'fractional bandwidth / nterms',
#        'aggregate bandwidth', 'score', 'Field', 'spwname',
#        'centre frequency of cube', 'non-pbcor image RMS / RMSmin / RMSmax',
#        'channels', 'PL image filename'],
#       dtype='object')

                    if 'cell' in df_PLweblog.columns:

                        cellsize = re.findall('[0-9E\.]+ *arcsec', df_PLweblog.at[ij[0], 'cell'], re.IGNORECASE)
                        if len(cellsize) > 0:
                            df_dict['CELL'] = [float(cellsize[0].upper().replace('ARCSEC', ''))]

                    if 'ROBUST' not in df_dict.keys() and 'robust' in df_PLweblog.columns:

                        robustfactor = re.findall('-?[0-2](\.[0-9]+)?', str(df_PLweblog.at[ij[0], 'robust']))
                        if len(robustfactor) > 0:
                            df_dict['ROBUST'] = [float(robustfactor[0])]

                    if 'final theoretical sensitivity' in df_PLweblog.columns:

                        rms0 = re.findall('[0-9E\.\-]+ *m?jy/beam', df_PLweblog.at[ij[0], 'final theoretical sensitivity'], re.IGNORECASE)
                        if len(rms0) > 0:
                            rms0 = rms0[0]
                            if re.search('mjy/beam', rms0, re.IGNORECASE) != None:
                                df_dict['RMS0'] = [float(rms0.upper().replace('MJY/BEAM', ''))/1000.]
                            else:
                                if re.search('jy/beam', rms0, re.IGNORECASE) != None:
                                    df_dict['RMS0'] = [float(rms0.upper().replace('JY/BEAM', ''))]

                    if 'non-pbcor image RMS' in df_PLweblog.columns:

                        rms1 = re.findall('[0-9E\.\-]+ *m?jy/beam', df_PLweblog.at[ij[0], 'non-pbcor image RMS'], re.IGNORECASE)
                        if len(rms1) > 0:
                            rms1 = rms1[0]
                            if re.search('mjy/beam', rms1, re.IGNORECASE) != None:
                                df_dict['RMS1'] = [float(rms1.upper().replace('MJY/BEAM', ''))/1000.]
                            else:
                                if re.search('jy/beam', rms1, re.IGNORECASE) != None:
                                    df_dict['RMS1'] = [float(rms1.upper().replace('JY/BEAM', ''))]

            dff = pd.DataFrame(df_dict)
            df = pd.concat([df, dff], ignore_index=True)

    os.chdir(topdirname)

    return(df)
