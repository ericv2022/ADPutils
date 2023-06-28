import os
import re
import sys
import glob
import json
import shutil
import tarfile
import subprocess
import numpy as np
from multiprocessing.pool import ThreadPool

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

try:
    import requests
except:
    print('requests module could not be loaded')




def downloadImagesPerMOUS(MOUSuid, mfsOnly = True, sciOnly = True, PLweblog = True):

    datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(f"https://almascience.eso.org/datalink/sync?ID={MOUSuid}")
    data_info = datalink.to_table().to_pandas()

    for i in data_info.index:

        if re.search('^DataLink', data_info['service_def'][i], re.I) != None:

            uid = data_info['service_def'][i].replace('DataLink.', '')
            datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(f"https://almascience.eso.org/datalink/sync?ID={uid}")
            data_info1 = datalink.to_table().to_pandas()

            data_info = pd.concat([data_info, data_info1], ignore_index=True)

    fnamePattern = '\.pbcor\.fits$'
    if mfsOnly == True: fnamePattern = '\.mfs\..*' + fnamePattern
    if sciOnly == True: fnamePattern = '_sci\..*' + fnamePattern

    urlList = []
    fitsImageList = []

    for url in data_info['access_url']:

        fname = os.path.basename(url)

        if re.search(fnamePattern, fname) != None:

            fname = fname.replace('.pbcor.fits', '.pb.fits.gz')

            urlList1 = []

            for url1 in data_info['access_url']:

                if re.search(fname, os.path.basename(url1)) != None:
                    urlList1.append(url1)

                if len(urlList1) == 1:
                    urlList.append(url)
                    urlList.append(urlList1[0])

        if PLweblog == True and re.search('weblog', fname) != None:
            urlList.append(url)

    urlList = np.unique(urlList)

    urlList1 = []
    for url in urlList:
        urlList1.append(tuple([os.path.basename(url), url]))

    fitsImageList = []

    for i in range(len(urlList1)):

        fetch_url(urlList1[i])

        if re.search('\.pbcor\.fits$', urlList1[i][0]) != None:
            fitsImageList.append(urlList1[i][0])

        if re.search('\.gz$', urlList1[i][0]) != None:
            os.system('gunzip '+urlList1[i][0])

        if re.search('\.tgz$', urlList1[i][0]) != None:

            tar = tarfile.open(os.path.basename(urlList1[i][0]))
            weblog_path = tar.getnames()[0].split('/')[0]
            tar.extractall()
            tar.close()

            os.rename(weblog_path, os.path.basename(urlList1[i][0]).replace('.tgz', '')+'.'+weblog_path)
            os.remove(os.path.basename(urlList1[i][0]))

    return(fitsImageList)

def fetch_url(entry):
    path, uri = entry
    if not os.path.exists(path):
        r = requests.get(uri, stream=True)
        if r.status_code == 200:
            with open(path, 'wb') as f:
                for chunk in r:
                    f.write(chunk)
    return path

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


def extractSources(uidListFile=''):

    f = open(uidListFile)
    uidList = f.read().splitlines()
    f.close()

    result = ThreadPool(8).imap_unordered(extractSourcesPerMOUS, uidList)

#     for result in ThreadPool(8).imap_unordered(extractSources, uidList):
#         for ij in result:
#             fitsImageList.append(ij)


def extractSourcesPerMOUS(MOUSuid, runFilterSources=True, runMakeSourceList=True):

    MOUSuid1 = MOUSuid.replace('://', '___').replace('/', '_')

    fitsImageList = downloadImagesPerMOUS(MOUSuid)

    f = open('member.'+MOUSuid1+'.casascript.py', 'w')
    f.write("import os\n")
    f.write("import sys\n")
    f.write("import json\n")
    f.write("sys.path.append('/Users/evillard/software/ADPutils')\n")
    f.write("import ADPutils as au\n")
    f.write("\n")
    f.write("image_fname = ['")
    f.write("',\n'".join(fitsImageList))
    f.write("']\n")
    f.write("\n")
    f.write("src_dict = {}\n")
    f.write("\n")
    f.write("for i in range(len(image_fname)):\n")
    f.write("    img_pb = image_fname[i].replace('.pbcor.fits', '.pb.fits')\n")
    f.write("    if os.path.exists(img_pb) == True:\n")
    f.write("        src_dict1 = au.extractSourcesFromImage(image_fname[i], img_pb)\n")
    f.write("        if len(src_dict1.keys()) != 0:\n")
    f.write("            if len(src_dict.keys()) == 0:\n")
    f.write("                ij = 0\n")
    f.write("            else:\n")
    f.write("                ij = max(src_dict.keys()) + 1\n")
    f.write("            for j in src_dict1.keys():\n")
    f.write("                src_dict[ij+j] = src_dict1[j]\n")
    f.write("\n")
    f.write("f = open('member."+MOUSuid1+".sources.json', 'w')\n")
    f.write("json.dump(src_dict, f)\n")
    f.write("f.close()\n")
    f.close()

    subprocess.run('/usr/local/bin/casa --agg --nologger -c "member.'+MOUSuid1+'.casascript.py"', shell=True)

    if runFilterSources == True:

        f = open('member.'+MOUSuid1+'.sources.json')
        src_dict = json.load(f)
        f.close()

        src_dict1 = filterSources(src_dict)

        f = open('member.'+MOUSuid1+'.sources.json', 'w')
        json.dump(src_dict1, f)
        f.close()

    if runMakeSourceList == True:

        f = open('member.'+MOUSuid1+'.sources.json')
        src_dict = json.load(f)
        f.close()

        src_dict1 = makeSourceList(src_dict)

        f = open('member.'+MOUSuid1+'.sourcelist.json', 'w')
        json.dump(src_dict1, f)
        f.close()


def extractSourcesFromImage(img_sky_pbcor, img_pb, min_snr=6, fit_width=10, nmax_src=1000, zerooff_snr=3):

    myia = iatool()
    myqa = qatool()

    prefix1 = img_sky_pbcor.replace('pbcor.fits', '')

    img_sky_pbuncor = prefix1+'pbuncor.image'

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

            with open(prefix1+'estimates.txt', 'w') as f:
                f.write(f"{peakI}, {peakX}, {peakY}, {bMaj}{bMajUnit}, {bMin}{bMinUnit}, {bPa}{bPaUnit}\n")
            f.close()

            myia.open(img_sky_pbuncor)
        
            img_shape = myia.shape()

            found = False

            fit_width1 = round(fit_width/2.)
        
            if peakX-fit_width1 > 0 and peakY-fit_width1 > 0 and peakX+fit_width1 < img_shape[0]-1 and peakY+fit_width1 < img_shape[1]-1:
        
                src_box = str(peakX-fit_width1)+','+str(peakY-fit_width1)+','+str(peakX+fit_width1)+','+str(peakY+fit_width1)
                src_dict1 = myia.fitcomponents(box=src_box, rms=img_rms, dooff=True, estimates=prefix1+'estimates.txt', logfile=prefix1+'logfile.txt')

                if src_dict1['converged'][0] == True:

                    fit_width1 = round( 3 * src_dict1['results']['component0']['shape']['majoraxis']['value'] * max(src_dict1['pixelsperarcsec']) / 2.) # convert value of 3 into an input parameter

                    if peakX-fit_width1 > 0 and peakY-fit_width1 > 0 and peakX+fit_width1 < img_shape[0]-1 and peakY+fit_width1 < img_shape[1]-1:

                        src_box = str(peakX-fit_width1)+','+str(peakY-fit_width1)+','+str(peakX+fit_width1)+','+str(peakY+fit_width1)
                        src_dict1 = myia.fitcomponents(box=src_box, rms=img_rms, dooff=True, estimates=prefix1+'estimates.txt', logfile=prefix1+'logfile.txt')
                
                        found = True

            myia.close()

            if found == False: continue

            myia.open(img_sky_pbcor)
            if abs(src_dict1['zerooff']['value']/img_rms) < zerooff_snr:
                src_dict['component'+str(i)]['fitcomp'] = myia.fitcomponents(box=src_box, rms=img_rms, dooff=False, estimates=prefix1+'estimates.txt', logfile=prefix1+'logfile.txt')
            else:
                src_dict['component'+str(i)]['fitcomp'] = myia.fitcomponents(box=src_box, rms=img_rms, dooff=True, estimates=prefix1+'estimates.txt', logfile=prefix1+'logfile.txt')
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
            src_dict1['shape']['direction']['DEC'] = myqa.angle(src_dict1['shape']['direction']['m1'], prec=9)

            src_dict2[i] = src_dict1

        if src_dict2 != {}:
            for j in src_dict2.keys():
                src_dict2[j]['pixelcoords'] = src_dict2[j]['pixelcoords'].tolist()
                src_dict2[j]['pixelsperarcsec'] = src_dict2[j]['pixelsperarcsec'].tolist()
                src_dict2[j]['flux']['error'] = src_dict2[j]['flux']['error'].tolist()
                src_dict2[j]['flux']['value'] = src_dict2[j]['flux']['value'].tolist()
                if 'zerooff' in src_dict2[j].keys(): src_dict2[j]['zerooff']['value'] = src_dict2[j]['zerooff']['value'].tolist()

    os.system('rm -Rf '+img_sky_pbuncor)

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

def makeSourceList(src_dict, addKeywords=True):

#     keywds = ['BMAJ', 'BMIN', 'BPA', 'BUNIT', 'RADESYS', 'SPECSYS', \
#         'TELESCOP', 'DATE-OBS', 'TIMESYS', 'CASAVER', 'MEMBER', \
#         'PIPEVER', 'PROPCODE', 'SPECMODE', 'WEIGHT', 'DATE', 'ORIGIN', \
#         'ROBUST', 'IMGURL', 'CTRFREQ', 'BNDWID', 'FREQUNIT']

    src_dict1 = {}
    ij = 0

    for i in src_dict.keys():

        src_dict1[ij] = {}

        src_dict1[ij]['RA'] = src_dict[i]['shape']['direction']['RA'][0]
        src_dict1[ij]['DEC'] = src_dict[i]['shape']['direction']['DEC'][0]
        src_dict1[ij]['FLUX'] = src_dict[i]['flux']['value'][0]
        src_dict1[ij]['FLUXERR'] = src_dict[i]['flux']['error'][0]

#         if addKeywords == True:
# 
#             keyw_dict = extractMetadata(src_dict[i]['image filename'])
# 
#             for kk in keywds:
#                 if kk in keyw_dict.keys():
#                     src_dict1[ij][kk] = keyw_dict[kk]
#                 else:
#                     src_dict1[ij][kk] = ''

        ij += 1

    return(src_dict1)

def extractMetadata(fitsImageName, usePLweblog=True):

    keyw_dict = {}

    with fits.open(fitsImageName) as hdul:

        hdr_info = hdul[0].header.cards

        for j in range(len(hdr_info)):
            keyw_dict[hdr_info[j][0]] = hdr_info[j][1]

    if usePLweblog == True:

        PLweblog = glob.glob('pipeline-*')

        if len(PLweblog) == 1:
        
            df_PLweblog = readPLweblog(path=PLweblog[0])

            if 'PL image filename' in df_PLweblog.columns:

                ij = np.where(df_PLweblog['PL image filename'].str.contains(fitsImageName.replace('member.', '')) == True)[0]
                if len(ij) == 1:

                    if 'cell' in df_PLweblog.columns:

                        cellsize = re.findall('[0-9E\.]+ *arcsec', df_PLweblog.at[ij[0], 'cell'], re.IGNORECASE)
                        if len(cellsize) > 0:
                            keyw_dict['CELL'] = float(cellsize[0].upper().replace('ARCSEC', ''))

                    if 'ROBUST' not in keyw_dict.keys() and 'robust' in df_PLweblog.columns:

                        robustfactor = re.findall('-?[0-2](\.[0-9]+)?', str(df_PLweblog.at[ij[0], 'robust']))
                        if len(robustfactor) > 0:
                            keyw_dict['ROBUST'] = float(robustfactor[0])

                    if 'final theoretical sensitivity' in df_PLweblog.columns:

                        rms0 = re.findall('[0-9E\.\-]+ *m?jy/beam', df_PLweblog.at[ij[0], 'final theoretical sensitivity'], re.IGNORECASE)
                        if len(rms0) > 0:
                            rms0 = rms0[0]
                            if re.search('mjy/beam', rms0, re.IGNORECASE) != None:
                                keyw_dict['RMS0'] = float(rms0.upper().replace('MJY/BEAM', ''))/1000.
                            else:
                                if re.search('jy/beam', rms0, re.IGNORECASE) != None:
                                    keyw_dict['RMS0'] = float(rms0.upper().replace('JY/BEAM', ''))

                    if 'non-pbcor image RMS' in df_PLweblog.columns:

                        rms1 = re.findall('[0-9E\.\-]+ *m?jy/beam', df_PLweblog.at[ij[0], 'non-pbcor image RMS'], re.IGNORECASE)
                        if len(rms1) > 0:
                            rms1 = rms1[0]
                            if re.search('mjy/beam', rms1, re.IGNORECASE) != None:
                                keyw_dict['RMS1'] = float(rms1.upper().replace('MJY/BEAM', ''))/1000.
                            else:
                                if re.search('jy/beam', rms1, re.IGNORECASE) != None:
                                    keyw_dict['RMS1'] = float(rms1.upper().replace('JY/BEAM', ''))

    return(keyw_dict)
