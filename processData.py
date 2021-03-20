import os
from math import sqrt
import random
import numpy as np
from scipy.spatial import cKDTree

# peakid latiCo LongiCo eleMeter Lati Longi ProMeter saddleId saddleLati longiLati isoKM
# peakid latiCo LongiCo eleMeter Lati Longi ProMeter isoKM dom


def processData(dvtPath, txtPath, outputPath, name):
    outputFile = outputPath + name + ".txt"
    peaks = []
    saddles = []
    edges = []
    runoffedges = []

    dvtfile = open(dvtPath, 'r')

    # skip comment line and G line
    dvtfile.readline()
    dvtfile.readline()

    # lat lng
    latlist = []
    lnglist = []
    while True:
        line = dvtfile.readline()
        if not line:
            break
        line = line.strip('\n')
        content = line.split(',')
        if content[0] == 'P':
            latlist.append(int(content[2]))
            lnglist.append(int(content[3]))
            peaks.append([content[1], int(content[2]), int(content[3]),
                          str(int(0.3048*int(content[4]))), '-1', '-1', '-1', '0', '1', '0', '0', '0.1'])
        elif content[0] == 'S':
            if content[2] == 'p':
                saddleid = content[1]
                lati = int(content[3])
                longi = int(content[4])
                ele = int(0.3048*int(content[5]))
                saddles.append([saddleid, lati, longi, ele])
        elif content[0] == 'R':
            continue
        elif content[0] == 'N':
            if int(content[2]) != -1:
                ridgeid = content[1]
                parent = content[2]
                saddle = content[3]
                edges.append([ridgeid, parent, saddle])
        elif content[0] == 'E':
            runoffedgesid = int(content[1])
            peakid = content[2]
            runoffedges.append([runoffedgesid, peakid])

    dvtfile.close()
    # normalization
    length = max(max(latlist)-min(latlist), max(lnglist) - min(lnglist))
    minLat = min(latlist)
    minLng = min(lnglist)
    for i in range(len(peaks)):
        peaks[i][1] = (peaks[i][1]-minLat)/length
        peaks[i][2] = (peaks[i][2]-minLng)/length

    for i in range(len(saddles)):
        saddles[i][1] = (saddles[i][1]-minLat)/length
        saddles[i][2] = (saddles[i][2]-minLng)/length
        
    txtfile = open(txtPath, 'r')
    existPeakIDs = []
    peakIDChange = dict()
    newPeaks = []
    while True:
        line = txtfile.readline()
        if not line:
            break
        peakid, lati, longi, ele, pro, saddleid, slat, slong, iso = line.split()
        newPeaks.append([len(newPeaks)+1] + peaks[int(peakid)-1][1:4] + [lati, longi, pro, saddleid, slat, slong, iso])
        peakIDChange[peakid] = len(newPeaks)
        existPeakIDs.append(peakid)
    newEdge = []
    for e in edges:
        if e[1] not in existPeakIDs or e[0] not in existPeakIDs: continue
        else:
            e[0] = peakIDChange[e[0]]
            e[1] = peakIDChange[e[1]]
            newEdge.append(e)
    txtfile.close()

    writetofile(outputFile, newPeaks, saddles, newEdge, runoffedges)


def writetofile(outputFile, peaks, saddles, edges, runoffedges):
    output = open(outputFile, 'w')
    output.write("Peaks %d\n" % (len(peaks)))
    for peak in peaks:
        output.write("%s %.4f %.4f %s %s %s %s %s %.4f\n" % (peak[0], float(peak[1]), float(
            peak[2]), peak[3], peak[4], peak[5], int(0.3048*int(peak[6])), peak[10], 0.3048*int(peak[6])/int(peak[3])))

    output.write("PromSaddles %d\n" % (len(saddles)))
    for saddle in saddles:
        output.write("%s %f %f %d\n" %
                     (saddle[0], saddle[1], saddle[2], saddle[3]))

    output.write("BasinSaddles 0\n")
    output.write("Runoffs 0\n")

    output.write("Edges %d\n" % (len(edges)))
    for edge in edges:
        output.write("%d %s %s\n" % (edge[0], edge[1], edge[2]))

    output.write("RunoffEdges %d\n" % (len(runoffedges)))
    for runoffedge in runoffedges:
        output.write("%d %s\n" % (runoffedge[0], runoffedge[1]))

    output.close()


def mergeIsolation(isoPath, promPath, outFile):

    fiso = open(isoPath)
    fpro = open(promPath)
    fout = open(outFile, 'w')

    print('Reading isolations')
    isos = []
    ni = 0
    for line in fiso:
        isos.append([float(x) for x in line.split(',')])
        ni += 1

    isoPoints = np.zeros((ni, 2))
    for i in range(ni):
        isoPoints[i, 0] = isos[i][0]
        isoPoints[i, 1] = isos[i][1]
    kdIsos = cKDTree(isoPoints)

    print('Merging isolation and prominence lists')
    numPeaks = 0
    numIsolated = 0
    for line in fpro:
        vals = [float(x) for x in line.split()]
        minDist = (3/3600)*3.3  # about 300m
        match = None

        nndist, nnid = kdIsos.query(np.array([vals[1], vals[2]]), k=1)
        if nndist < 0.25*minDist or nndist < minDist and abs(vals[3] - 0.3048*isos[nnid][3]) < 800:
            match = isos[nnid]

        if match:
            fout.write('%s %.4f\n' % (line.strip(), match[5]))
        else:
            fout.write('%s %.4f\n' % (line.strip(),  0.1))
            numIsolated += 1

        numPeaks += 1

    fiso.close()
    fpro.close()
    fout.close()

    print('Processed %d peaks, low isolation on %d' % (numPeaks, numIsolated))


if __name__ == "__main__":
    name = "J42"
    dvtPath = name + ".dvt"
    isoPath = name + "_iso.txt"
    proPath = name + "_pro.txt"
    txtPath = name + "_merge.txt"
    mergeIsolation(isoPath, proPath, txtPath)
    outputPath = ""
    processData(dvtPath, txtPath, outputPath, name)
