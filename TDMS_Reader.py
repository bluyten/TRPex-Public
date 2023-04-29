import csv
import numpy as np
from nptdms import TdmsFile
from matplotlib import pyplot as plt

g0 = 9.80665            # m/s^2

# FORMAT = ([Volts/Amps], [Newton/Pascal])
pressureCalibration = ([], [1.01325*10**5, 8.01325*10**5, 1.01325*10**5])
thrustCalibration = ([], [0, 3.938*g0, 0, 17.593*g0, 0])
calibration = {"PS0": pressureCalibration, "LC1": thrustCalibration}

def readChannel(file, channel, start, stop, dt, n, calibration):

    group = channel[:2] + " data"

    # Read Raw Data From TDMS File
    with TdmsFile.open("Data\\" + file) as tdmsFile:
        chan = tdmsFile[group][channel]
        dataRaw = np.array(chan[start:stop])

    tRaw = [i*dt for i in range(len(dataRaw))]

    # Filter Data
    tFiltered, dataFiltered = filterData(tRaw, dataRaw, n)
    
    # Calibration Data
    try:
        # Force to perform new calibration
        if calibration:
            raise Exception
        readCalibrationData(channel)
        print("Calibration Data Read")
    except: 
        print("Performing Calibration")
        calibStop = {"PS0": 13*10**5, "LC1": 20*10**5}
        calibM = {"PS0": 10**4, "LC1": 10**4}
        calibK = {"PS0": 10**4, "LC1": 2*10**4}
        calibThreshold = {"PS0": 9*10**(-5), "LC1": 1*10**(-5)}
        calibrate(channel, 0, calibStop[channel], dt, n, calibM[channel], calibK[channel], calibThreshold[channel], plot = False)
        writeCalibrationData(channel)
        print("Calibration Complete")
    
    # Calibrate Data
    tScaled, dataScaled = tFiltered, scaleData(channel, dataFiltered)

    return tScaled, dataScaled

    
# Filter Data Using N Nearest Neighbours
def filterData(tRaw, dataRaw, n):
    if n%2 != 0:
        n += 1

    dataFiltered = np.array([np.sum(dataRaw[i - n//2 : i + n//2])/n for i in range(n//2, len(dataRaw) - n//2)])   
    tFiltered = tRaw[n//2 : -n//2]

    return tFiltered, dataFiltered

def scaleData(channel, data):
    calib = calibration[channel]

    dataScaled = []

    for i in range(len(data)):
        # Find Corresponding Two Interpolation Points (Also Works for Extrapolation)
        j = 1
        while data[i] > calib[0][j]:
            if j + 1 == len(calib[0]):
                break
            j += 1
        
        dataScaled += [interpolate(data[i], calib[0][j-1], calib[0][j], calib[1][j-1], calib[1][j])]
    
    return np.array(dataScaled)

def interpolate(point, x1, x2, y1, y2):
    # x1 < point < x2
    return y1 + (point - x1)/(x2 - x1) * (y2 - y1)

def calibrate(channel, start, stop, dt, n, m, k, threshold, plot = True):

    group = channel[:2] + " data"

    # Read Raw Data From TDMS File
    with TdmsFile.open("Data\Calibration.tdms") as tdmsFile:
        chan = tdmsFile[group][channel[:2] + "0"]
        dataRaw = np.array(chan[start:stop])

    tRaw = [i*dt for i in range(len(dataRaw))]

    # Divide Data in Step-wise Pieces
    pieces = []
    pieceData = [dataRaw[0]]
    pieceTime = [tRaw[0]]

    for i in range(len(dataRaw) - m):
        diff = np.abs(dataRaw[i + m] - dataRaw[i])
        if diff < threshold:
            # Still In Same Piece
            pieceData += [dataRaw[i + m]]
            pieceTime += [tRaw[i + m]]
        else:
            if len(pieceData[k:-k]) > m:
                # New Piece Created
                pieces += [(pieceTime[k:-k], pieceData[k:-k])]
            # Piece too short --> Discarded
            pieceData = []
            pieceTime = []
    pieces += [(pieceTime[k:-k], pieceData[k:-k])]

    # Plot Pieces
    if plot:
        # Filter Data
        for i in range(len(pieces)):
            tFiltered, dataFiltered = filterData(pieces[i][0], pieces[i][1], n)
            pieces[i] = (tFiltered, dataFiltered)

        conversion = 10**3 if channel == "PS0" else 10**4

        # Plot Data
        plt.figure()
        ax = plt.axes()
        ax.plot(np.array(tRaw)*4, np.array(dataRaw)*conversion, label = "Raw")
        for i in range(len(pieces) - 1):
            t = pieces[i][0]
            data = pieces[i][1]
            ax.plot(np.array(t)*4, np.array(data)*conversion, 'C1')
        t = pieces[-1][0]
        data = pieces[-1][1]
        ax.plot(np.array(t)*4, np.array(data)*conversion, 'C1', label = "Calibration Steps")
        ax.set_xlabel(r'Ticks [$10^5$]')
        if channel == "PS0":
            ax.set_ylabel(r'Current [$mA$]')
        else:
            ax.set_ylabel(r'Voltage Ratio [$10^{-4}$]')
        ax.legend()
        plt.show()
    
    global calibration, thrustCalibration, pressureCalibration
    # Calculate Average of Pieces
    averages = [np.average(piece[1]) for piece in pieces]
    data = calibration[channel][1]

    # Sort Calibration Points
    averagesSorted = [0]
    dataSorted = [0]
    for j, average in enumerate(averages):
        i = 0
        while average > averagesSorted[i] and i < len(averagesSorted) - 1:
            i += 1
        averagesSorted.insert(i, average)
        dataSorted.insert(i, data[j])
    
    calibration[channel] = (averagesSorted[:-1], dataSorted[:-1])

# Write Calibration Data Externally
def writeCalibrationData(channel):
    name = 'Data\CalibrationData' + channel + '.csv'
    with open(name, 'w', newline = '') as file:
        writer = csv.writer(file)
        writer.writerow([channel[:-1] + " Channel"])
        a = calibration[channel][0]
        b = calibration[channel][1]
        writer.writerow(a)
        writer.writerow(b)
# Read Calibration Data Externally
def readCalibrationData(channel):
    name = 'Data\CalibrationData' + channel + '.csv'
    with open(name, newline = '') as file:
        reader = csv.reader(file, delimiter = ',', quotechar = '|')
        for i, row in enumerate(reader):
            if i == 1:
                a = [float(val) for val in row]
            elif i == 2:
                b = [float(val) for val in row]
        calibration[channel] = (a, b)        
