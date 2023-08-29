import sys
import time

import struct

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.colors import LinearSegmentedColormap

from pytrios.TClasses import TProtocolError, TPackMeasKeyError, TPacket, TSerial, TCommandSend, TChannel
from ramses_calibrate import *

g_currentRamses = 0
g_Ram1Success = 0


class Ini(object):
    def __init__(self, DeviceType=None, SensorName=None,
                 SAMDevice=None, DeviceTypeSub1=None, DeviceTypeSub2=None,
                 DarkPixelStart=None, DarkPixelStop=None, Reverse=None,
                 WavelengthRange=None, c0s=None,
                 c1s=None, c2s=None, c3s=None, cs=None):
        pass

class Cal(object):
    def __init__(self, SAMDevice_Aqua=None, SAMDevice_Air=None,
                 SAMDevice_Back=None, SAMDateTime_Aqua=None,
                 SAMDateTime_Air=None, SAMDateTime_Back=None,
                 SAMspectrum_Aqua=None, SAMspectrum_Air=None,
                 SAMspectrum_Back0=None, SAMspectrum_Back1=None, ini=Ini):
        self.ini = ini()



tchannels = TChannel()
g_calOut_513 = Cal()
g_out1 = None
g_spec1 = None
g_out1_RAW = None

g_Roll = None
g_Pitch = None
g_Yaw = None

g_Time = None
g_Date = None

g_Lat = None
g_Lon = None

def wavelength_to_rgb(wavelength, gamma=0.8):
    ''' taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range
    '''
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 750:
        A = 1.
    else:
        A = 0.5
    if wavelength < 380:
        wavelength = 380.
    if wavelength > 750:
        wavelength = 750.
    if 380 <= wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif 440 <= wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif 490 <= wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif 510 <= wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif 580 <= wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif 645 <= wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R, G, B, A)


def SAMInterpreter(regch, packet):
    global g_currentRamses
    global g_out1
    global g_spec1

    formatstring = '<'+'H'*int(packet.id1_databytes/2)
    rawdata = bytearray(y for y in packet.databytes)
    LEdata = struct.unpack(formatstring, rawdata)
    """ sloppy code comment:
    In the following, if we place LEdata directly into the dataframes slice
    it will be overwritten upon prompt arrival of a new packet, even if this
    concerns a different entry of the tchannels dictionary. Odd!
    Reading and writing back the list of dataframes first circumvents this.
    Packets arriving in rapid succession re-use memory blocks. So this
    possibly points to sloppy garbage collection?
    """
    dataframes = regch.TSAM.dataframes[:]
    dataframes[packet.framebyte] = LEdata
    regch.TSAM.dataframes = dataframes
   # print("SAMInterpreter: Spectrum framebyte {0}"
   #           .format(packet.framebyte))
    if packet.framebyte == 0:
        frames = regch.TSAM.dataframes
        if sum(y is None for y in frames) == 0:
            outspec = []
            for sublist in frames:
                sl = list(sublist)
                sl.reverse()
                outspec = outspec+sl
            outspec.reverse()  # assuming this is not a UV sensor..
            regch.TSAM.lastRawSAM = outspec
            regch.TSAM.lastRawSAMTime = packet.timeStampPC
            msintt = 2*2**(outspec[0] & 0b1111)  # integration time
            regch.TSAM.lastIntTime = msintt
            # reset to receive the next spectrum
            regch.TSAM.dataframes = [[None]]*8
            
            print("INT_TIME: ", msintt, "\t DATA")
           # print(*regch.TSAM.lastRawSAM, sep='\t')
            x = np.arange(350, 951, 1)
            #plt.plot(x, regch.TSAM.lastRawSAM)
            #x = np.arange(346, 955, 3.3)

            global g_Ram1Success
            global g_out1_RAW
            
            if(g_currentRamses == 1):
                try:
                    g_out1_RAW = regch.TSAM.lastRawSAM
                    g_out1 = raw2cal_Air(regch.TSAM.lastRawSAM, msintt, "513D", g_calOut_513, x)
                    g_Ram1Success=1
                except Exception as e: 
                    print(e)
                    print("\t\t\t\t\t\t\t\tRAM 1 FAILURE")
            if(g_currentRamses == 2) and g_Ram1Success==1:
                g_Ram1Success=0
                out = None
                spec = None
                try:
                    out = raw2cal_Air(regch.TSAM.lastRawSAM, msintt, "8027", g_calOut_513, x)

                    print(x.shape)
                    print(len(out))
                    print(len(g_out1))

                    fig, axs = plt.subplots(3,1, height_ratios=[2,1,1], figsize=(7,7)) 
                    clim = (380, 750)
                    norm = plt.Normalize(*clim)
                    wl = np.arange(clim[0], clim[1] + 1, 2)
                    colorlist = list(zip(norm(wl), [wavelength_to_rgb(w) for w in wl]))
                    spectralmap = LinearSegmentedColormap.from_list("spectrum", colorlist)

                    wavelengths = x
                    spectrum = out
                    axs[1].plot(wavelengths, spectrum, color='black', linewidth=1)

                    #wavelengths = g_spec1 #x

                    maxRad = 1.3 * max(g_out1)
                    maxIrRad = 1.1 * max(out) 

                    y = np.linspace(0, maxIrRad, 100)
                    X, Y = np.meshgrid(wavelengths, y)

                    y_low = np.linspace(0, maxRad, 100)
                    X_low, Y_low = np.meshgrid(wavelengths, y)

                    extent = (np.min(wavelengths), np.max(wavelengths), np.min(y), np.max(y))
                    extent_low = (np.min(wavelengths), np.max(wavelengths), np.min(y_low), np.max(y_low))

                    axs[0].imshow(X, clim=clim, extent=extent_low, cmap=spectralmap, aspect='auto')
                    axs[0].plot(wavelengths, g_out1, color='black', linewidth=1)
                    axs[0].set_title('Radiance, sensor ID: 513D')
                    axs[0].set_xlabel('Wavelength [nm]')
                    axs[0].set_ylabel('mW m^-2 nm Sr^1')     
                    #axs[0].set_xlim(400,700)           
                    axs[0].set_ylim(0)
                    
                    axs[0].fill_between(wavelengths, g_out1, maxRad, color='w')


                    axs[1].imshow(X, clim=clim, extent=extent, cmap=spectralmap, aspect='auto')                
                    axs[1].set_title('Irradiance, sensor ID: 8027')
                    axs[1].set_xlabel('Wavelength [nm]')
                    axs[1].set_ylabel('mW m^-2 nm^1')  
                    #axs[1].set_xlim(400,700)         
                    axs[1].set_ylim(0)

                    axs[1].fill_between(wavelengths, spectrum, maxIrRad, color='w')

                    global g_Time, g_Date, g_Lat, g_Lon, g_Roll, g_Pitch, g_Yaw

                    text = 'Time (local): ' + str(g_Time)+"\n"
                    text += 'Date: ' + str(g_Date) +"\n"+"\n"
                    text += 'Lat: ' + str(g_Lat) +"\n"
                    text += 'Long: ' + str(g_Lon) +"\n"+"\n"
                    text += 'Hdg: ' + str(g_Yaw) 
                    text += 'Roll: ' + str(g_Roll)  +"\n"
                    text += 'Pitch: ' + str(g_Pitch)

                    axs[0].text(770, 0.2*maxRad, text, style='italic', 
                            bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
                    
                    plt.annotate("artur.zolich@ntnu.no", (0,0), (300,-35), xycoords='axes fraction', textcoords='offset points', va='top')

                    fig.tight_layout()
                    #plt.show()


                    filename = str(g_Time) + '_' + str(g_Date)
                    filename = filename.replace(' ','_') 
                    filename = filename.replace(':','_') 
                    filename2 = filename.replace('.','_') + ".txt"
                    filename3 = "RAW_" + filename.replace('.','_') + ".txt"
                    filename = filename.replace('.','_') + ".png"

                    radiance = []
                    for i in range(601):
                        radiance.append(g_out1[i] / out[i])
                    
                    #print("RADIANCE: ")
                    #print(radiance)
                    
                   # txt = axs[2].text(630, max(radiance)*3, "DRAFT",
                     #           ha="center", va="center", rotation=45, size=125,
                      #          bbox=dict(boxstyle="round4,pad=0.5",
                       #         fc="lightcoral", ec="red", lw=4, alpha=0.5))
                    #txt.set_alpha(0.3)

                    axs[2].plot(x, radiance, color='black', linewidth=1)
                    axs[2].set_title('Water leaving reflectance')
                    axs[2].set_xlabel('Wavelength [nm]')
                    axs[2].set_ylabel('Sr^-1')                    
                    #axs[2].set_xlim(400,700)
                    plt.savefig(os.path.join("C:\\Users\\arturz\\Downloads\\RAMSES_data\\Graphs\\" , filename))
                    
                    with open(os.path.join("C:\\Users\\arturz\\Downloads\\RAMSES_data\\DATA\\" , filename2), 'w') as f:   
                        f.write("Pamela + RAMSES data for hyperspectral satellite calibration")
                        f.write("\n")
                        f.write("Artur Zolich, artur.zolich@ntnu.no, v0.2, June 2023")
                        f.write("\n")
                        f.write(text)  
                        f.write("\n")
                        f.write("Irradiance sensor ID: 8027")
                        f.write("\n")
                        f.write("Radiance sensor ID: 513D")
                        f.write("\n")
                        f.write("Wavelength\tIrradiance\tRadiance\tReflectance")
                        f.write("\n")
                        f.write("[nm]\t[mW m^-2 nm Sr^1]\t[mW m^-2 nm^1]\t[Sr^-1]")
                        f.write("\n")
                        for i in range(601):
                            line = str(x[i]) + '\t' + str(out[i]) + '\t' + str(g_out1[i]) + '\t' + str(radiance[i])
                            f.write(line)
                            f.write("\n")
                
				
                    with open(os.path.join("C:\\Users\\arturz\\Downloads\\RAMSES_data\\DATA\\" , filename3), 'w') as f:   
                        f.write("Pamela + RAMSES data for hyperspectral satellite calibration")
                        f.write("\n")
                        f.write("Artur Zolich, artur.zolich@ntnu.no, v0.2, June 2023")
                        f.write("\n")
                        f.write(text)
                        f.write("\n")
                        f.write("Irradiance sensor ID: 8027")
                        f.write("\n")
                        f.write("Radiance sensor ID: 513D")
                        f.write("\n")
                        f.write("Wavelength\tIrradiance\tRadiance\tReflectance")
                        f.write("\n")
                        f.write("[nm]\t[counts]\t[counts]")
                        f.write("\n")
                        for i in range(255):
                            line = str(i) + '\t' + str(regch.TSAM.lastRawSAM[i]) + '\t' + str(g_out1_RAW[i])
                            #print(line)
                            f.write(line)
                            f.write("\n")

                  

                


                except Exception as e: 
                    print(e)
                    print("\t\t\t\t\t\t\t\tRAM 2 FAILURE")

                
          

           # print("SAMInterpreter: Spectrum:{0}, ({1} s)"
            #          .format(regch.TSAM.lastRawSAM, msintt))
        else:
            emsg = "SAM Interpreter: Incomplete spectrum, discarded"
            print(emsg)
            raise TProtocolError(emsg)
            # reset to receive the next spectrum
            regch.TSAM.dataframes = [[None]]*8
    return regch


def handlePacket(packet):
    global tchannels
    "Directs incoming packets to appropriate interpreters, updating tchannels"
    p = packet  # shorten

    if p.packetType == 'measurement':
        try:
            ch = SAMInterpreter(tchannels, p)
            tchannels = ch
        except:
            print("BAD PACKET")





def TStrRepl(s):
    s = s.replace(b'@g', b'\x13')  # correct for escape chars (xOFF)
    s = s.replace(b'@f', b'\x11')  # correct for escape chars (xOn)
    s = s.replace(b'@e', b'\x23')  # correct for escape cha,rs (data start #)
    s = s.replace(b'@d', b'\x40')  # correct for escape chars (escape char @)
    return s



def _get_s2parse(s):
    global g_currentRamses
    "extract data blocks from serial buffer"

    first, last = s.find(b'#'), s.rfind(b'#')
    s = TStrRepl(s)  # correct replacement chars
    if first < 0 or not last >= first:
        return None, None

    s = s[s.find(b'#', 0):]  # omit incomplete sequence at start
    if len(s) <= 1:  # 1st byte after # = size
        return None, None
    
    if(s[1] != 43):
        #print("Zolich: {0}", s[1])
        ndatabytes = 2*2**(s[1] >> 5)
        blocklength = 8+ndatabytes
        if len(s) >= blocklength:
            s2parse = s[1:blocklength]  # block to parse
        
            prettyhex = ":".join("{0:x}".format(c) for c in s2parse)
          #  print("TListen: {0}".format(prettyhex))
            s = s[blocklength:]  # remainder to next cycle
            return s, s2parse
        else:
            return s, None
    else:
        s3 = str(s[:s.find(b'#', 1)])
        if "CGNSINF" in s3:
            tmp = s3.split(",")
            tmp[20]=tmp[20].replace('\"',',')
            tmp[20]=tmp[20].replace(' ',',')
            tmp[20]=tmp[20].replace('+',',')
            time = tmp[20].split(',')
            datmp = tmp[19].split('\"')
            data = datmp[1].split('/')
            data_txt = str(data[2] + "/" + data[1] + "/20" + data[0])
            print(tmp[3], " ", tmp[4], " ", time[0], " ", data_txt)

            global g_Time, g_Date, g_Lat, g_Lon, g_Roll, g_Pitch, g_Yaw
            g_Time = time[0] #tmp3[0]
            g_Date = data_txt.replace('/','.')
            g_Lon = np.nan
            g_Lat = np.nan
            g_Roll = np.nan
            g_Pitch = np.nan
            g_Yaw = np.nan

            with open("C:\\Users\\arturz\\Downloads\\PAMELA_GPS_RPY_UTC.txt") as file:
                while((data := file.readline())):
                    if(time[0] in data) and (data_txt in data):
                        tmp3 = data.split("\t")
                        g_Lon = tmp3[4]
                        g_Lat = tmp3[3]
                        g_Roll = tmp3[5]
                        g_Pitch = tmp3[6]
                        g_Yaw = tmp3[7]
                        print(data)

        if "RAM_1" in s3:
            print("RAM_1")
            global g_Ram1Success 
            g_Ram1Success = 0
            g_currentRamses = 1

        if "RAM_2" in s3: 
            print("RAM_2")
            g_currentRamses = 2

           # print(*tmp, end="\t")

        s = s[s.find(b'#', 1):]  # omit incomplete sequence at start
        return s, None
  
def main():
     
    #foldername = "C:\\Users\\arturz\\Downloads\\RAMSES_data\\"
    #foldername = "C:\\Users\\arturz\\Downloads\\RAMSES_data\\DIRECT\\DAT\\"
    #foldername = "C:\\Users\\arturz\\Downloads\\RAMSES_data\\VERIFICATION\\INPUT"
    foldername = "C:\\Users\\arturz\\Downloads\\PAM1807\\"
    
    files = os.listdir(foldername)
    for f in files:

        if f.endswith('.txt'):


        #    with open("C:\\Users\\arturz\\Downloads\\RAMSES_8027\\ramses2_startup_query_single_run.txt", mode='rb') as file: # b is important -> binary       
        # with open("C:\\23_05_23_21_27_28_08.txt", mode='rb') as file: # b is important -> binary
        #  with open("C:\\23_06_03_19_44_45_08.txt", mode='rb') as file: # b is important -> binary
            #with open("C:\\Users\\arturz\\Downloads\\MOD_23_06_15_12_21_49_08.txt", mode='rb') as file: # b is important -> binary       
            #with open("C:\\Users\\arturz\\Downloads\\23_06_15_11_18_13_08.txt", mode='rb') as file: # b is important -> binary       
            
            with open(os.path.join(foldername, f), mode='rb') as file: 
                fileContent = file.read()
                s, s2 =_get_s2parse(fileContent)
                s, s2 =_get_s2parse(s)
            #  while(s2 != None):
                while(s != None):
                    packet = TPacket(s2)
                    handlePacket(packet)
                    s, s2 =_get_s2parse(s)
                # print(packet.packetType)

                
                #fig.tight_layout()
                #plt.show()

if __name__ == '__main__':
    main()