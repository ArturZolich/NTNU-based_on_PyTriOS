from ramses_calibrate import *
from python import *

import sys
import time

import struct

import numpy as np


from pytrios.TClasses import TProtocolError, TPackMeasKeyError, TPacket, TSerial, TCommandSend, TChannel


def main():
    with open("C:\\Users\\arturz\\Downloads\\23_06_14_11_25_49_08.txt", mode='rb') as file: # b is important -> binary       
      
        fileContent = file.read()
        s, s2 =_get_s2parse(fileContent)
        s, s2 =_get_s2parse(s)
      #  while(s2 != None):
        while(s != None):
            packet = TPacket(s2)
            handlePacket(packet)
            s, s2 =_get_s2parse(s)
           # print(packet.packetType)

if __name__ == '__main__':
    main()