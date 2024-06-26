import struct, sys
from optparse import OptionParser

parser = OptionParser("data_checker.py input output")
parser.add_option("--orbitsPerLumi", dest="orbitsPerLumi", type=int, default=(1<<10), help="orbits per lumisection")
parser.add_option("--run", dest="run", type=int, default=37, help="run number")
options, args = parser.parse_args()

def writeFRDFileHeader_v2(fout, run, lumi, events, size, data_type=20):
  fout.write("RAW_0002".encode())
  fout.write(struct.pack('H',32)) # header_size
  fout.write(struct.pack('H',data_type))
  fout.write(struct.pack('I',events))
  fout.write(struct.pack('I',run))
  fout.write(struct.pack('I',lumi))
  fout.write(struct.pack('Q',size)) 
fileHeaderSize = 32

def writeFRDEventHeader_V6(fout, run, lumi, event, size, flags=0, crc32c=0, version=6):
  fout.write(struct.pack('H',version))
  fout.write(struct.pack('H',flags))
  fout.write(struct.pack('I',run))
  fout.write(struct.pack('I',lumi))
  fout.write(struct.pack('I',event)) 
  fout.write(struct.pack('I',size))
  fout.write(struct.pack('I',crc32c))
orbitHeaderSize = 2*2+5*4  

def readEventHeader(infile, checkBits=False):
    data = infile.read(8)
    if len(data) < 8: raise EOFError()
    evhead = int.from_bytes(data, sys.byteorder)
    if evhead == 0:
        return evhead, 0, 0, 0
    if checkBits:
      bits = evhead >> 62
      if bits != 0b10:
          raise RuntimeError("Bad event header %016x, bits 63-62 are %1d%1d, expected 10" % (evhead,bits/2,bits%2))
    orbit = (evhead >> 24) & 0xFFFFFFFF
    bx = (evhead >> 12) & 0xFFF
    nwords = evhead & 0xFFF;
    return evhead, orbit, bx, nwords 

infile = open(args[0], 'rb')
outfile = open(args[1], 'wb')
orbits, events, bytes = 0, 0, fileHeaderSize
run = options.run
firstlumi = None

## write an empty header
writeFRDFileHeader_v2(outfile, 0, 0, 0, 0)
oldorbit = -1
orbitbuffer = []
try:
  while infile.readable():
    evhead, orbit, bx, nwords = 0, 0, 0, 0
    while infile.readable() and evhead == 0:
      evhead, orbit, bx, nwords = readEventHeader(infile)
    if orbit != oldorbit:
      if oldorbit != -1:
        orbitSize = sum(len(x) for x in orbitbuffer)
        lumi = oldorbit//options.orbitsPerLumi
        writeFRDEventHeader_V6(outfile, run, lumi, oldorbit, orbitSize)
        for data in orbitbuffer:
          outfile.write(data)
        orbits += 1
        bytes += orbitSize + orbitHeaderSize
        #print("Wrote orbit %d (%u) of size %u (payload) + %u (header) = %u" % (orbits-1, oldorbit, orbitSize, orbitHeaderSize, orbitSize +  orbitHeaderSize))
      else:
        firstlumi = 1 + orbit//options.orbitsPerLumi
      orbitbuffer = []
      oldorbit = orbit
    events += 1
    orbitbuffer.append(evhead.to_bytes(8, sys.byteorder))
    if nwords > 0:
      orbitbuffer.append(infile.read(nwords * 8))
except EOFError:
  if orbitbuffer != []:
    orbitSize = sum(len(x) for x in orbitbuffer)
    lumi = 1 + oldorbit//options.orbitsPerLumi
    writeFRDEventHeader_V6(outfile, run, lumi, oldorbit, orbitSize)
    for data in orbitbuffer:
      outfile.write(data)
    orbits += 1
    bytes += orbitSize + orbitHeaderSize
    #print("Wrote orbit %d (%u) of size %u (payload) + %u (header) = %u" % (orbits-1, oldorbit, orbitSize, orbitHeaderSize, orbitSize +  orbitHeaderSize))
outfile.seek(0)
writeFRDFileHeader_v2(outfile, run, firstlumi, orbits, bytes)
print("Wrote %u orbits, %u events, %u bytes" % (orbits, events, bytes))