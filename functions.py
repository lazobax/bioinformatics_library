

#reads FASTA file and returns a list of sequences
def readFastaFile(fileName):

  fileObj = open(fileName, 'rU')
  sequences = []
  seqFragments = []

  for line in fileObj:
    #if line starts with >, it is a header line
    if line.startswith('>'):
      #if there is a sequence in seqFragments, join the fragments together and add to sequences and reset the fragments
      if len(seqFragments) > 0:
        sequences.append(''.join(seqFragments))
        seqFragments = []
    #if line does not start with >, it is a sequence fragment
    else:
      seqFragments.append(line.rstrip()) #rstrip removes trailing whitespace

  #add the last sequence to sequences
  if len(seqFragments) > 0:
    sequences.append(''.join(seqFragments))

  fileObj.close()

  return sequences

#reads a PDB file and calculates the average x, y, z coordinates of the atoms in the protein
def calculateCentroid(pdbFile):
  fileObj = open(pdbFile, 'rU')

  nAtoms = xsum = ysum = zsum = 0  #initialize variables to 0

  for line in fileObj:
    if line[:6] == 'ATOM  ':#ATOM is the first 6 characters of the line
      nAtoms += 1
      xsum += float(line[30:38]) #x coordinate is characters 31-38
      ysum += float(line[38:46]) #y coordinate is characters 39-46
      zsum += float(line[46:54]) #z coordinate is characters 47-54
    fileObj.close()

    if nAtoms = 0:
      xavg = yavg = zavg = 0
    else:
      xavg = xsum / nAtoms
      yavg = ysum / nAtoms
      zavg = zsum / nAtoms
      
  return xavg, yavg, zavg

def alignSequence(Q, S):

  best_score_index = 0

  for j in range(len(S)-len(Q)):
    
    for i in range(len(Q)):


