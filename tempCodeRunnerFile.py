
  matrix = initializeMatrix(Seq1, Seq2, PENALTY)
  tracker = [[None for i in range(len(Seq1))] for j in range(len(Seq2))]

  optimalCoordinates = [0,0]
  optimal = 0

  for j in range(1, len(matrix)):
    for i in range(1, len(matrix[0])):
      diagonal = matrix[j-1][i-1] + PAM250(Seq1[i], Seq2[j])
      right = matrix[j][i-1] + PENALTY
      down = matrix[j-1][i] + PENALTY

      choice = max(diagonal, right, down)

      if choice == diagonal:
        tracker[j][i] = 'Diagonal'
      elif choice == right:
        tracker[j][i] = 'Right'
      else:
        tracker[j][i] = 'Down'

      if(choice > optimal):
        optimal = choice
        optimalCoordinates = [j, i]

      matrix[j][i] = choice


  aligned1 = ""
  aligned2 = ""
  cD = tracker[optimalCoordinates[0]][optimalCoordinates[1]]
  cC = optimalCoordinates

  while cD != None:
    if cD == 'Diagonal':
      aligned1 += Seq1[cC[1]]
      aligned2 += Seq2[cC[0]]
      cC = [cC[0]-1, cC[1]-1]
    elif cD == 'Right':
      aligned1 += Seq1[cC[1]]
      aligned2 += "-"
      cC = [cC[0], cC[1]-1]
    else:
      aligned1 += "-"
      aligned2 += Seq2[cC[0]]
      cC = [cC[0]-1, cC[1]]
    
    cD = tracker[cC[0]][cC[1]]

  visualizeMatrix(matrix)

  print(aligned1[::-1])
  print(aligned2[::-1])

