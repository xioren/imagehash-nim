import std/math


const threshold: float = 1e-13

##########################################################################

proc roundToSignificantDigits(value: float, digits: int): float =
  if value == 0.0:
    return 0.0
  let d = ceil(log10(abs(value)))
  let power = digits - int(d)
  let magnitude = pow(10.0, power.float)
  round(value * magnitude) / magnitude


proc transpose(matrix: seq[seq[float]]): seq[seq[float]] =
  let numRows = matrix.len
  let numCols = matrix[0].len
  var transposed = newSeq[seq[float]](numCols)

  for i in 0 ..< numCols:
    transposed[i] = newSeq[float](numRows)
    for j in 0 ..< numRows:
      transposed[i][j] = matrix[j][i]

  return transposed


proc roundDown(n: float): float =
  ## values very close to zero are returned as zero
  if abs(n) < threshold:
    return 0.0
  return n
  

##########################################################################


proc dct1d*[T](vector: seq[T], ortho: bool = false): seq[float] =
  ## DCT II
  let N = vector.len.float
  result = newSeq[float](vector.len)

  for k in 0 ..< vector.len:
    var sum = 0.0
    for n in 0 ..< vector.len:
      sum += vector[n].float * cos(PI * k.float * (2.0 * n.float + 1.0) / (2.0 * N))
    
    if ortho:
      if k == 0:
        result[k] = roundDown((sqrt(1.0 / (4 * N)) * 2.0) * sum)
      else:
        result[k] = roundDown((sqrt(1.0 / (2 * N)) * 2.0) * sum)
    else:
      result[k] = roundDown(2.0 * sum)

  return result


proc idct1d*[T](vector: seq[T], ortho: bool = false): seq[float] =
  ## DCT III
  let N = vector.len.float
  result = newSeq[float](vector.len)

  for k in 0 ..< vector.len:
    var sum = 0.0
    for n in 1 ..< vector.len:
      sum += vector[n] * cos(PI * n.float * (2.0 * k.float + 1.0) / (2.0 * N))
    
    if ortho:
      result[k] = roundDown((vector[0]/sqrt(N)) + (sqrt(2.0/N)) * sum)
    else:
      result[k] = roundDown(vector[0] + (2.0 * sum))

  return result


proc dct2d*[T](matrix: seq[seq[T]], ortho: bool = false): seq[seq[float]] =
  ## DCT II
  let N = matrix.len
  var intermediate = newSeq[seq[float]](N)
  result = newSeq[seq[float]](N)

  # NOTE: apply 1D DCT to each row
  for i in 0 ..< N:
    intermediate[i] = dct1d(matrix[i], ortho)

  # NOTE: transpose the matrix
  let transposed = transpose(intermediate)

  # NOTE: apply 1D DCT to each transposed row (original column)
  for i in 0 ..< N:
    result[i] = dct1d(transposed[i], ortho)

  return transpose(result)


proc idct2d*[T](matrix: seq[seq[T]], ortho: bool = false): seq[seq[float]] =
  ## DCT III
  let N = matrix.len
  var intermediate = newSeq[seq[float]](N)
  result = newSeq[seq[float]](N)

  # NOTE: apply 1D IDCT to each row
  for i in 0 ..< N:
    intermediate[i] = idct1d(matrix[i], ortho)

  # NOTE: transpose the matrix
  let transposed = transpose(intermediate)

  # NOTE: apply 1D IDCT to each transposed row (original column)
  for i in 0 ..< N:
    result[i] = idct1d(transposed[i], ortho)

  # NOTE: transpose back to get the final result
  return transpose(result)

##########################################################################
