import std/[algorithm, math, strutils]

import imageman/imageman
import dct

export imageman


const hashSize: int = 8
const hashBits: int = 64

type
  ImageHash* = array[hashBits, int]
  Matrix[T] = seq[seq[T]]

#################################################################################

proc toUInt64(bits: openArray[int]): uint64 =
  if bits.len != 64:
    raise newException(ValueError, "sequence must be 64 bits long")
  
  for i in 0 ..< bits.len:
    if bits[i] == 1:
      result = result or (uint64(1) shl (63 - i))
  return result


proc popCount(x: uint64): int =
  ## counts the number of set bits (1s) in x
  var y = x
  while y != 0:
    inc result
    # NOTE: clear the least significant bit set
    y = y and (y - 1)


proc hammingDistance*(x, y: ImageHash): int =
  ## compute the Hamming distance between two uint64 hashes
  return popCount(toUInt64(x) xor toUInt64(y))


proc hammingDistanceN*(x, y: ImageHash): float =
  ## compute the normalized Hamming distance between two uint64 hashes
  return 1 / popCount(toUInt64(x) xor toUInt64(y))


proc slice[T](arr: Matrix[T], startRow, endRow, startCol, endCol: int): Matrix[T] =
  if not (startRow <= endRow and startCol <= endCol):
    raise newException(ValueError, "invalid slice indices")
  
  result = newSeq[seq[T]](endRow - startRow)
  for i in startRow ..< endRow:
    result[i - startRow] = newSeq[T](endCol - startCol)
    for j in startCol ..< endCol:
      result[i - startRow][j - startCol] = arr[i][j]


proc asArray(image: Image): Matrix[uint8] =
  ## convert greyscale image into 2d array of uint8
  if not (image.data.len == image.width * image.height):
    raise newException(ValueError, "pixel count does not match the image dimensions.")
  
  result = newSeq[seq[uint8]](image.height)
  for i in 0 ..< image.height:
    result[i] = newSeq[uint8](image.width)
  
  for y in 0 ..< image.height:
    for x in 0 ..< image.width:
      let idx = y * image.width + x
      # QUESTION: should we convert to float?
      result[y][x] = image.data[idx].r


proc haarWaveletTransform[T](data: Matrix[T], depth: int = 1): tuple[LL, LH, HL, HH: Matrix[float]] =
  ## LL: Low-Low   (approximation)
  ## LH: Low-High  (horizontal detail)
  ## HL: High-Low  (vertical   detail)
  ## HH: High-High (diagonal   detail)
  let rows = data.len
  let cols = data[0].len
  var LL = newSeq[seq[float]](rows div 2)
  var LH = newSeq[seq[float]](rows div 2)
  var HL = newSeq[seq[float]](rows div 2)
  var HH = newSeq[seq[float]](rows div 2)

  # NOTE: initialize matrices
  for i in 0 ..< rows div 2:
    LL[i] = newSeq[float](cols div 2)
    LH[i] = newSeq[float](cols div 2)
    HL[i] = newSeq[float](cols div 2)
    HH[i] = newSeq[float](cols div 2)

  # NOTE: compute Haar wavelet transform
  for i in 0 ..< rows div 2:
    for j in 0 ..< cols div 2:
      let a = data[i * 2    ][j * 2    ]
      let b = data[i * 2    ][j * 2 + 1]
      let c = data[i * 2 + 1][j * 2    ]
      let d = data[i * 2 + 1][j * 2 + 1]

      LL[i][j] = float((a + b + c + d)) / 4.0
      LH[i][j] = float((a - b + c - d)) / 4.0
      HL[i][j] = float((a + b - c - d)) / 4.0
      HH[i][j] = float((a - b - c + d)) / 4.0

  if depth > 1:
    return haarWaveletTransform(LL, depth.pred)
  else:
    return (LL, LH, HL, HH)


proc hex(x: ImageHash): string =
  return toUInt64(x).toHex().toLowerAscii()


proc `$`*(x: ImageHash): string =
  return x.hex()


proc `==`*(x, y: ImageHash): bool =
  # QUESTION: should we convert it uint64 or just compare the bits?
  return toUInt64(x) == toUInt64(y)


proc `-`*(x, y: ImageHash): int =
  return hammingDistance(x, y)


proc `--`*(x, y: ImageHash): float =
  return hammingDistanceN(x, y)


proc flatten[T](array2D: Matrix[T]): seq[T] =
  ## convert 2d matrix into 1d sequence
  for subArray in array2D:
    for item in subArray:
      result.add(item)


proc sumIt[T](data: openArray[T]): int =
  for n in data:
    result.inc(n.int)


proc mean[T](data: openArray[T]): float =
  return data.sumIt / data.len


proc median[T](data: openArray[T]): float =
  if data.len == 0:
    raise newException(ValueError, "data must not be empty")
  let sortedData = data.sorted
  let midIndex = data.len div 2
  if data.len mod 2 == 0:
    # NOTE: even number of elements -> return the average of the two middle elements
    return (sortedData[midIndex.pred] + sortedData[midIndex]) / 2.0
  else:
    # NOTE: odd number of elements -> return the middle element
    return sortedData[midIndex].float

#################################################################################

proc aHash*(image: Image): ImageHash =
  ## Average Hash: compares images based on their overall brightness levels
  var reducedImage = image
  reducedImage.filterGreyscale()
  reducedImage = reducedImage.resizedNN(hashSize, hashSize)

  let data = reducedImage.asArray()
  let avg = mean(data.flatten())

  var diff: array[hashBits, int]
  var idx: int
  for i in 0 ..< hashSize:
    for j in 0 ..< hashSize:
      idx = i * hashSize + j
      if data[i][j].float > avg:
        diff[idx] = 1
      else:
        diff[idx] = 0

  return ImageHash(diff)


proc dHash*(image: Image): ImageHash =
  ## Difference Hash: detects changes in edge information by comparing the brightness of adjacent pixels,
  ## suitable for identifying images with similar structures or orientations.

  var reducedImage = image
  reducedImage.filterGreyscale()
  reducedImage = reducedImage.resizedNN(hashSize.succ, hashSize)
  
  let data = reducedImage.asArray()

  var diff: array[hashBits, int]
  var idx: int
  for i in 0 ..< hashSize:
    for j in 0 ..< hashSize:
      idx = i * hashSize + j
      if data[i][j] > data[i][j.succ]:
        diff[idx] = 1
      else:
        diff[idx] = 0

  return ImageHash(diff)


proc pHash*(image: Image, highFreqFactor: Positive = 4): ImageHash =
  ## Perceptual Hash: uses the Discrete Cosine Transform (DCT) to analyze the frequency components of an image,
  ## focusing on low-frequency components for generating a hash that captures visual essence and textures

  let imgSize = hashSize * highFreqFactor
  var reducedImage = image
  reducedImage.filterGreyscale()
  reducedImage = reducedImage.resizedNN(imgSize, imgSize)
  
  let data = reducedImage.asArray()
  let dct = dct2d(data)
  let dctLowFreq = dct.slice(0, hashSize, 0, hashSize)
  let med = median(dctLowFreq.flatten())

  var diff: array[hashBits, int]
  var idx: int
  for i in 0 ..< hashSize:
    for j in 0 ..< hashSize:
      idx = i * hashSize + j
      if dctLowFreq[i][j] > med:
        diff[idx] = 1
      else:
        diff[idx] = 0

  return ImageHash(diff)


proc wHash*(image: Image, scalingFactor: int = 0): ImageHash =
  ## Wavelet Hash: applies wavelet transforms to capture both spatial and frequency details at multiple resolutions,
  ## using wavelet coefficients to generate a hash that balances texture, detail, and layout.
  if not (scalingFactor mod 2 == 0):
    raise newException(ValueError, "scalingFactor must be a power of 2 < smallest image dimension")
  
  var imgSize: int
  if scalingFactor == 0:
    let imageMinDimension = min(image.width, image.height)
    let imageNaturalScale = pow(2.0, floor(log2(float(imageMinDimension))))
    imgSize = max(imageNaturalScale.int, hashSize)
  else:
    imgSize = scalingFactor
  
  let llMaxLevel = int(log2(imgSize.float))
  let level = int(log2(hashSize.float))
  let dwtLevel = llMaxLevel - level

  var reducedImage = image
  reducedImage.filterGreyscale()
  reducedImage = reducedImage.resizedNN(imgSize, imgSize)
  
  let data = reducedImage.asArray()
  let waveletTransform = haarWaveletTransform(data, dwtLevel)
  let dwtLow = waveletTransform.LL
  let med = median(dwtLow.flatten())

  var diff: array[hashBits, int]
  var idx: int
  for i in 0 ..< hashSize:
    for j in 0 ..< hashSize:
      idx = i * hashSize + j
      if dwtLow[i][j] > med:
        diff[idx] = 1
      else:
        diff[idx] = 0

  return ImageHash(diff)

#################################################################################
