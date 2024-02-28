Based on [imagehash](https://github.com/JohannesBuchner/imagehash). Provides Average, Difference, Perceptual, and Wavelet hash functions. Comes packaged with [imageman](https://github.com/SolitudeSF/imageman) for basic image support.

```nim
let imgx = loadImage[ColorRGBU]("/path/to/imagex.jpeg")
let imgy = loadImage[ColorRGBU]("/path/to/imagey.png")

let hx = pHash(imgx)
let hy = pHash(imgy)

echo hx == hy
echo hx - hy # hamming distance

echo $hx
echo $hy
```
