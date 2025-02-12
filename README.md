# eccodecs

Poject for compiling and archiving implementations of codecs for different linear codes (mostly Reed-Solomon, BCH, and Reed-Muller linear codes). 

If I'm not the author of the code, credit will be given in the README.md file inside each folder.

In the future: Website for exploring linear codes (eccodecs.com).
- ECCodecs (Error Correction Codecs)
- Have a textbox where you can paste UTF-8/ASCII Code 7/8 text or JPEG/PNG Images (copy to base64 like in google colab) or directly bits and apply an ECC encoding to them. Also an option to decode. High error rate resistant codecs are preferable. But also have simple CRC, Hamming.
- Compare speed of different codecs both with benchmarks (dubious results) and strace of the encoding/decoding algorithms.
- Everything local
- Sanitize input if not local?
