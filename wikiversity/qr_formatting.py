def qr_check_format(fmt: int, debug: int = 1) -> int:  
    """
    This function performs euclidean division between fmt (the dividend) and g (the divisor)
    where g is the generator polynomial for qr format codes. 
    It returns the remainder of the operation.
    It can also be seen as a modulo operation: fmt mod(g).
    It only works for (15, 10) BCH codes, meaning 15 bit block length with 10 redundant bits.
    It optionally prints debug information to the terminal 
    depending on the value of the 'debug' argument.

    Args:
        fmt (int): 15 bit binary number.
        debug (int, optional): A flag that controls whether debug information is printed. Default is 1 (enabled).

    Returns:
        int: Remainder of euclidean division between fmt and 0x537.
    """

    g = 0x537 # = 0b10100110111, generator of the code

    for i in range(4,-1,-1):
        if debug:
            print(f"Format input :\t\t {fmt:015b}")
        # Check for fmt > (2^14, 2^13, 2^12, 2^11, 2^10)
        # if fmt is bigger then XOR with i left-shifted generator
        # ensuring result is AT LEAST 1 bit shorter
        # because you line up the first set bit of the generator 
        # with the first set bit of fmt.
        if fmt & (1 << (i+10)): # if there is a 1 in the 10+i+1 bit [or fmt & 1000...00 (10+i zeros)]
            """
            It might not seem intuitive at first glance but this is simply long division.
            Think of this example: 423 / 2
            1. Subtract 2 100 times: 423 - 200 = 223 (implicit check 2000 > 423 > 200, you can subtract 200 but not 2000)
            2. Subtract 2 100 times: 223 - 200 =  23 (implicit check 2000 > 223 > 200, you can subtract 200 but not 2000)
            3. Subtract 2  10 times:  23 -  20 =   3 (implicit check  200 >  23 >  20, you can subtract  20 but not  200)
            4. Subtract 2   1  time:   3 -   2 =   1 (implicit check   20 >   3 >   2, you can subtract   2 but not   20)
            423 / 2 = 100 + 100 + 10 + 1 remainder 1 = 211 remainder 1
            
            In binary, subtraction is addition since 1 - 1 = 0 = 1 + 1.
            So fmt ^ (g << i) subtracts g 2^i times from fmt. (^ is XOR which is bit-wise addition, with no carry)
            You can think of fmt ^ g as fmt subtracting g 1 time simply because 
            (fmt ^ g) ^ g = fmt, so adding g back negates subtracting g.
            """
            if debug:
                print(f"Bit XOR      :\t\t {(g << i):015b}")
            fmt ^= g << i # subract g 2^i times, and set fmt to that new value
            if debug:
                print(f"Result       :\t\t {fmt:015b}")
    return fmt

def hamming_weight(x):
   weight = 0
   while x > 0:
      weight += x & 1
      x >>= 1
   return weight

def qr_decode_format(fmt):
   best_fmt = -1
   best_dist = 15
   for test_fmt in range(0,32):
      test_code = (test_fmt<<10) ^ qr_check_format(test_fmt<<10)
      test_dist = hamming_weight(fmt ^ test_code)
      if test_dist < best_dist:
         best_dist = test_dist
         best_fmt = test_fmt
      elif test_dist == best_dist:
         best_fmt = -1
   return best_fmt


print("Example format checked: " + bin(qr_check_format(0b000111101011001, debug=0)))

print("Encoding 5 bit formatting information:")
format=0b00011
print(f"=>Format: {format:05b}")
print("")
remainder_format = qr_check_format(format<<10)
print("")
print(f"=>Check format:\t\t {remainder_format:015b}")
encode_format=(format<<10)+remainder_format
qr_check_format(encode_format)
print("")
print(f"=>Final valid format:\t {encode_format:015b}")

