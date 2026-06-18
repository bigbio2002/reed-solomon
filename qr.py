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
