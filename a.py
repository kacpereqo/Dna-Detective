def sex(x, k): return "".join([chr((((ord(n) % 97) + k) % 26) + 97) for n in x])


print(sex("sex", 3))
