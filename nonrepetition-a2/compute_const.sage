
# Bound 2^omega(N)
def c1(p):
    return  ( 2 / p^(1/4) ).n()

for p in prime_range(40):
    print(p, c1(p))

print(product(c1(p) for p in prime_range(16)))
print(product(prime_range(16)))
# 2^omega(N) <= 4.862 N^(1/4)


print('#'*30)


# Bound sigma0(N)
def c2(p):
    return max([((r+1)/(p^(r/4))).n() for r in range(1,100)])

for p in prime_range(30):
    print(p, c2(p))

print(product(c2(p) for p in prime_range(16)))
print(product(prime_range(16)))
# sigma0(N) <= 8.447 N^(1/4)


