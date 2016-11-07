pi=3.141597623

def gcd(a, b):		#Greatest Common Divisor
    while b:
        a, b = b, a%b
    return a
    
def factorial(n):
    result = 1
    for i in range(1, n+1):
        result *= i
    return result
    
if __name__ == '__main__':
    print("as main program")
else:
    print("as module")