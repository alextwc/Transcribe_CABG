import mymath

print('---Hello Python 3.5.2---')
print('pi is ' + str(mymath.pi))
print('gcd of 24 and 16 is ' + str(mymath.gcd(24,16)))
print('factoirial of 6 is ' + str(mymath.factorial(6)))
print('---Bye Python 3.5.2---')

if __name__ == '__main__':
    print("as main program")
else:
    print("as module")