import matlab.engine

eng = matlab.engine.start_matlab()
t = eng.isprime(37)
c = eng.get_W(0.05,float(2))
print(c)
print(t)
