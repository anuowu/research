max_press = 1.0E6
pres_start = 1.0E-8

pres = []
#pres.pop(pres_start)
deltx = (max_press - pres_start)/120.0
print(deltx)
for i in range(120):
    pres.append(pres_start + i*deltx)
    #pres.append(1)

c = 1
picard_factor = 0.22
IFDEN = 0
ITERA = 0

f = open('input_pressure.dat', 'w')
f.write("{} {} {} {} {} {}\n".format(repr('nb').ljust(4), repr('pres(pa)').ljust(20), repr('c').ljust(4), repr('picard_factor').ljust(4), repr('IFDEN').ljust(4), repr('ITERA').ljust(4) ))
for i in range(120):
    f.write("{} {} {} {} {} {}\n".format(repr(i+1).ljust(4), repr(pres[i]).ljust(20), repr(c).ljust(4), repr(picard_factor).ljust(8), repr(IFDEN).ljust(4), repr(ITERA).ljust(4) ))
f.close()
