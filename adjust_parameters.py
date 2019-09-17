import os
for filename in os.listdir('family'):
    f = open('family/{}'.format(filename))
    f2 = open('family2/{}'.format(filename), 'w')
    for i, line in enumerate(f):
        if i == 3:
            f2.write("20 28 1 10 70 0.0 15\n")
        elif i == 4:
            f2.write("400 300 0 0\n")
        else:
            f2.write(line)
    f2.close()
       