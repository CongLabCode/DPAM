import sys

fp = open(sys.argv[1] + '_check','r')
check1 = 0
check2 = 0
check3 = 0
check4 = 0
for line in fp:
    words = line.split()
    check1 += int(words[1])
    check2 += int(words[2])
    check3 += int(words[3])
    check4 += int(words[4])
fp.close()

print (sys.argv[1] + '\t' + str(check1) + '\t' + str(check2) + '\t' + str(check3) + '\t' + str(check4))
