# calculate several basic stats in terms of subgenome
import sys

def calcComposition(file):
    Nsyl, Ntom, Unassigned, Ambiguous = 0,0,0,0
    with open(file) as f:
        line = f.readline()
        while line:
            _, start, end, label = line.strip().split('\t')
            span = int(end) - int(start)
            if (span <= 0):
                print(f"span smaller than 1 with start={start} and end={end} ")

            if label == 'Nsyl':
                Nsyl += span
            elif label == 'Ntom':
                Ntom += span
            elif label == 'unassigned':
                Unassigned += span
            else:
                Ambiguous += span
            line = f.readline()
    total = Nsyl + Ntom + Unassigned + Ambiguous
    print("Nsyl: {:.4f}".format(Nsyl/total))
    print("Ntom: {:.4f}".format(Ntom/total))
    print("Unassigned: {:4f}".format(Unassigned/total))
    print("Ambiguous: {:4f}".format(Ambiguous/total))

# calculate 
calcComposition(sys.argv[1])
