import sys
# this script tests whether the singlebase assignment file and the block assignment file
# are consistent. Used to test the correctness of assignBaseV2.py.

# python testAssign.py <block assignment file> <singlebase assignment file>

block_file = open(sys.argv[1],'r')
single_file = open(sys.argv[2],'r')
block_line = block_file.readline()
#single_line = single_file.readline()

while block_line:
    content_block = block_line.strip().split('\t')
    start, end = int(content_block[1]),int(content_block[2])
    assign_block = content_block[3]
    for i in list(range(start,end)):
        content_single = single_file.readline().strip().split('\t')
        assign_single = content_single[2]
        if assign_single != assign_block:
            print("Inconsistent Assignment at base {} on chromosome {}".format(int(content_single[1]),content_single[0]))
            print("Block Assignment: {}".format(assign_block))
            print("Singlebase Assignment: {}".format(assign_single))

    block_line = block_file.readline()

block_file.close()
single_file.close()
print("Checking is done!")
