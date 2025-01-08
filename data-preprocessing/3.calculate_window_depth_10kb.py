import sys,pickle

with open('Human_hg38_chrom_ends.pkl','rb') as b:
    E = pickle.load(b)

infile = sys.argv[1]
outfile = sys.argv[2]

window = 10000
COUNT={}

data = open(infile)
for lines in data:
    lines = lines.rstrip()
    content = lines.split('\t')
    chrom = content[0]
    rbound = E[chrom]
    pos = int(content[1])
    count = int(content[2])
    left = ((pos-1) // window)*window+1
    right = ((pos-1) // window + 1)*window
    if rbound > right and rbound - right < window:
        right = rbound
    elif rbound < right:
        left -= window
        right = rbound
    key = chrom+':'+str(left)+'-'+str(right)
    if key not in COUNT:
        COUNT[key]=0
    COUNT[key]+=count
data.close()

pipeline = open(outfile,'w')
pipeline.write('chrom\tstart\tend\twindow_size\tdepth\n')
for keys in COUNT:
    chrom = keys.split(':')[0]
    start = (keys.split(':')[1]).split('-')[0]
    end = keys.split('-')[1]
    size = str(int(end) - int(start) + 1)
    depth = str(COUNT[keys])
    pipeline.write(chrom+'\t'+start+'\t'+end+'\t'+size+'\t'+depth+'\n')
pipeline.close()